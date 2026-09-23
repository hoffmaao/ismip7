r"""What each calving law does to the front, and what one step removes.

The level-set unit tests were lost (issue #35). These rebuild the part that
belongs here: the laws as ISMIP7 configures them through
``icepack2_tools.levelset`` (the ``extent`` anchoring default, the thresholds
from ``runconfig``), the masks the transport applies afterwards, the drag gate
the momentum residual reads, and the step-size behaviour of the retreat. The
shared class's own operators are tested in ``icepack_tools``; what is checked
here is the behaviour this repository depends on.

Serial, one small unit mesh, no data files and no momentum solve: a velocity
and a thickness are imposed, so every number below follows from the law and
the geometry alone.
"""

import numpy as np
import pytest

from firedrake import (Constant, Function, FunctionSpace, SpatialCoordinate,
                       UnitSquareMesh, VectorFunctionSpace, as_vector)

from icepack2_tools.levelset import LevelSet, initial_distance

#: Ice on the left half of a unit domain, open water beyond, so the extent has
#: an interior front. Cells are 1/8 wide; on the default diagonal every front
#: cell is a triangle with one 1/8 facet on x = 0.5 and area 1/128.
N = 8
H_ICE = 300.0
HMIN = 1.0
FRONT_LEN = 1.0 / N
CELL_AREA = 1.0 / (2 * N * N)
#: An extending flow u = (STRAIN x, 0), so e_xx = STRAIN and the von Mises
#: stress is set by the strain rate rather than by a regulariser.
STRAIN = 0.1
#: A bed deep enough that every ice cell floats, and one shallow enough that
#: every ice cell is grounded.
BED_FLOATING = -2000.0
BED_GROUNDED = -100.0
SIGMA_MAX_GROUNDED = 1.0
SIGMA_MAX_FLOATING = 0.15


def _state(law, bed=BED_FLOATING, **kwargs):
    mesh = UnitSquareMesh(N, N)
    Q0 = FunctionSpace(mesh, "DG", 0)
    x, _ = SpatialCoordinate(mesh)
    h = Function(Q0).interpolate(H_ICE * (x < 0.5))
    b = Function(Q0).interpolate(Constant(bed))
    u = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        as_vector((STRAIN * x, Constant(0.0))))
    ls = LevelSet(mesh, h, law=law, h_min=HMIN,
                  sigma_max_grounded=SIGMA_MAX_GROUNDED,
                  sigma_max_floating=SIGMA_MAX_FLOATING, **kwargs)
    return mesh, Q0, h, b, u, ls


def _xc(mesh):
    return Function(VectorFunctionSpace(mesh, "DG", 0)).interpolate(
        SpatialCoordinate(mesh)).dat.data_ro[:, 0]


def test_the_extent_anchor_is_the_default():
    r"""ISMIP7 anchors on the ice extent rather than on a fixed contour,
    because on a sheet with melt and a buffered ocean the thickness is the
    authority on where ice is."""
    _, _, _, _, _, ls = _state("none")
    assert ls.anchor == "extent"


def test_a_front_cell_is_one_with_an_ice_free_neighbour():
    mesh, _, h, _, _, ls = _state("none")
    front = ls.front_len.dat.data_ro > 0.0
    xc = _xc(mesh)
    # the front runs at x = 0.5, so every front cell sits in the column
    # immediately inside it
    assert front.any()
    assert np.all(xc[front] < 0.5)
    assert np.all(xc[front] > 0.5 - 2.0 / N)


def test_fixed_holds_the_front_and_removes_nothing_from_it():
    r"""`fixed` is the control configuration: the front cannot move, so the
    rate is zero and only the cells beyond the t=0 extent are flagged."""
    mesh, Q0, h, b, u, ls = _state("fixed")
    rate = ls.advance(0.1, u, h, b)
    assert rate == 0.0
    beyond, frac = ls.calving_masks()
    assert frac is None, "a held front sheds no sub-cell mass"
    xc = _xc(mesh)
    assert np.all(beyond == (xc > 0.5))


@pytest.mark.parametrize("bed, sigma_max", [
    (BED_FLOATING, SIGMA_MAX_FLOATING),
    (BED_GROUNDED, SIGMA_MAX_GROUNDED),
])
def test_vonmises_rate_is_speed_times_stress_over_threshold(bed, sigma_max):
    r"""``c = |u| sqrt(3) A^(-1/n) e~^(1/n) / sigma_max``. Under uniaxial
    extension ``e~ = e_xx / sqrt(2)``, and the threshold is the floating one
    over the deep bed and the grounded one over the shallow bed."""
    mesh, Q0, h, b, u, ls = _state("vonmises", bed=bed)
    A, n, dt = 20.0, 3.0, 0.1
    rate = ls.advance(dt, u, h, b, Constant(A), n)
    front = ls.front_len.dat.data_ro > 0.0
    xc = _xc(mesh)
    sigma = np.sqrt(3.0) * A ** (-1.0 / n) * (STRAIN / np.sqrt(2.0)) ** (1.0 / n)
    expected = STRAIN * xc[front] * sigma / sigma_max
    assert np.allclose(ls.c_cell.dat.data_ro[front], expected, rtol=1e-6)
    assert rate == pytest.approx(expected.mean(), rel=1e-6)
    beyond, frac = ls.calving_masks()
    assert np.allclose(frac[front], expected * dt * FRONT_LEN / CELL_AREA,
                       rtol=1e-6)
    assert np.all(frac[~front] == 0.0)
    assert not beyond[h.dat.data_ro > HMIN].any(), \
        "one small step must not carry the front past a whole ice cell"


def test_the_shed_fraction_is_the_swept_area_over_the_cell():
    r"""``min(1, c dt L / A)``: the fraction a front cell loses is the area the
    front sweeps through it, which is what makes the tally a volume. Every
    front cell here has L = 1/8 and A = 1/128, so ``0.05 * 0.5 * 16 = 0.4``."""
    mesh, Q0, h, b, u, ls = _state("prescribed")
    c = Function(Q0).interpolate(Constant(0.05))
    ls.advance(0.5, u, h, b, rate=c)
    beyond, frac = ls.calving_masks()
    front = ls.front_len.dat.data_ro > 0.0
    assert front.sum() == N
    assert np.allclose(ls.front_len.dat.data_ro[front], FRONT_LEN)
    assert np.allclose(ls.cell_area[front], CELL_AREA)
    assert np.allclose(ls.c_cell.dat.data_ro[front], 0.05)
    assert np.allclose(frac[front], 0.4)
    assert np.all(frac[~front] == 0.0)
    assert not beyond.any()


def test_retreat_scales_with_the_step():
    r"""Two half steps remove what one full step does, to the order of the
    scheme: a run that halves its timestep must not calve a different amount.

    The masks are applied the way the transport applies them: cells beyond
    the front are emptied, then each front cell loses ``frac`` of what it
    still holds. One step of fraction ``f`` sheds ``f H``; two of ``f / 2``
    shed ``(f / 2)(2 - f / 2) H``, which differs from it at ``O(f^2)``.
    """
    c0, dt = 0.02, 0.4

    def shed(steps):
        mesh, Q0, h, b, u, ls = _state("prescribed")
        c = Function(Q0).interpolate(Constant(c0))
        total = 0.0
        for _ in range(steps):
            ls.advance(dt / steps, u, h, b, rate=c)
            beyond, frac = ls.calving_masks()
            data = h.dat.data
            total += float((data[beyond] * ls.cell_area[beyond]).sum())
            data[beyond] = 0.0
            removed = data * frac
            total += float((removed * ls.cell_area).sum())
            data -= removed
        return total

    f = c0 * dt * FRONT_LEN / CELL_AREA
    one, two = shed(1), shed(2)
    assert one == pytest.approx(N * H_ICE * CELL_AREA * f, rel=1e-9)
    assert two == pytest.approx(
        N * H_ICE * CELL_AREA * (f / 2) * (2 - f / 2), rel=1e-9)
    assert two == pytest.approx(one, rel=0.05)


def test_a_zero_rate_law_leaves_the_front_where_it_was():
    mesh, Q0, h, b, u, ls = _state("prescribed")
    before = ls.phi.dat.data_ro.copy()
    ls.advance(0.1, u, h, b, rate=Function(Q0).interpolate(Constant(0.0)))
    assert np.allclose(ls.phi.dat.data_ro, before, atol=1e-9)


def test_the_drag_gate_is_off_at_the_front_and_on_beyond_it():
    r"""The floor-cell ocean drag must not reach the cells the front runs
    through, or it damps the ice rather than the water. The gate is what the
    momentum residual holds, so it is checked against the distance itself."""
    mesh, Q0, h, b, u, ls = _state("vonmises")
    ls.advance(0.1, u, h, b, Constant(20.0), 3.0)
    gate = ls.drag_mask.dat.data_ro
    phi = ls.phi.dat.data_ro
    diam = ls.cell_diam.dat.data_ro
    assert np.all(gate[phi <= diam] == 0.0)
    assert np.all(gate[phi > diam] == 1.0)


def test_an_unknown_law_is_refused():
    with pytest.raises(ValueError, match="front law"):
        _state("nonesuch")


def test_a_law_without_its_input_is_refused():
    r"""Both refusals matter: a missing rate or a missing rheology would
    otherwise be read as a zero rate, and the front would silently hold."""
    _, _, h, b, u, ls = _state("prescribed")
    with pytest.raises(ValueError, match="rate"):
        ls.advance(0.1, u, h, b)
    _, _, h, b, u, vm = _state("vonmises")
    with pytest.raises(ValueError, match="vonmises"):
        vm.advance(0.1, u, h, b)


def test_fixed_anchors_on_the_t0_extent_not_the_restarted_one():
    r"""A resumed run must not re-freeze the front where it had already
    retreated to: that would permanently bar cells a continuous run of the
    same length would keep. ``phi_init`` from the t=0 thickness is what the
    forward passes for exactly this reason.
    """
    mesh = UnitSquareMesh(N, N)
    Q0 = FunctionSpace(mesh, "DG", 0)
    x, _ = SpatialCoordinate(mesh)
    h0 = Function(Q0).interpolate(H_ICE * (x < 0.5))       # the t=0 extent
    h_now = Function(Q0).interpolate(H_ICE * (x < 0.25))   # already retreated
    phi0 = initial_distance(mesh, h0, h_min=HMIN)
    ls = LevelSet(mesh, h_now, law="fixed", h_min=HMIN, phi_init=phi0)
    # the forward advances before it asks for the masks, and the advance
    # re-solves phi from the current extent, so the bar must survive it
    b = Function(Q0).interpolate(Constant(BED_FLOATING))
    u = Function(VectorFunctionSpace(mesh, "CG", 1))
    assert ls.advance(0.1, u, h_now, b) == 0.0
    beyond, _ = ls.calving_masks()
    xc = _xc(mesh)
    # the bar is the ORIGINAL front, so the band between them is still open
    assert np.all(beyond == (xc > 0.5))
    assert not beyond[(xc > 0.25) & (xc < 0.5)].any()


# The two laws of ``icepack2_tools.calving_laws``, driven through the level set
# the way the forward wires them: ``law="prescribed"`` with the UFL rate.

#: One mesh, several margins. Ice fills x < 0.5 except for a nunatak hole; the
#: front at x = 0.5 is cut into bands by y, each a different margin, and every
#: band edge falls on a cell row so no front cell straddles two bands.
HC = 150.0
U_FRONT = 0.01
_BANDS = {
    # name: (y range, thickness, bed)
    "marine_at_hc": ((0.0, 0.25), HC, -500.0),
    "marine_thin": ((0.25, 0.5), HC / 2, -500.0),
    # 100 m on a 50 m deep bed has 44 m above flotation, so it is grounded
    "grounded_cliff": ((0.5, 0.75), 100.0, -50.0),
    "land": ((0.75, 1.0), 50.0, 200.0),
}
#: The nunatak: one ice-free cell square well inside the ice, with the bed
#: above sea level under it and under the ice ringing it.
_NUNATAK_HOLE = (0.125, 0.25, 0.375, 0.5)
_NUNATAK_BED = (0.0, 0.375, 0.25, 0.625)


def _inside(box, x, y):
    x0, x1, y0, y1 = box
    return (x > x0) & (x < x1) & (y > y0) & (y < y1)


def _margins():
    mesh = UnitSquareMesh(N, N)
    Q0 = FunctionSpace(mesh, "DG", 0)
    xy = Function(VectorFunctionSpace(mesh, "DG", 0)).interpolate(
        SpatialCoordinate(mesh)).dat.data_ro
    x, y = xy[:, 0], xy[:, 1]
    band = np.empty(len(x), dtype=object)
    h = Function(Q0)
    b = Function(Q0)
    for name, ((y0, y1), thickness, bed) in _BANDS.items():
        rows = (y > y0) & (y < y1)
        band[rows] = name
        h.dat.data[rows] = np.where(x[rows] < 0.5, thickness, 0.0)
        b.dat.data[rows] = bed
    hole = _inside(_NUNATAK_HOLE, x, y)
    h.dat.data[hole] = 0.0
    b.dat.data[_inside(_NUNATAK_BED, x, y)] = 20.0
    u = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        as_vector((Constant(U_FRONT), Constant(0.0))))
    ls = LevelSet(mesh, h, law="prescribed", h_min=HMIN)
    return mesh, h, b, u, ls, x, y, band


def test_the_thickness_law_calves_marine_fronts_and_spares_land_ones():
    r"""The bed gate, in one mesh with every kind of margin the level set
    anchors on. At ``H = Hc`` the removal is the ice's own arrival speed, so
    the front holds; thinner marine ice, a grounded cliff on a bed below sea
    level included, is removed faster than it arrives; a land margin and the
    ice ringing a nunatak, both on a bed above sea level, shed nothing, so
    nothing of theirs reaches the calving tally."""
    from icepack2_tools.calving_laws import thickness_calving_rate
    mesh, h, b, u, ls, x, y, band = _margins()
    ls.advance(1.0, u, h, b, rate=thickness_calving_rate(u, h, b, HC))
    beyond, frac = ls.calving_masks()
    c = ls.c_cell.dat.data_ro
    front = ls.front_len.dat.data_ro > 0.0
    ring = front & _inside(_NUNATAK_BED, x, y)
    outer = front & ~ring
    shed = h.dat.data_ro * frac * ls.cell_area
    assert not beyond.any()

    def at(name):
        cells = outer & (band == name)
        assert cells.any(), name
        return cells

    assert c[at("marine_at_hc")] == pytest.approx(U_FRONT, rel=1e-9)
    assert np.all(c[at("marine_thin")] == pytest.approx(1.5 * U_FRONT, rel=1e-9))
    assert np.all(c[at("grounded_cliff")] == pytest.approx(
        (2.0 - 100.0 / HC) * U_FRONT, rel=1e-9))
    for name in ("marine_thin", "grounded_cliff"):
        assert np.all(shed[at(name)] > 0.0), name
    assert ring.sum() >= 4
    for spared in (at("land"), ring):
        assert np.all(c[spared] == 0.0)
        assert np.all(frac[spared] == 0.0)
        assert shed[spared].sum() == 0.0
    # every calving cell is a marine one
    assert np.all(b.dat.data_ro[shed > 0.0] < 0.0)


def test_the_hfb_law_on_the_level_set_follows_the_stress_ratio(monkeypatch):
    r"""The resistive-stress law through the level set with the knobs as
    ``runconfig`` resolves them: a floating front at Buck's unbuttressed
    threshold is removed at the ice speed, one at twice it at twice the speed,
    one far past it at ``ratio_max`` times, and a buttressed or compressive
    front not at all. A tensile strength raises the threshold, so the same
    stress calves more slowly."""
    from firedrake import TensorFunctionSpace, as_matrix
    from icepack2.constants import (gravity as G, ice_density as RHO_I,
                                    water_density as RHO_W)
    from icepack2_tools.calving_laws import hfb_calving_rate
    from icepack2_tools.runconfig import calving_hfb_parameters
    for k in ("SIGMA_MAX", "RHO_C", "HFB_EXPONENT", "HFB_RATIO_MAX"):
        monkeypatch.delenv(f"ISMIP7_CALVING_{k}", raising=False)
    buck = (float(RHO_I) * float(G) * H_ICE * (1.0 - float(RHO_I / RHO_W)) / 2.0)

    def front_rate(stress):
        mesh, Q0, h, b, _, ls = _state("prescribed")
        u = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
            as_vector((Constant(U_FRONT), Constant(0.0))))
        M = Function(TensorFunctionSpace(mesh, "DG", 0, symmetry=True)
                     ).interpolate(as_matrix(((stress, 0.0), (0.0, 0.0))))
        rate = hfb_calving_rate(u, M, h, b, ls.ghat, **calving_hfb_parameters())
        # short enough that even the capped rate stays inside the front cell
        ls.advance(0.1, u, h, b, rate=rate)
        front = ls.front_len.dat.data_ro > 0.0
        _, frac = ls.calving_masks()
        return ls.c_cell.dat.data_ro[front], frac[front]

    for stress, ratio in ((buck, 1.0), (2.0 * buck, 2.0), (100.0 * buck, 5.0),
                          (0.0, 0.0), (-buck, 0.0)):
        c, frac = front_rate(stress)
        assert np.all(c == pytest.approx(ratio * U_FRONT, rel=1e-6, abs=1e-15))
        assert np.all((frac > 0.0) == (ratio > 0.0))
    monkeypatch.setenv("ISMIP7_CALVING_SIGMA_MAX", "0.15")
    c, _ = front_rate(buck)
    assert np.all(c < U_FRONT)
