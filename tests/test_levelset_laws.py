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
#: an interior front. Cells are 1/8 wide, and the imposed speed is small
#: against that, so one step retreats the front by a fraction of a cell.
N = 8
H_ICE = 300.0
U_X = 1.0e-3
HMIN = 1.0


def _state(law, **kwargs):
    mesh = UnitSquareMesh(N, N)
    Q0 = FunctionSpace(mesh, "DG", 0)
    x, _ = SpatialCoordinate(mesh)
    h = Function(Q0).interpolate(H_ICE * (x < 0.5))
    # a bed deep enough that every ice cell floats, so the floating threshold
    # is the one the laws use
    b = Function(Q0).interpolate(Constant(-2000.0))
    u = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        as_vector((Constant(U_X), Constant(0.0))))
    ls = LevelSet(mesh, h, law=law, h_min=HMIN, **kwargs)
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


def test_the_signed_distance_is_negative_in_ice():
    _, _, h, _, _, ls = _state("none")
    phi = ls.phi.dat.data_ro
    ice = h.dat.data_ro > HMIN
    assert np.all(phi[ice] < 0.0)
    assert np.all(phi[~ice] > 0.0)


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


def test_vonmises_retreats_the_front_and_sheds_sub_cell_mass():
    mesh, Q0, h, b, u, ls = _state("vonmises")
    A, n = Constant(20.0), 3.0
    rate = ls.advance(0.1, u, h, b, A, n)
    assert rate > 0.0, "a stressed floating front must calve"
    beyond, frac = ls.calving_masks()
    assert frac is not None
    front = ls.front_len.dat.data_ro > 0.0
    assert np.all(frac[front] > 0.0), "every front cell sheds something"
    assert np.all(frac <= 1.0), "a cell cannot shed more than it holds"
    assert not beyond[h.dat.data_ro > HMIN].any(), \
        "one small step must not carry the front past a whole ice cell"


def test_the_shed_fraction_is_the_swept_area_over_the_cell():
    r"""``min(1, c dt L / A)``: the fraction a front cell loses is the area the
    front sweeps through it, which is what makes the tally a volume."""
    mesh, Q0, h, b, u, ls = _state("prescribed")
    c = Function(Q0).interpolate(Constant(0.05))
    dt = 0.5
    ls.advance(dt, u, h, b, rate=c)
    _, frac = ls.calving_masks()
    front = ls.front_len.dat.data_ro > 0.0
    expected = np.clip(
        ls.c_cell.dat.data_ro * dt * ls.front_len.dat.data_ro / ls.cell_area,
        0.0, 1.0)
    assert np.allclose(frac[front], expected[front])


def test_retreat_scales_with_the_step():
    r"""Two half steps remove what one full step does, to the order of the
    scheme: a run that halves its timestep must not calve a different amount.
    """
    def shed(dt, steps):
        mesh, Q0, h, b, u, ls = _state("prescribed")
        c = Function(Q0).interpolate(Constant(0.02))
        total = 0.0
        for _ in range(steps):
            ls.advance(dt, u, h, b, rate=c)
            _, frac = ls.calving_masks()
            total += float((frac * ls.cell_area).sum())
        return total

    one = shed(0.4, 1)
    two = shed(0.2, 2)
    assert one > 0.0
    assert two == pytest.approx(one, rel=0.25)


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
    beyond, _ = ls.calving_masks()
    xc = _xc(mesh)
    # the bar is the ORIGINAL front, so the band between them is still open
    assert np.all(beyond == (xc > 0.5))
    assert not beyond[(xc > 0.25) & (xc < 0.5)].any()
