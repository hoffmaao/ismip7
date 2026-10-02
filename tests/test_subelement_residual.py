r"""The sub-element grounding scheme reaches ISMIP7 through
icepack_tools.momentum.dual_residual (icepack2_tools.subelement).  Before an
inversion is run on it: the shared residual is ours on the blocks the scheme
does not touch; on fully grounded ice the sub-element friction is the
ordinary one; on a partly grounded cell the basal force scales with the
grounded fraction; and the quadrature follows the geometry the forward would
hand it.

Serial, a 30 km x 10 km slab with a grounding line, no data files.
"""
# dual_friction first: it pulls icepack2 -> irksome, which must be imported
# before any UFL form is assembled.
from icepack2_tools.dual_friction import build_rc_residual, grounded_mask  # noqa: E402

import numpy as np                                              # noqa: E402
import pytest                                                   # noqa: E402

from firedrake import (                                         # noqa: E402
    Constant, FiniteElement, Function, FunctionSpace, RectangleMesh,
    SpatialCoordinate, TensorFunctionSpace, VectorFunctionSpace, as_vector,
    assemble, conditional,
)

from icepack2_tools.subelement import (                         # noqa: E402
    build_subelement_residual, ice_indicator, subelement_from_geometry,
)

RHO_I, RHO_W = 917.0, 1024.0
X_GL = 20e3


def _slab(gl=True):
    mesh = RectangleMesh(30, 10, 30e3, 10e3)
    Q = FunctionSpace(mesh, "CG", 1)
    V = VectorFunctionSpace(mesh, "CG", 1)
    dg0 = FiniteElement("DG", "triangle", 0)
    Z = V * TensorFunctionSpace(mesh, dg0, symmetry=True) * VectorFunctionSpace(mesh, dg0)
    Q0 = FunctionSpace(mesh, "DG", 0)
    x, y = SpatialCoordinate(mesh)
    if gl:
        # grounded for x < X_GL through a thinning wedge over a deepening bed,
        # floating beyond: the grounding line runs across cells
        H = Function(Q0).interpolate(1200.0 - 0.03 * x)
        b = Function(Q0).interpolate(-100.0 - 0.04 * x)
    else:
        H = Function(Q0).interpolate(Constant(1200.0))
        b = Function(Q0).interpolate(Constant(100.0))
    haf = H - Constant(RHO_W / RHO_I) * conditional(b < 0.0, -b, 0.0)
    s = Function(Q0).interpolate(conditional(haf > 0.0, b + H, (1.0 - RHO_I / RHO_W) * H))
    C_w0 = Function(Q, name="C_w0").interpolate(0.05 + 0.02 * y / 10e3)
    theta = Function(Q).interpolate(0.1 * x / 30e3)
    phi = Function(Q).interpolate(-0.05 + 0.1 * y / 10e3)
    A_prior = Function(Q).interpolate(Constant(20.0) + 5.0 * x / 30e3)
    z = Function(Z)
    z.subfunctions[0].interpolate(as_vector([200.0 + 0.02 * x, 5.0 * y / 10e3]))
    z.subfunctions[1].interpolate(Constant(((0.05, 0.01), (0.01, 0.02))))
    z.subfunctions[2].interpolate(as_vector([0.03, 0.005]))
    return dict(mesh=mesh, Q=Q, Q0=Q0, Z=Z, H=H, b=b, s=s, C_w0=C_w0, theta=theta,
                phi=phi, A_prior=A_prior, z=z, n_flow=Constant(3.0), m_slide=Constant(3.0))


COMMON = dict(n_flow_val=3.0, tau_c=Constant(0.1), alpha=Constant(1e-2),
              H_ref=Constant(100.0), alpha_gl=0.5, c_w0_floor=0.0,
              h_visc_floor=10.0, ocean_drag=1e-2, h_ocean=10.0, u_lim=0.0)


def _ours(f, **over):
    kw = dict(COMMON, **over)
    F = build_rc_residual(
        f["z"], f["theta"], f["phi"], H=f["H"], s=f["s"], b=f["b"], C_w0=f["C_w0"],
        A4_base=f["A_prior"], n_flow=f["n_flow"], m_slide=f["m_slide"],
        fric_law="budd", N_ref=None, nhat_floor=0.0, nhat_cap=3.0, k_lim=0.0, **kw)
    return _blocks(f, F)


def _shared(f, subelement=None, exact_front=False, **over):
    kw = dict(COMMON, **over)
    F = build_subelement_residual(
        f["z"], f["theta"], f["phi"], H=f["H"], s=f["s"], b=f["b"], C_w0=f["C_w0"],
        A4_base=f["A_prior"], n_flow=f["n_flow"], m_slide=f["m_slide"],
        subelement=subelement, fric_law="budd", nhat_cap=3.0, k_lim=0.0,
        exact_front=exact_front, **kw)
    return _blocks(f, F)


def _blocks(f, F):
    r = assemble(F)
    return [np.array(sub.dat.data_ro, copy=True) for sub in r.subfunctions]


def _rel(a, b):
    return np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300)


def test_the_shared_residual_is_ours_where_the_scheme_does_not_act():
    f = _slab(gl=True)
    ours = _ours(f)
    shared = _shared(f)
    # velocity block: identical form
    assert _rel(shared[0], ours[0]) < 1e-12
    # membrane block: the shared form regularises |M_dev|^2 by 1e-6 MPa^2
    assert _rel(shared[1], ours[1]) < 1e-4
    # basal-stress block: identical outside the smooth grounding band, where
    # the shared cell-wise Budd gate also multiplies N_hat by He
    He = Function(f["Q0"]).interpolate(grounded_mask(f["H"], f["b"])).dat.data_ro
    outside = (He > 1 - 1e-9) | (He < 1e-9)
    assert outside.sum() > 0.8 * outside.size
    assert _rel(shared[2][outside], ours[2][outside]) < 1e-12


def test_fully_grounded_ice_gets_the_ordinary_friction():
    f = _slab(gl=False)
    sub = subelement_from_geometry(f["mesh"], f["H"], f["b"])
    assert np.all(sub.full.dat.data_ro == 1.0)
    assert np.all(sub.fraction.dat.data_ro == 1.0)
    plain = _shared(f)
    withsub = _shared(f, subelement=sub)
    for k in range(3):
        assert _rel(withsub[k], plain[k]) < 1e-12, k


def test_a_partly_grounded_cell_carries_its_grounded_fraction_of_drag():
    f = _slab(gl=True)
    sub = subelement_from_geometry(f["mesh"], f["H"], f["b"])
    frac = sub.fraction.dat.data_ro
    part = (frac > 0.0) & (frac < 1.0)
    assert part.any()
    # the grounded area is continuous in the geometry: it lies between the
    # cell-wise counts of grounded and of not-floating cells
    area = assemble(Function(f["Q0"]).interpolate(Constant(1.0)) * __import__("firedrake").dx)
    n_cells = frac.size
    cell_area = area / n_cells
    haf_cell = Function(f["Q0"]).interpolate(
        f["H"] - Constant(RHO_W / RHO_I) * conditional(f["b"] < 0.0, -f["b"], 0.0)).dat.data_ro
    assert (haf_cell > 0).sum() * cell_area * 0.9 < frac.sum() * cell_area < (haf_cell > 0).sum() * cell_area * 1.1
    # on a partly grounded cell the basal-stress closure carries the grounded
    # fraction of the Weertman drag, so the residual's tau block there is
    # the whole-cell drag scaled by the fraction (velocity nearly uniform)
    full_law = _shared(f)                          # cell-wise gate, N_hat 1 grounded
    withsub = _shared(f, subelement=sub)
    tau_block_full = full_law[2].reshape(-1, 2)
    tau_block_sub = withsub[2].reshape(-1, 2)
    z_tau = f["z"].subfunctions[2].dat.data_ro
    # the closure residual is tau + drag: subtract tau to compare the drags
    Mtau = assemble(__import__("firedrake").TestFunction(f["Q0"]) * __import__("firedrake").dx).dat.data_ro
    drag_full = tau_block_full - z_tau * Mtau[:, None]
    drag_sub = tau_block_sub - z_tau * Mtau[:, None]
    # fully grounded cells outside the cell-wise gate's smooth band, where
    # that gate is exactly 1: the two closures coincide
    He = Function(f["Q0"]).interpolate(grounded_mask(f["H"], f["b"])).dat.data_ro
    grounded_full = (frac == 1.0) & (He > 1.0 - 1e-9)
    assert grounded_full.sum() > 0
    assert np.allclose(drag_sub[grounded_full], drag_full[grounded_full], rtol=1e-10, atol=0.0)
    ratio = np.linalg.norm(drag_sub[part], axis=1) / np.maximum(
        np.linalg.norm(drag_full[part], axis=1), 1e-300)
    # where the cell-wise gate counts the cell as grounded, the sub-element
    # drag is the fraction of it (to the in-cell variation of the drag)
    gated_on = haf_cell[part] > 0
    assert np.allclose(ratio[gated_on], frac[part][gated_on], rtol=0.05, atol=0.0)


def test_ice_free_cells_get_no_grounded_part():
    f = _slab(gl=False)
    ice = ice_indicator(f["H"], 1.0)
    ice.dat.data[:10] = 0.0
    sub = subelement_from_geometry(f["mesh"], f["H"], f["b"], ice=ice)
    assert np.all(sub.fraction.dat.data_ro[:10] == 0.0)
    assert np.all(sub.fraction.dat.data_ro[10:] == 1.0)
