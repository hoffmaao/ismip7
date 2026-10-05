r"""subelement.py - the momentum residual with the sub-element grounding scheme.

ISSM's SEP2 (Seroussi et al. 2014): height above flotation is taken as
linear on each triangle, the grounding line inside it is a straight
segment, and the unchanged basal friction is integrated over the grounded
part only, so the grounding line migrates continuously instead of cell by
cell (the staircase ``GEOMETRY_DISCRETIZATION.md`` names).  The scheme lives
in :mod:`icepack_tools.grounding` (``SubelementGrounding``) and is applied
by :func:`icepack_tools.momentum.dual_residual`; this module is what ISMIP7
adds around it:

* the quadrature is rebuilt from ISMIP7's own geometry
  (:func:`subelement_from_geometry`): height above flotation from the CG1
  lifts of the DG0 thickness and bed, and no grounded part in cells that
  hold no ice;
* the residual is the shared one plus ISMIP7's floor-cell ocean drag, which
  the shared builder does not carry (:func:`build_subelement_residual`).

What the scheme changes in the friction law: the grounded stress must need
nothing but the velocity and the controls, so the Budd law runs with
``N_ref=None`` and no delta floor (``N_hat = 1`` on the grounded part, zero
afloat, the grounding line exact); the smooth ``He`` band of the cell-wise
gate plays no part.  The forward therefore has no effective-pressure
feedback under this scheme until the shared builder evaluates ``N`` at the
grounded-part points.  Andrew (26 Sep 2026): "start the process of setting
up the inversions with these grounding zone fixes".

``scheme="sep1"`` is ISSM's default ``SubelementFriction1`` instead: the
cell's ordinary quadrature with the drag scaled by its grounded fraction
(``alpha2 = phi * alpha2`` in ISSM's ``CreateKMatrixSSAFriction``).  The
drag is the same Weertman stress, ``exp(theta)`` unscaled by the smooth
``He`` band as under SEP2, zero on floating cells, but it varies
with the grounded fraction alone, without SEP2's quadrature points that
jump across the grounding line inside a cell.  At 32 km SEP2 alone turned
a run with no failed line-search trials into one with twelve (30 Sep 2026).
"""
from firedrake import (Constant, Function, FunctionSpace, TestFunction, dx,
                       exp, inner, max_value, split, sqrt)

from icepack2_tools.fssa import fssa_term
from icepack2_tools.geometry import cg1_lift
from icepack2_tools.grounding import height_above_flotation


def subelement_from_geometry(mesh, H, b, ice=None, sub=None):
    r"""A :class:`icepack_tools.grounding.SubelementGrounding` current for
    the geometry ``(H, b)`` (DG0 or CG1 Functions).

    The flotation field at the vertices is the CG1 lift of the cell fields
    (a convex combination of the cell values around each vertex, so it
    cannot overshoot).  ``ice`` (DG0, 1 in ice cells) removes the grounded
    part of cells holding no ice.  Pass ``sub`` to update an existing object
    in place (the forward, every step) instead of building one."""
    from icepack_tools.grounding import SubelementGrounding
    if sub is None:
        sub = SubelementGrounding(mesh)
    H_cg = H if H.function_space().ufl_element().family() != "Discontinuous Lagrange" else cg1_lift(H)
    b_cg = b if b.function_space().ufl_element().family() != "Discontinuous Lagrange" else cg1_lift(b)
    sub.update(height_above_flotation(H_cg, b_cg), ice=ice)
    return sub


def ice_indicator(H, h_min):
    r"""DG0 indicator of the cells holding ice (``H > h_min``)."""
    Q0 = FunctionSpace(H.function_space().mesh(), "DG", 0)
    ice = Function(Q0, name="ice_cell")
    H0 = Function(Q0).project(H) if H.function_space() != Q0 else H
    ice.dat.data[:] = (H0.dat.data_ro > h_min).astype(float)
    return ice


def ocean_drag_closure(z, H, ocean_drag, h_ocean, drag_mask=None, u_min=1.0):
    r"""ISMIP7's floor-cell drag as an addition to the basal-stress closure:
    ``ocean_drag * gate * max(0, 1 - H/h_ocean) * |u|`` along ``u``, exactly
    zero for ``H >= h_ocean`` (:func:`icepack2_tools.dual_friction.build_rc_residual`)."""
    if not ocean_drag:
        return 0
    u, _, _ = split(z)
    _, _, sig = split(TestFunction(z.function_space()))
    u_reg = sqrt(inner(u, u) + Constant(u_min) ** 2)
    gate = Constant(1.0) if drag_mask is None else drag_mask
    tau_o = (Constant(ocean_drag) * gate
             * max_value(Constant(0.0), Constant(1.0) - H / Constant(h_ocean)) * u_reg)
    return inner(tau_o * u / u_reg, sig) * dx


SCHEMES = ("sep2", "sep1")


def build_subelement_residual(z, theta, phi, *, H, s, b, C_w0, A4_base, n_flow,
                              n_flow_val, m_slide, tau_c, alpha, H_ref,
                              subelement, scheme="sep2", fric_law="budd", nhat_cap=3.0,
                              alpha_gl=0.0, c_w0_floor=0.0, h_visc_floor=0.0,
                              ocean_drag=0.0, h_ocean=10.0, drag_mask=None,
                              u_lim=0.0, k_lim=0.0, gl_width=10.0,
                              calving_ids=None, exact_front=True, u_min=1.0,
                              fssa_tau=None, u_ref=None):
    r"""The single-layer dual residual with the sub-element grounded friction.

    Same controls and fields as ``build_rc_residual``; the friction law is
    Budd with ``N_ref=None`` and no floor (or Weertman), which is what the
    scheme supports.  ``exact_front`` adds the exact depth-integrated push on
    a calving front inside the mesh (``icepack_tools.momentum.front_cliff_correction``).
    ``scheme`` is ``"sep2"`` (grounded-part quadrature) or ``"sep1"``
    (whole-cell quadrature, drag times the grounded fraction).  ``fssa_tau``
    (a Constant, theta times the step) and ``u_ref`` add the free-surface
    stabilization of :func:`icepack2_tools.fssa.fssa_term`.
    """
    from icepack_tools.momentum import dual_residual
    if scheme not in SCHEMES:
        raise ValueError(f"sub-element scheme must be one of {SCHEMES}, not {scheme!r}")
    if fric_law not in ("budd", "weertman"):
        raise ValueError(
            f"the sub-element grounding scheme supports budd or weertman, not {fric_law!r}")
    # the geometry's mesh: a mixed space's may come back as a MeshSequence
    mesh = H.function_space().mesh()
    if isinstance(theta, Constant) and theta.ufl_domain() is None:
        # the scheme evaluates exp(theta) at points through grad(theta); a
        # domainless Constant (the sqrt control's zero deviation) has none
        theta = Function(FunctionSpace(mesh, "R", 0)).assign(float(theta))
    A_eff = A4_base * exp(phi)
    if scheme == "sep1":
        # Budd with N_hat = 1 on the grounded part is the Weertman stress;
        # the floor applies to the coefficient, never to floating cells.
        # exp(theta) rides in the coefficient, unscaled as under SEP2: the
        # shared Weertman closure would damp theta by the smooth He band
        C = max_value(C_w0, Constant(c_w0_floor)) if c_w0_floor else C_w0
        F = dual_residual(
            z, Constant(0.0), phi, H=H, s=s, b=b, h_layers=[H],
            C_w0=C * subelement.fraction * exp(theta),
            A_layers=[A_eff], n_consts=[n_flow], n_vals=[n_flow_val],
            m_slide=m_slide, mesh=mesh, law="weertman", tau_c=float(tau_c),
            alpha=float(alpha), H_ref=float(H_ref), u_min=u_min,
            N_ref=None, nhat_floor=0.0, nhat_cap=nhat_cap,
            c_w0_floor=0.0, h_visc_floor=h_visc_floor, alpha_gl=alpha_gl,
            u_lim=u_lim, k_lim=k_lim, gl_width=gl_width,
            outflow_ids=tuple(calving_ids) if calving_ids else None,
            subelement=None, exact_front=exact_front,
        )
        return (F + ocean_drag_closure(z, H, ocean_drag, h_ocean, drag_mask, u_min=u_min)
                + fssa_term(z, u_ref, fssa_tau, H, b, gl_width=gl_width))
    F = dual_residual(
        z, theta, phi, H=H, s=s, b=b, h_layers=[H], C_w0=C_w0,
        A_layers=[A_eff], n_consts=[n_flow], n_vals=[n_flow_val],
        m_slide=m_slide, mesh=mesh, law=fric_law, tau_c=float(tau_c),
        alpha=float(alpha), H_ref=float(H_ref), u_min=u_min,
        N_ref=None, nhat_floor=0.0, nhat_cap=nhat_cap,
        c_w0_floor=c_w0_floor, h_visc_floor=h_visc_floor, alpha_gl=alpha_gl,
        u_lim=u_lim, k_lim=k_lim, gl_width=gl_width,
        outflow_ids=tuple(calving_ids) if calving_ids else None,
        subelement=subelement, exact_front=exact_front,
    )
    return (F + ocean_drag_closure(z, H, ocean_drag, h_ocean, drag_mask, u_min=u_min)
            + fssa_term(z, u_ref, fssa_tau, H, b, gl_width=gl_width))
