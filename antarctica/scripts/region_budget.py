#!/usr/bin/env python3
r"""Where the mass goes: the budget split into grounded and floating ice.

The run's own timeseries closes the budget for the ice sheet as a whole
(``SMB - melt - calving - outflux = dM/dt``), which is what a conservation
check needs and exactly the wrong resolution for the question "why is the
volume above flotation rising?".  VAF is grounded ice, and grounded ice does
not see shelf melt or the calving front at all: it gains from surface mass
balance and loses across the grounding line.  So a control that gains VAF is
saying one of two things,

.. code::

    dM_gr/dt = SMB_gr - Q_gl        (grounded)
    dM_fl/dt = SMB_fl + Q_gl - melt - calving - outflux    (floating)

either its grounded accumulation is too large or its grounding-line flux is
too small, and the ~2000-2200 Gt/yr of observed Antarctic discharge (Rignot
et al. 2019; Gardner et al. 2018; IMBIE) is the number that settles it.

This script computes both from a checkpoint pair: the flux across the
grounding line with the model's own DG0 upwind operator (the one the
transport uses, so the number is the flux the run actually applied), the
flux across the calving front the same way, and the SMB integrated over each
region with the run's own RACMO climatology.  ``dM/dt`` per region comes from
differencing the two checkpoints, so the residual of each line is a check on
the diagnosis rather than an assumption.

    python antarctica/scripts/region_budget.py A.h5 [B.h5] [--csv run_timeseries.csv]

With one checkpoint the fluxes and the SMB are reported at that state; with
two, the mass tendencies are differenced over the interval as well.  The
environment knobs are the run's (``ISMIP7_CLIM_START/END``, ``ISMIP7_LC``
and so on), so the SMB is the field the run was forced with.

Observed comparisons printed alongside:

===========================  ==================  ==========================
quantity                     observed            source
===========================  ==================  ==========================
grounded SMB                 2100 +/- 100 Gt/yr  van Wessem et al. 2018 RACMO2
ice-sheet SMB (incl shelves) 2500 +/- 150 Gt/yr  van Wessem et al. 2018
grounding-line discharge     2050 +/- 100 Gt/yr  Rignot et al. 2019
calving flux                 1300 +/- 140 Gt/yr  Rignot et al. 2013
shelf basal melt             1100-1300 Gt/yr     Rignot 2013 / Adusumilli 2020
===========================  ==================  ==========================
"""

import argparse
import csv
import os
import sys

import numpy as np
import firedrake as fd
from firedrake import (
    Constant, Function, assemble, conditional, dS, ds, dx, gt, lt, max_value,
)

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(_ROOT))

RHO_I = 917.0
RHO_W = 1024.0
RHO_GT = RHO_I / 1e12                     # m^3 of ice -> Gt
MM_SLE = RHO_I / 1e12 / 362.5             # m^3 of ice above flotation -> mm SLE

OBSERVED = {
    "grounded SMB": (2100.0, 100.0, "van Wessem 2018 RACMO2"),
    "total SMB": (2500.0, 150.0, "van Wessem 2018"),
    "grounding-line discharge": (2050.0, 100.0, "Rignot 2019"),
    "calving flux": (1300.0, 140.0, "Rignot 2013"),
    "shelf basal melt": (1200.0, 200.0, "Rignot 2013 / Adusumilli 2020"),
}


def load_state(path):
    r"""Mesh, DG0 geometry and velocity from an ISMIP7 checkpoint."""
    with fd.CheckpointFile(path, "r") as chk:
        mesh = chk.load_mesh()
        out = {"mesh": mesh, "t_yr": (float(chk.get_attr("/", "t_yr"))
                                      if chk.has_attr("/", "t_yr") else None)}
        for name in ("thickness", "surface", "bed", "velocity"):
            out[name] = chk.load_function(mesh, name=name)
        # A MAP also carries the observations the inversion fitted.
        try:
            out["velocity_obs"] = chk.load_function(mesh, name="velocity_obs")
        except Exception:
            out["velocity_obs"] = None
    return out


def regions(state):
    r"""Grounded / floating / ice indicators as DG0 functions, and the height
    above flotation.

    The space is DG0 whatever the checkpoint's geometry space is
    (``ISMIP7_GEOMETRY_SPACE=cg1`` stores h/s/b in CG1). The facet integrals in
    :func:`crossing_flux` gate on ``source('+') * sink('-')``, which is
    identically zero for a continuous indicator, so a CG1 indicator would
    report every flux as exactly zero without an error."""
    mesh = state["mesh"]
    Q0 = fd.FunctionSpace(mesh, "DG", 0)
    h = Function(Q0, name="h_dg0").project(state["thickness"])
    b = Function(Q0, name="b_dg0").project(state["bed"])
    haf = Function(Q0, name="haf").interpolate(
        h - Constant(RHO_W / RHO_I) * max_value(-b, Constant(0.0)))
    ice = Function(Q0, name="ice")
    ice.dat.data[:] = np.where(h.dat.data_ro > 1.0, 1.0, 0.0)
    gr = Function(Q0, name="grounded")
    gr.dat.data[:] = np.where((haf.dat.data_ro > 0.0)
                              & (h.dat.data_ro > 1.0), 1.0, 0.0)
    fl = Function(Q0, name="floating")
    fl.dat.data[:] = ice.dat.data_ro - gr.dat.data_ro
    return {"haf": haf, "ice": ice, "grounded": gr, "floating": fl}


def crossing_flux(state, reg, source, sink):
    r"""Upwind volume flux [Gt/yr] across the facets from ``source`` cells to
    ``sink`` cells, with the transport's own operator: ``h u.n`` taken from
    the upwind cell.  ``source``/``sink`` are DG0 indicator Functions."""
    mesh = state["mesh"]
    h, u = state["thickness"], state["velocity"]
    n = fd.FacetNormal(mesh)
    a, b = reg[source], reg[sink]
    un = fd.dot(u, n)
    # '+' side is the source: flux leaving it through this facet
    f_p = conditional(gt(a('+') * b('-'), 0.0), 1.0, 0.0) \
        * max_value(un('+'), 0.0) * h('+')
    f_m = conditional(gt(a('-') * b('+'), 0.0), 1.0, 0.0) \
        * max_value(un('-'), 0.0) * h('-')
    return float(assemble((f_p + f_m) * dS)) * RHO_GT


def crossing_flux_binned(state, reg, source, sink, speed_cell, bins,
                         velocity=None):
    r"""``crossing_flux`` split by the speed of the SOURCE cell, so a flux
    deficit can be attributed to the fast outlets or to the slow margins.
    ``speed_cell`` is a DG0 speed field, ``bins`` its edges [m/yr]."""
    mesh = state["mesh"]
    h = state["thickness"]
    u = velocity if velocity is not None else state["velocity"]
    n = fd.FacetNormal(mesh)
    a, b = reg[source], reg[sink]
    un = fd.dot(u, n)
    out = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        def band(side):
            return (conditional(gt(speed_cell(side), Constant(lo)), 1.0, 0.0)
                    * conditional(lt(speed_cell(side), Constant(hi)), 1.0, 0.0))
        f_p = conditional(gt(a('+') * b('-'), 0.0), 1.0, 0.0) * band('+') \
            * max_value(un('+'), 0.0) * h('+')
        f_m = conditional(gt(a('-') * b('+'), 0.0), 1.0, 0.0) * band('-') \
            * max_value(un('-'), 0.0) * h('-')
        out.append(float(assemble((f_p + f_m) * dS)) * RHO_GT)
    return out


def boundary_flux(state, reg, source):
    r"""Outflow [Gt/yr] across the mesh boundary from ``source`` cells."""
    mesh = state["mesh"]
    h, u = state["thickness"], state["velocity"]
    n = fd.FacetNormal(mesh)
    un = fd.dot(u, n)
    return float(assemble(reg[source] * max_value(un, 0.0) * h * ds)) * RHO_GT


def integrate(field, reg=None, key=None):
    r"""``Gt/yr`` (or Gt) of a DG0 field over a region."""
    w = field if reg is None else reg[key] * field
    return float(assemble(w * dx)) * RHO_GT


def smb_field(state):
    r"""The run's own SMB [m/yr] as a DG0 field: the RACMO climatology the
    forward is forced with, over ``ISMIP7_CLIM_START/END``."""
    from icepack2_tools.forcing import load_racmo_smb_climatology
    Q0 = state["thickness"].function_space()
    y0 = int(os.environ.get("ISMIP7_CLIM_START", "2000"))
    y1 = int(os.environ.get("ISMIP7_CLIM_END", "2029"))
    return load_racmo_smb_climatology(Q0, y0, y1)


def report(states, csv_path=None):
    a = states[0]
    reg = regions(a)
    area_gr = float(assemble(reg["grounded"] * dx)) / 1e6      # km^2
    area_fl = float(assemble(reg["floating"] * dx)) / 1e6
    mass = {k: integrate(a["thickness"], reg, k) for k in ("grounded", "floating")}
    vaf = float(assemble(max_value(reg["haf"], Constant(0.0)) * dx)) * MM_SLE

    print(f"state {os.path.basename(a['path'])}"
          + (f"  t = {a['t_yr']}" if a["t_yr"] else ""))
    print(f"  grounded {area_gr / 1e6:6.2f} M km2, {mass['grounded'] / 1e6:7.3f} "
          f"x 10^6 Gt      floating {area_fl / 1e6:5.2f} M km2, "
          f"{mass['floating'] / 1e6:6.3f} x 10^6 Gt")
    print(f"  VAF {vaf:.1f} mm SLE")

    q_gl = crossing_flux(a, reg, "grounded", "floating")
    # Grounded ice that flows straight into open water (a grounded front)
    open_water = Function(reg["ice"].function_space())
    open_water.dat.data[:] = 1.0 - reg["ice"].dat.data_ro
    reg["open"] = open_water
    q_gr_front = crossing_flux(a, reg, "grounded", "open")
    q_fl_front = crossing_flux(a, reg, "floating", "open")
    q_bnd = boundary_flux(a, reg, "ice")

    try:
        smb = smb_field(a)
        smb_gr = integrate(smb, reg, "grounded")
        smb_fl = integrate(smb, reg, "floating")
    except Exception as e:                        # forcing data may be absent
        print(f"  [SMB unavailable: {type(e).__name__}: {e}]")
        smb_gr = smb_fl = float("nan")

    print("\n  fluxes [Gt/yr]")
    print(f"    grounding line, grounded -> floating   {q_gl:8.0f}")
    print(f"    grounded front -> open water           {q_gr_front:8.0f}")
    print(f"    floating front -> open water           {q_fl_front:8.0f}")
    print(f"    across the mesh boundary               {q_bnd:8.0f}")
    print(f"    discharge (grounding line + grounded front) "
          f"{q_gl + q_gr_front:8.0f}   observed "
          f"{OBSERVED['grounding-line discharge'][0]:.0f} +/- "
          f"{OBSERVED['grounding-line discharge'][1]:.0f}")
    if a.get("velocity_obs") is not None:
        # The same facets, the same thickness, the observed velocity: this
        # separates "the model flows too slowly" from "the grounding line or
        # the thickness is in the wrong place".
        obs_state = dict(a, velocity=a["velocity_obs"])
        q_gl_o = crossing_flux(obs_state, reg, "grounded", "floating")
        q_fr_o = crossing_flux(obs_state, reg, "grounded", "open")
        print(f"    the same facets with the OBSERVED velocity "
              f"{q_gl_o + q_fr_o:8.0f}")
        print(f"      -> the model carries "
              f"{100 * (q_gl + q_gr_front) / max(q_gl_o + q_fr_o, 1e-9):.0f}% "
              f"of the observed-velocity flux across its own grounding line")

        # Where the deficit lives: the same fluxes binned by the observed
        # speed of the cell they leave.
        Q0 = reg["ice"].function_space()
        sp = Function(Q0)
        uo = a["velocity_obs"]
        sp.interpolate(fd.sqrt(fd.dot(uo, uo) + Constant(1e-12)))
        bins = [0.0, 100.0, 500.0, 1500.0, 1e9]
        m_bins = crossing_flux_binned(a, reg, "grounded", "floating", sp, bins)
        o_bins = crossing_flux_binned(a, reg, "grounded", "floating", sp, bins,
                                      velocity=uo)
        print("\n    grounding-line flux by observed speed of the source cell "
              "[Gt/yr]")
        print(f"      {'speed [m/yr]':>16s} {'model':>9s} {'observed':>9s} "
              f"{'model/obs':>10s}")
        labels = ["< 100", "100 - 500", "500 - 1500", "> 1500"]
        for lab, mv, ov in zip(labels, m_bins, o_bins):
            r = mv / ov if ov > 1e-9 else float("nan")
            print(f"      {lab:>16s} {mv:9.0f} {ov:9.0f} {r:10.2f}")

    print(f"\n  SMB [Gt/yr]   grounded {smb_gr:8.0f}   floating {smb_fl:8.0f}"
          f"   total {smb_gr + smb_fl:8.0f}")
    print(f"    observed grounded {OBSERVED['grounded SMB'][0]:.0f} +/- "
          f"{OBSERVED['grounded SMB'][1]:.0f}, total "
          f"{OBSERVED['total SMB'][0]:.0f} +/- {OBSERVED['total SMB'][1]:.0f}")

    print("\n  grounded budget [Gt/yr]")
    print(f"    SMB_gr - discharge = {smb_gr:.0f} - {q_gl + q_gr_front:.0f} "
          f"= {smb_gr - q_gl - q_gr_front:+.0f}")
    print(f"    that tendency is {(smb_gr - q_gl - q_gr_front) / 362.5:+.2f} "
          f"mm SLE/yr of VAF if it all sits above flotation")

    if len(states) > 1:
        b = states[1]
        reg_b = regions(b)
        dt = (b["t_yr"] or 0.0) - (a["t_yr"] or 0.0)
        if dt > 0:
            dm_gr = (integrate(b["thickness"], reg_b, "grounded")
                     - mass["grounded"]) / dt
            dm_fl = (integrate(b["thickness"], reg_b, "floating")
                     - mass["floating"]) / dt
            vaf_b = float(assemble(max_value(reg_b["haf"], Constant(0.0)) * dx)) * MM_SLE
            print(f"\n  differenced over {dt:.1f} yr to "
                  f"{os.path.basename(b['path'])}")
            print(f"    dM_gr/dt {dm_gr:+8.0f} Gt/yr    dM_fl/dt {dm_fl:+8.0f} Gt/yr"
                  f"    dVAF/dt {(vaf_b - vaf) / dt:+.2f} mm SLE/yr")
            print(f"    grounded residual (dM_gr/dt - [SMB_gr - discharge]) "
                  f"{dm_gr - (smb_gr - q_gl - q_gr_front):+.0f} Gt/yr")

    if csv_path and os.path.exists(csv_path):
        with open(csv_path) as f:
            rows = list(csv.DictReader(f))
        if rows:
            first, last = rows[0], rows[-1]
            dt = float(last["year"]) - float(first["year"])
            if dt > 0:
                dv = (float(last["vaf_mm_sle"]) - float(first["vaf_mm_sle"])) / dt
                dm = (float(last["mass_gt"]) - float(first["mass_gt"])) / dt
                print(f"\n  the run's own timeseries over {dt:.1f} yr: "
                      f"dM/dt {dm:+.0f} Gt/yr, dVAF/dt {dv:+.2f} mm SLE/yr")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("checkpoints", nargs="+")
    ap.add_argument("--csv", default=None, help="the run's timeseries, for context")
    ap.add_argument("--obs-velocity", default=None,
                    help="a MAP checkpoint whose velocity_obs is compared across "
                         "the same facets (default: the state's own, if present)")
    args = ap.parse_args()
    states = []
    for p in args.checkpoints[:2]:
        s = load_state(p)
        s["path"] = p
        states.append(s)
    if args.obs_velocity:
        with fd.CheckpointFile(args.obs_velocity, "r") as chk:
            mesh_o = chk.load_mesh()
            u_obs = chk.load_function(mesh_o, name="velocity_obs")
        # The dof copy below is raw and positional, so equal cell counts are
        # not enough: the two meshes must be the same mesh, or the observed
        # velocity lands on the wrong cells and every flux below is garbage.
        coords_o = mesh_o.coordinates.dat.data_ro
        coords_s = states[0]["mesh"].coordinates.dat.data_ro
        if (coords_o.shape != coords_s.shape
                or not np.allclose(coords_o, coords_s)):
            raise SystemExit(
                f"--obs-velocity {args.obs_velocity} is on a different mesh "
                f"than the state {states[0]['path']}: coordinates differ"
            )
        target = Function(states[0]["velocity"].function_space())
        target.dat.data[:] = u_obs.dat.data_ro
        states[0]["velocity_obs"] = target
    report(states, args.csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
