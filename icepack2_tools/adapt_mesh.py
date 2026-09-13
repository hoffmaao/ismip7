r"""Úa-style adaptive remeshing (``AdaptMesh``) for the icepack2 forward.

This is a port of the scheme Úa actually runs between run-steps, read from
UaSource (``UaMain/AdaptMesh.m``, ``NewDesiredEleSizesAndElementsToRefineOr
Coarsen2.m``, ``Error2EleSize.m``, ``GlobalRemeshing.m``,
``MapFbetweenMeshes.m``; Sep 2026). Úa's ``GLmorphing`` mesh-deformation code
was NOT ported: it has no caller in ``Ua.m``/``Ua2D.m``, its remeshing hook is
commented out as "broken anyhow", and the Compendium never mentions it.

The scheme, step by step, with Úa's names:

1. **Desired element size at the nodes of the current mesh** (Úa
   ``EleSizeDesired``). Start at ``MeshSizeMax``. For every enabled entry of
   ``ExplicitMeshRefinementCriteria`` compute a nodal error proxy ``e`` and map
   it with ``Error2EleSize``,
   ``h = hMin + (e0 / (e + e0))**(1/p) * (hMax - hMin)``  (``e0`` = ``Scale``),
   then take the minimum over criteria. If no criterion fired anywhere, use
   ``MeshSize``. Relax toward the current size with ``W = 0.5`` and clip the
   ratio to ``[MinRatioOfChange, MaxRatioOfChange]`` (1/5 .. 5). THEN the
   absolute bands: nodes within ``d_i`` of the grounding line get
   ``min(h, s_i)`` for each ``(d_i, s_i)`` row of ``MeshAdapt.GLrange``,
   floored at ``MeshSizeMin``; same for ``CFrange`` at the calving front.
2. **Global remeshing** (``explicit:global`` with gmsh): the nodal size field
   of the OLD mesh becomes gmsh's background scalar view, the domain outline is
   rebuilt exactly as ``mesh_antarctica.py`` builds it, and the mesh is
   regenerated. Úa then rescales ``MeshSizeMin`` up to four times so the
   element count lands within ``[LowerLimitFactor, UpperLimitFactor] *
   MaxNumberOfElements``.
3. **Transfer** (``MapFbetweenMeshes``): fields are interpolated old -> new by
   point evaluation (Úa: FE shape functions, nearest outside), thickness
   outside the old mesh is ``ThickMin``, the bed is re-sampled from data on
   the new mesh (Úa's ``DefineGeometry`` route) and the surface recomputed by
   flotation. ``transfer="project"`` swaps the DG0 thickness transfer for the
   conservative supermesh projection, which Úa (nodal) has no analogue of.

The forward's frozen anchors (``a_ref_mb``, ``N_ref``, ``C_w0``, ``H_init``,
the level set) are transferred too. A remapped ``a_ref_mb`` no longer cancels
the NEW mesh's discrete flux divergence exactly, which is a real cost of
remeshing under a frozen mass-balance correction; the transfer reports the
integrated mass change so it can be audited.
"""

import json
import os
from dataclasses import dataclass, field

import numpy as np
from firedrake import (
    And,
    Constant,
    Function,
    FunctionSpace,
    Mesh,
    SpatialCoordinate,
    TestFunction,
    VectorFunctionSpace,
    assemble,
    conditional,
    dS,
    ds,
    dx,
    grad,
    inner,
    jump,
    max_value,
    project,
    sqrt,
)
from firedrake.petsc import PETSc
from mpi4py import MPI

from .geometry import cg1_lift

RHO_I = 917.0
RHO_W = 1024.0

CRITERIA = (
    "effective strain rates",
    "effective strain rates gradient",
    "flotation",
    "thickness gradient",
    "upper surface gradient",
    "lower surface gradient",
    "|dhdt|",
    "dhdt gradient",
)


@dataclass
class Criterion:
    r"""One entry of Úa's ``CtrlVar.ExplicitMeshRefinementCriteria``."""
    name: str
    scale: float
    ele_min: float | None = None
    ele_max: float | None = None
    p: float = 1.0
    use: bool = True


@dataclass
class AdaptMeshConfig:
    r"""Úa's ``CtrlVar`` fields that shape one adaptation, same defaults."""
    mesh_size: float = 10e3            # CtrlVar.MeshSize
    mesh_size_min: float = 1e3         # CtrlVar.MeshSizeMin
    mesh_size_max: float = 10e3        # CtrlVar.MeshSizeMax
    gl_range: list = field(default_factory=list)   # CtrlVar.MeshAdapt.GLrange rows (d, h)
    cf_range: list = field(default_factory=list)   # CtrlVar.MeshAdapt.CFrange rows (d, h)
    criteria: list = field(default_factory=list)   # ExplicitMeshRefinementCriteria
    relaxation_w: float = 0.5          # W (hard-coded in Úa)
    max_ratio_change: float = 5.0      # MaxRatioOfChangeInEleSizeDuringAdaptMeshing
    min_ratio_change: float = 0.2      # MinRatioOfChangeInEleSizeDuringAdaptMeshing
    max_number_of_elements: int = 0    # MaxNumberOfElements (0 = no rescaling loop)
    upper_limit_factor: float = 1.3    # MaxNumberOfElementsUpperLimitFactor
    lower_limit_factor: float = 0.0    # MaxNumberOfElementsLowerLimitFactor
    thick_min: float = 1.0             # CtrlVar.ThickMin (OutsideValue.h)
    gl_threshold: float = 0.5          # CtrlVar.GLthreshold (grounded fraction)
    front_hmin: float = 1.0            # ice / no-ice threshold for the front
    refine_dirac_width: float = 100.0  # CtrlVar.RefineDiracDeltaWidth [m]
    refine_dirac_offset: float = 0.0   # CtrlVar.RefineDiracDeltaOffset
    transfer: str = "interpolate"      # "interpolate" (Úa) | "project" (conservative DG0)
    front_preserve: bool = True        # boundary cells take the nearest old boundary cell's h
    # The two rules of Úa's PIG-TWG DefineDesiredEleSize.m (UaExamples):
    #   EleSizeIndicator(GF.node<0.1)=UserVar.MeshSizeIceShelves  (= MeshSizeMax/5)
    #   EleSizeIndicator(s<1500)=CtrlVar.MeshSizeMax/5
    shelf_size: float | None = None    # floating ice gets min(h, shelf_size)
    low_surface: tuple | None = None   # (elevation [m], size): s below it gets min(h, size)
    # CtrlVar.MapOldToNew.Transient.Geometry: Úa's default moves the SURFACE
    # and derives h from it and the re-sampled bed ("bh-FROM-sBS"); the
    # alternative moves the thickness and derives the surface ("bs-FROM-hBS").
    # Moving h onto a finer bed puts a coarse-mean thickness on deep troughs
    # and floats it (measured 32 -> 8 km: outflux 782 -> 14770 Gt/yr).
    geometry: str = "bh-FROM-sBS"
    # Null test: desired size = the current size everywhere (before the bands),
    # i.e. remesh at the resolution the mesh already has.
    keep_current: bool = False

    @staticmethod
    def _rows(text):
        rows = []
        for item in text.split(","):
            item = item.strip()
            if not item:
                continue
            d, h = item.split(":")
            rows.append((float(d), float(h)))
        return rows

    @classmethod
    def from_env(cls):
        r"""Build from ``ISMIP7_ADAPT_*`` variables.

        ``ISMIP7_ADAPT_GL_RANGE="5000:2000,1000:500"`` is Úa's
        ``GLrange=[5000 2000; 1000 500]``; ``ISMIP7_ADAPT_CRITERIA=
        "effective strain rates:0.01,thickness gradient:0.01"`` lists
        criteria with their ``Scale`` (``:p`` and ``:min:max`` optional, e.g.
        ``flotation:0.001:1:500:2000``).
        """
        env = os.environ.get
        preset = env("ISMIP7_ADAPT_PRESET", "").lower()
        if preset:
            if preset != "ua":
                raise ValueError("ISMIP7_ADAPT_PRESET must be 'ua' or unset")
            # Úa's own Antarctic sizes. Whole-continent Úa (Úa-FESOM coupling,
            # GMD 18, 2025): 180 km in the interior, 4 km where strain rates
            # are high, 2 km at the grounding line, ~250,000 elements, remeshed
            # every step. The PIG-TWG example (UaExamples): MeshSize = Max/2,
            # MeshSizeMin = Max/20, ice shelves and low ground (s < 1500 m)
            # at Max/5, 'effective strain rates' Scale 0.001, initial
            # adaptation iterated 5x, MaxNumberOfElements 70e3 for a regional
            # box. MISMIP+: GLrange [20000 5000; 10000 2000; 5000 500].
            # The GL bands below are the inferred pan-Antarctic form of that
            # ("down to 2 km at the grounding line"); every other value is
            # Úa's verbatim.
            base = dict(mesh_size_max=180e3, mesh_size=90e3, mesh_size_min=2e3,
                        shelf_size=4e3, low_surface=(1500.0, 36e3),
                        gl_range=[(10e3, 4e3), (5e3, 2e3)],
                        max_number_of_elements=250_000,
                        criteria=[Criterion("effective strain rates", 0.001, ele_min=4e3)])
        else:
            base = {}
        cfg = cls(
            mesh_size=float(env("ISMIP7_ADAPT_MESH_SIZE", base.get("mesh_size", 10e3))),
            mesh_size_min=float(env("ISMIP7_ADAPT_MESH_SIZE_MIN", base.get("mesh_size_min", 1e3))),
            mesh_size_max=float(env("ISMIP7_ADAPT_MESH_SIZE_MAX", base.get("mesh_size_max", 10e3))),
            gl_range=cls._rows(env("ISMIP7_ADAPT_GL_RANGE")) if env("ISMIP7_ADAPT_GL_RANGE") is not None else base.get("gl_range", []),
            cf_range=cls._rows(env("ISMIP7_ADAPT_CF_RANGE", "")),
            shelf_size=float(env("ISMIP7_ADAPT_SHELF_SIZE")) if env("ISMIP7_ADAPT_SHELF_SIZE") else base.get("shelf_size"),
            low_surface=(tuple(float(t) for t in env("ISMIP7_ADAPT_LOW_SURFACE").split(":"))
                         if env("ISMIP7_ADAPT_LOW_SURFACE") else base.get("low_surface")),
            relaxation_w=float(env("ISMIP7_ADAPT_RELAXATION_W", 0.5)),
            max_ratio_change=float(env("ISMIP7_ADAPT_MAX_RATIO_CHANGE", 5.0)),
            min_ratio_change=float(env("ISMIP7_ADAPT_MIN_RATIO_CHANGE", 0.2)),
            max_number_of_elements=int(float(env("ISMIP7_ADAPT_MAX_ELEMENTS", base.get("max_number_of_elements", 0)))),
            thick_min=float(env("ISMIP7_ADAPT_THICK_MIN", 1.0)),
            front_hmin=float(env("ISMIP7_FRONT_HMIN", 1.0)),
            refine_dirac_width=float(env("ISMIP7_ADAPT_DIRAC_WIDTH", 100.0)),
            transfer=env("ISMIP7_ADAPT_TRANSFER", "interpolate"),
            geometry=env("ISMIP7_ADAPT_GEOMETRY", "bh-FROM-sBS"),
            keep_current=env("ISMIP7_ADAPT_KEEP_CURRENT", "0") == "1",
            front_preserve=env("ISMIP7_ADAPT_FRONT_PRESERVE", "1") == "1",
        )
        if cfg.geometry not in ("bh-FROM-sBS", "bs-FROM-hBS"):
            raise ValueError("ISMIP7_ADAPT_GEOMETRY must be bh-FROM-sBS or bs-FROM-hBS")
        crit_env = env("ISMIP7_ADAPT_CRITERIA")
        if crit_env is None:
            cfg.criteria = list(base.get("criteria", []))
        for item in (crit_env or "").split(","):
            item = item.strip()
            if not item:
                continue
            parts = item.split(":")
            name = parts[0].strip()
            if name not in CRITERIA:
                raise ValueError(f"unknown Úa refinement criterion {name!r}; "
                                 f"one of {CRITERIA}")
            c = Criterion(name=name, scale=float(parts[1]))
            if len(parts) > 2 and parts[2]:
                c.p = float(parts[2])
            if len(parts) > 4:
                c.ele_min, c.ele_max = float(parts[3]), float(parts[4])
            cfg.criteria.append(c)
        return cfg


# ---------------------------------------------------------------------------
# Úa building blocks
# ---------------------------------------------------------------------------

def error_to_ele_size(e, e0, h_min, h_max, p=1.0):
    r"""Úa ``Error2EleSize``: equidistribute ``e h^p`` between ``h_min`` and
    ``h_max``: ``h = h_min + (e0/(e+e0))^(1/p) (h_max - h_min)``."""
    e = np.asarray(e, dtype=float)
    return h_min + (e0 / (e + e0)) ** (1.0 / p) * (h_max - h_min)


def dirac_delta(k, x, x0=0.0):
    r"""Úa ``DiracDelta(k, x, x0) = 0.5 k sech^2(k (x - x0))`` (UaMain/DiracDelta.m)."""
    return k / (2.0 * np.cosh(np.clip(k * (x - x0), -300, 300)) ** 2)


def _nodal(expr, Qc):
    r"""Úa ``ProjectFintOntoNodes``: L2-project a cellwise expression to the
    CG1 nodes."""
    return Function(Qc).project(expr)


def _grad_mag_nodal(f_cg1, Qc):
    g = grad(f_cg1)
    return _nodal(sqrt(inner(g, g) + Constant(1e-30)), Qc)


def effective_strain_rate_nodal(u, Qc):
    r"""Úa ``CalcHorizontalNodalStrainRates``: ``e = sqrt(exx^2 + eyy^2 +
    exx eyy + exy^2)`` from the nodal velocity, projected to nodes."""
    exx = _nodal(grad(u)[0, 0], Qc)
    eyy = _nodal(grad(u)[1, 1], Qc)
    exy = _nodal(0.5 * (grad(u)[0, 1] + grad(u)[1, 0]), Qc)
    return Function(Qc).interpolate(
        sqrt(exx ** 2 + eyy ** 2 + exx * eyy + exy ** 2 + Constant(1e-30)))


def current_element_size_nodal(mesh, Qc):
    r"""Úa ``EleSizeCurrent = sqrt(M * EleArea)``: the square root of the mean
    area of the elements around each node. Note this is 0.66 of an
    equilateral edge length; Úa mixes this with edge-length targets in its
    relaxation and ratio limits, and so does this port, deliberately."""
    from firedrake import CellVolume
    Q0 = FunctionSpace(mesh, "DG", 0)
    area = Function(Q0).interpolate(CellVolume(mesh))
    mean_area = cg1_lift(area)
    out = Function(Qc)
    out.dat.data[:] = np.sqrt(np.maximum(mean_area.dat.data_ro, 1e-30))
    return out


def facet_midpoints_and_jumps(mesh, flags):
    r"""Per-facet midpoints plus the jump of each DG0 flag across interior
    facets, via an HDiv-trace test function (one dof per facet).

    Returns ``(xm, ym, {name: jump_values}, is_exterior)`` for this rank's
    owned facets. Exterior facets get the flag's own value in ``jumps``.
    """
    T = FunctionSpace(mesh, "HDiv Trace", 0)
    vt = TestFunction(T)
    x = SpatialCoordinate(mesh)
    # A trace test function is single-valued on a facet, but UFL still needs
    # an explicit restriction inside dS.
    vp = vt("+")
    length = assemble(vp * dS + vt * ds).dat.data_ro
    xm = assemble(x[0] * vp * dS + x[0] * vt * ds).dat.data_ro / length
    ym = assemble(x[1] * vp * dS + x[1] * vt * ds).dat.data_ro / length
    ext = assemble(vt * ds).dat.data_ro > 0.5 * length
    jumps = {}
    for name, f in flags.items():
        j = assemble(jump(f) * vp * dS + f * vt * ds).dat.data_ro / length
        jumps[name] = j
    return xm, ym, jumps, ext


def gather_points(comm, xs, ys):
    pts = np.column_stack([xs, ys]) if len(xs) else np.zeros((0, 2))
    allp = comm.allgather(pts)
    return np.vstack(allp) if allp else np.zeros((0, 2))


def grounding_line_points(mesh, H, b, cfg):
    r"""Midpoints of the facets separating grounded from floating ice cells
    (Úa's GL from the ``GF`` = ``GLthreshold`` crossing), gathered globally."""
    Q0 = FunctionSpace(mesh, "DG", 0)
    rr = Constant(RHO_I / RHO_W)
    ice = conditional(H > cfg.front_hmin, 1.0, 0.0)
    grounded = Function(Q0).interpolate(conditional(rr * H + b > 0.0, 1.0, 0.0) * ice)
    floating = Function(Q0).interpolate(conditional(rr * H + b <= 0.0, 1.0, 0.0) * ice)
    xm, ym, jumps, ext = facet_midpoints_and_jumps(mesh, {"g": grounded, "f": floating})
    # interior facets where one side is grounded ice and the other floating ice
    sel = (~ext) & (np.abs(jumps["g"]) > 0.5) & (np.abs(jumps["f"]) > 0.5)
    return gather_points(mesh.comm, xm[sel], ym[sel])


def calving_front_points(mesh, H, cfg):
    r"""Midpoints of the facets separating ice from no-ice cells, plus exterior
    facets of ice cells (the front on an unbuffered mesh), gathered globally."""
    Q0 = FunctionSpace(mesh, "DG", 0)
    ice = Function(Q0).interpolate(conditional(H > cfg.front_hmin, 1.0, 0.0))
    xm, ym, jumps, ext = facet_midpoints_and_jumps(mesh, {"i": ice})
    sel = ((~ext) & (np.abs(jumps["i"]) > 0.5)) | (ext & (jumps["i"] > 0.5))
    return gather_points(mesh.comm, xm[sel], ym[sel])


def nodal_distance_to(mesh, Qc, points):
    r"""Distance from every CG1 node to the nearest of ``points`` (Úa uses a
    KD-tree range search; the nearest-distance form is equivalent for bands)."""
    from scipy.spatial import cKDTree
    X = mesh.coordinates.dat.data_ro[:, :2]
    d = Function(Qc)
    if points.shape[0] == 0:
        d.dat.data[:] = np.inf
        return d
    d.dat.data[:] = cKDTree(points).query(X, k=1)[0]
    return d


# ---------------------------------------------------------------------------
# Step 1: desired element size
# ---------------------------------------------------------------------------

def desired_element_size(mesh, cfg, H, b, u=None, dhdt=None, log=PETSc.Sys.Print):
    r"""Úa ``NewDesiredEleSizesAndElementsToRefineOrCoarsen2`` for the
    ``explicit:global`` method. Returns the CG1 field ``EleSizeDesired``
    on the current mesh and a dict of diagnostics."""
    Qc = FunctionSpace(mesh, "CG", 1)
    comm = mesh.comm
    n_nodes = Qc.dof_dset.size
    h_des = Function(Qc)
    h_des.dat.data[:] = cfg.mesh_size_max
    fired = False

    s = Function(H.function_space()).interpolate(
        max_value(b + H, (Constant(1.0) - Constant(RHO_I / RHO_W)) * H))
    H_cg, b_cg, s_cg = cg1_lift(H), cg1_lift(b), cg1_lift(s)

    for c in cfg.criteria:
        if not c.use:
            continue
        h_min = cfg.mesh_size_min if c.ele_min is None else c.ele_min
        h_max = cfg.mesh_size_max if c.ele_max is None else c.ele_max
        name = c.name
        if name in ("effective strain rates", "effective strain rates gradient"):
            if u is None:
                log(f"  adapt: criterion {name!r} skipped (no velocity given)")
                continue
            e = effective_strain_rate_nodal(u, Qc)
            if name.endswith("gradient"):
                e = _grad_mag_nodal(e, Qc)
        elif name == "flotation":
            hf = Function(Qc).interpolate(Constant(RHO_W / RHO_I) * max_value(-b_cg, Constant(0.0)))
            e = Function(Qc)
            e.dat.data[:] = dirac_delta(1.0 / cfg.refine_dirac_width,
                                        H_cg.dat.data_ro - hf.dat.data_ro,
                                        cfg.refine_dirac_offset)
        elif name == "thickness gradient":
            e = _grad_mag_nodal(H_cg, Qc)
        elif name == "upper surface gradient":
            e = _grad_mag_nodal(s_cg, Qc)
        elif name == "lower surface gradient":
            e = _grad_mag_nodal(b_cg, Qc)
        elif name in ("|dhdt|", "dhdt gradient"):
            if dhdt is None:
                log(f"  adapt: criterion {name!r} skipped (no dhdt given)")
                continue
            dh = dhdt if dhdt.function_space() == Qc else cg1_lift(dhdt)
            e = Function(Qc).interpolate(abs(dh))
            if name == "dhdt gradient":
                e = _grad_mag_nodal(dh, Qc)
            emax = comm.allreduce(float(np.abs(e.dat.data_ro).max()) if n_nodes else 0.0, op=MPI.MAX)
            if emax < 1e-5:
                log(f"  adapt: criterion {name!r} too small to be of use, discarded")
                continue
        else:
            raise ValueError(name)
        h_c = error_to_ele_size(e.dat.data_ro, c.scale, h_min, h_max, c.p)
        h_des.dat.data[:] = np.minimum(h_des.dat.data_ro, h_c)
        fired = True
        lo = comm.allreduce(float(h_c.min()) if n_nodes else np.inf, op=MPI.MIN)
        log(f"  adapt: criterion {name!r} scale={c.scale:g} -> min desired size {lo:.0f} m")

    all_max = comm.allreduce(bool(np.all(h_des.dat.data_ro >= cfg.mesh_size_max - 1e-9)), op=MPI.LAND)
    if not fired or all_max:
        h_des.dat.data[:] = cfg.mesh_size
        log(f"  adapt: no relative criterion fired; desired size = MeshSize {cfg.mesh_size:g} m")

    # Relaxation toward the current size, then the ratio-of-change limits.
    h_cur = current_element_size_nodal(mesh, Qc)
    if cfg.keep_current:
        # Úa's EleSizeCurrent is sqrt(area); gmsh wants an edge length.
        h_des.dat.data[:] = h_cur.dat.data_ro * (4.0 / np.sqrt(3.0)) ** 0.5
        log("  adapt: keep_current: desired size = current element size")
    W = cfg.relaxation_w if not cfg.keep_current else 1.0
    h_des.dat.data[:] = W * h_des.dat.data_ro + (1.0 - W) * h_cur.dat.data_ro
    ratio = h_des.dat.data_ro / h_cur.dat.data_ro
    hi = ratio > cfg.max_ratio_change
    lo_ = ratio < cfg.min_ratio_change
    h_des.dat.data[hi] = cfg.max_ratio_change * h_cur.dat.data_ro[hi]
    h_des.dat.data[lo_] = cfg.min_ratio_change * h_cur.dat.data_ro[lo_]

    # Absolute bands from the grounding line and the calving front.
    diag = {}
    if cfg.gl_range:
        gl = grounding_line_points(mesh, H, b, cfg)
        d_gl = nodal_distance_to(mesh, Qc, gl)
        for (dist, size) in cfg.gl_range:
            size_eff = max(size, cfg.mesh_size_min)
            if size < cfg.mesh_size_min:
                log(f"  adapt: GLrange size {size:g} < MeshSizeMin {cfg.mesh_size_min:g}, using the latter")
            sel = d_gl.dat.data_ro < dist
            h_des.dat.data[sel] = np.minimum(h_des.dat.data_ro[sel], size_eff)
        diag["n_gl_points"] = int(gl.shape[0])
    if cfg.cf_range:
        cf = calving_front_points(mesh, H, cfg)
        d_cf = nodal_distance_to(mesh, Qc, cf)
        for (dist, size) in cfg.cf_range:
            size_eff = max(size, cfg.mesh_size_min)
            sel = d_cf.dat.data_ro < dist
            h_des.dat.data[sel] = np.minimum(h_des.dat.data_ro[sel], size_eff)
        diag["n_cf_points"] = int(cf.shape[0])

    # Úa's user hook (DefineDesiredEleSize) runs after the bands; PIG-TWG's is
    # two rules: floating ice and low ground get a fixed size.
    if cfg.shelf_size is not None or cfg.low_surface is not None:
        rr = Constant(RHO_I / RHO_W)
        if cfg.shelf_size is not None:
            fl = Function(H.function_space()).interpolate(
                conditional(And(rr * H + b <= 0.0, H > cfg.front_hmin), 1.0, 0.0))
            fl_n = cg1_lift(fl).dat.data_ro > 0.9          # Úa: GF.node < 0.1
            h_des.dat.data[fl_n] = np.minimum(h_des.dat.data_ro[fl_n], cfg.shelf_size)
        if cfg.low_surface is not None:
            elev, size = cfg.low_surface
            low = s_cg.dat.data_ro < elev
            h_des.dat.data[low] = np.minimum(h_des.dat.data_ro[low], size)
    h_des.dat.data[:] = np.clip(h_des.dat.data_ro, cfg.mesh_size_min, cfg.mesh_size_max)
    lo = comm.allreduce(float(h_des.dat.data_ro.min()) if n_nodes else np.inf, op=MPI.MIN)
    hi = comm.allreduce(float(h_des.dat.data_ro.max()) if n_nodes else -np.inf, op=MPI.MAX)
    diag.update(size_min=lo, size_max=hi)
    log(f"  adapt: desired element size in [{lo:.0f}, {hi:.0f}] m"
        + (f", GL points {diag.get('n_gl_points')}" if "n_gl_points" in diag else "")
        + (f", front points {diag.get('n_cf_points')}" if "n_cf_points" in diag else ""))
    return h_des, diag


# ---------------------------------------------------------------------------
# Step 2: global remeshing with gmsh (rank 0)
# ---------------------------------------------------------------------------

def _gather_nodal_field(mesh, f):
    comm = mesh.comm
    X = mesh.coordinates.dat.data_ro[:, :2]
    v = f.dat.data_ro
    allX = comm.gather(np.asarray(X), root=0)
    allv = comm.gather(np.asarray(v), root=0)
    if comm.rank == 0:
        return np.vstack(allX), np.concatenate(allv)
    return None, None


def _physical_groups(msh_path):
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 1)
    gmsh.open(msh_path)
    out = {}
    for dim, tag in gmsh.model.getPhysicalGroups(1):
        out[gmsh.model.getPhysicalName(dim, tag)] = tag
    gmsh.finalize()
    return out


def remesh_global(mesh, h_des, cfg, out_msh, old_msh, old_sidecar, log=PETSc.Sys.Print):
    r"""Úa ``GlobalRemeshing`` with the gmsh generator: background size field
    from the old mesh's nodes, the same domain outline and physical groups as
    ``mesh_antarctica.py``, and the ``MaxNumberOfElements`` rescaling loop.

    Writes ``out_msh`` (gmsh 2.2) and its boundary-id sidecar. Returns the
    element count. Collective; gmsh runs on rank 0.
    """
    comm = mesh.comm
    X, v = _gather_nodal_field(mesh, h_des)
    n_ele = None
    if comm.rank == 0:
        import gmsh
        from scipy.spatial import Delaunay
        from .mesh import (build_gmsh_geometry, classify_boundaries,
                           extract_ice_outline, load_bedmachine_mask)

        # mesh.py's own DATA_DIR points at <repo>/data; the BedMachine copy
        # every other script uses lives under antarctica/data.
        data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "antarctica", "data")
        mask, x, y = load_bedmachine_mask(data_dir)
        outline = extract_ice_outline(mask, x, y)
        boundaries, names = classify_boundaries(outline, mask, x, y)

        tri = Delaunay(X).simplices
        sizes = v.copy()
        size_min = cfg.mesh_size_min

        def generate(sizes):
            gmsh.initialize()
            gmsh.option.setNumber("General.Verbosity", 1)
            gmsh.model.add("adapt")
            build_gmsh_geometry(boundaries, names, cfg.mesh_size_min, cfg.mesh_size_max)
            view = gmsh.view.add("ua size field")
            P = X[tri]                      # (ntri, 3, 2)
            S = sizes[tri]                  # (ntri, 3)
            data = np.column_stack([P[:, :, 0], P[:, :, 1], np.zeros_like(S), S]).ravel()
            gmsh.view.addListData(view, "ST", tri.shape[0], data.tolist())
            bg = gmsh.model.mesh.field.add("PostView")
            gmsh.model.mesh.field.setNumber(bg, "ViewTag", view)
            gmsh.model.mesh.field.setAsBackgroundMesh(bg)
            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
            gmsh.model.mesh.generate(2)
            tags, _ = gmsh.model.mesh.getElementsByType(2)
            gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
            gmsh.write(out_msh)
            groups = {gmsh.model.getPhysicalName(d, t): t for d, t in gmsh.model.getPhysicalGroups(1)}
            gmsh.finalize()
            return len(tags), groups

        n_ele, groups = generate(sizes)
        log(f"  adapt: remeshed -> {n_ele} elements")
        # Úa's element-count control: rescale MeshSizeMin, up to 4 times.
        if cfg.max_number_of_elements > 0:
            it = 0
            N = cfg.max_number_of_elements
            while ((n_ele > cfg.upper_limit_factor * N or n_ele < cfg.lower_limit_factor * N) and it < 4):
                scaling = np.sqrt(n_ele / N)
                max_e, min_e = float(sizes.max()), float(sizes.min())
                new_min = min(min_e * scaling, 0.9 * cfg.mesh_size_max)
                if max_e != min_e:
                    sizes = new_min + (cfg.mesh_size_max - new_min) * (sizes - min_e) / (max_e - min_e)
                else:
                    sizes = np.full_like(sizes, new_min)
                it += 1
                if new_min != size_min:
                    log(f"  adapt: MeshSizeMin rescaled {size_min:.0f} -> {new_min:.0f} m "
                        f"(Nele {n_ele} vs MaxNumberOfElements {N})")
                size_min = new_min
                n_ele, groups = generate(sizes)
                log(f"  adapt: remeshed -> {n_ele} elements")
        # Boundary ids: the new mesh must expose the same physical groups as
        # the old one, or the calving BC lands on the wrong facets.
        old_groups = _physical_groups(old_msh)
        if groups != old_groups:
            raise RuntimeError("remesh produced different physical groups than the "
                               f"old mesh: {len(groups)} vs {len(old_groups)}")
        with open(old_sidecar) as f:
            side = json.load(f)
        new_sidecar = os.path.join(os.path.dirname(out_msh),
                                   "boundary_ids_" + os.path.splitext(os.path.basename(out_msh))[0] + ".json")
        with open(new_sidecar, "w") as f:
            json.dump(side, f)
    comm.barrier()
    return comm.bcast(n_ele, root=0)


# ---------------------------------------------------------------------------
# Step 3: transfer (MapFbetweenMeshes)
# ---------------------------------------------------------------------------

def _boundary_cells(mesh, Q0):
    r"""Owned DG0 dofs of cells with an exterior facet, plus the midpoints of
    those facets (one entry per boundary facet; a corner cell appears twice)."""
    from firedrake import CellVolume
    T = FunctionSpace(mesh, "HDiv Trace", 0)
    vt = TestFunction(T)
    x = SpatialCoordinate(mesh)
    length = assemble(vt * ds).dat.data_ro
    ext = length > 0.0
    xm = assemble(x[0] * vt * ds).dat.data_ro[ext] / length[ext]
    ym = assemble(x[1] * vt * ds).dat.data_ro[ext] / length[ext]
    # which cell owns each exterior facet: mark cells by an exterior-facet
    # integral of a DG0 test function
    phi = TestFunction(Q0)
    on_bnd = assemble(phi * ds).dat.data_ro > 0.0
    cells = np.flatnonzero(on_bnd)
    # cell centroids for matching facets to cells
    xc = Function(Q0).interpolate(x[0]).dat.data_ro[cells]
    yc = Function(Q0).interpolate(x[1]).dat.data_ro[cells]
    return cells, xc, yc, xm, ym


def _preserve_front(mesh_old, H_old, mesh_new, H_new, cfg):
    r"""Overwrite the thickness of every new boundary cell with that of the
    nearest old boundary cell (nearest by cell centroid, globally)."""
    from scipy.spatial import cKDTree
    comm = mesh_new.comm
    Q_old, Q_new = H_old.function_space(), H_new.function_space()
    oc, oxc, oyc, _, _ = _boundary_cells(mesh_old, Q_old)
    old_pts = gather_points(comm, oxc, oyc)
    old_vals = np.concatenate(comm.allgather(np.asarray(H_old.dat.data_ro[oc])))
    if old_pts.shape[0] == 0:
        return 0
    nc, nxc, nyc, _, _ = _boundary_cells(mesh_new, Q_new)
    if nc.size:
        _, j = cKDTree(old_pts).query(np.column_stack([nxc, nyc]), k=1)
        H_new.dat.data[nc] = old_vals[j]
    return int(comm.allreduce(nc.size))


def fv_flux_divergence(mesh, h_dg, u):
    r"""The forward's DG0 upwind facet-flux divergence, per unit area
    (simulation.py's ``flux0 / cell_area``): the discrete ``div(h u)``."""
    from firedrake import CellVolume, FacetNormal, dot
    Q0 = h_dg.function_space()
    phi = TestFunction(Q0)
    n = FacetNormal(mesh)
    un = dot(u, n)
    unp = (un + abs(un)) / 2
    flux = assemble((unp("+") * h_dg("+") - unp("-") * h_dg("-")) * jump(phi) * dS
                    + unp * h_dg * phi * ds)
    area = Function(Q0).interpolate(CellVolume(mesh))
    out = Function(Q0, name="flux_div")
    out.dat.data[:] = flux.dat.data_ro / area.dat.data_ro
    return out


def physical_divergence(mesh, h_dg, u, a_ref):
    r"""``P = flux/area - a_ref``: the flux divergence net of the frozen
    correction, i.e. the divergence the run was physically experiencing."""
    d = fv_flux_divergence(mesh, h_dg, u)
    out = Function(h_dg.function_space(), name="phys_div")
    out.dat.data[:] = d.dat.data_ro - a_ref.dat.data_ro
    return out


_DEFAULTS = {"thickness": None, "thickness_dg": None, "H_init": None}


def _xfer(f_old, mesh_new, how="interpolate", default=0.0):
    r"""Move ``f_old`` onto ``mesh_new`` in a space of the same element.
    ``interpolate`` is Úa's point evaluation (``default`` outside the old
    mesh); ``project`` is the conservative supermesh projection."""
    from firedrake import interpolate as _interp
    V_new = FunctionSpace(mesh_new, f_old.function_space().ufl_element())
    if how == "project":
        # Firedrake's supermesh projection refuses two independently
        # partitioned meshes in parallel ("Whoever made mesh_B should
        # explicitly mark mesh_A as having a compatible parallel layout"), and
        # forcing the mark would silently miss overlaps across rank
        # boundaries. Conservative transfer therefore runs on one rank.
        if mesh_new.comm.size > 1:
            raise RuntimeError("ISMIP7_ADAPT_TRANSFER=project needs the adapt step on ONE "
                               "rank (mpiexec -n 1); the forward can still run in parallel")
        return project(f_old, V_new)
    f_new = assemble(_interp(f_old, V_new, allow_missing_dofs=True,
                             default_missing_val=default))
    return f_new


def transfer_state(chk_in, mesh_new, cfg, chk_out, new_msh_basename, bed_sampler,
                   rebuild_aref=False, thickness_sampler=None, log=PETSc.Sys.Print):
    r"""Úa ``MapFbetweenMeshes`` onto ``mesh_new`` and write ``chk_out``.

    ``bed_sampler(Q_g_new, Q_cg_new)`` returns the bed on the new mesh from
    data (Úa's ``DefineGeometry`` route); ``thickness_sampler`` likewise, used
    only for the initial adaptation (``rebuild_aref``), where Úa takes ALL
    geometry from data. Otherwise ``cfg.geometry`` selects Úa's
    ``MapOldToNew.Transient.Geometry`` route; everything else is interpolated.
    Returns a dict of mass audits.
    """
    import firedrake as fd
    with fd.CheckpointFile(chk_in, "r") as chk:
        mesh_old = chk.load_mesh()
        attrs = {k: chk.get_attr("/", k) for k in
                 ("t_yr", "friction", "geometry_space", "mesh_basename", "lc",
                  "lc_coarse", "buffer_m", "raster_sample", "adapt_count")
                 if chk.has_attr("/", k)}
        names = ["log_friction", "log_fluidity", "fluidity_prior", "thickness",
                 "bed", "surface", "velocity", "membrane_stress", "basal_stress",
                 "H_init", "phi_eff", "C_w0", "N_ref", "a_ref_mb", "levelset",
                 "thickness_dg"]
        old = {}
        for name in names:
            try:
                old[name] = chk.load_function(mesh_old, name=name)
            except Exception:
                pass
    log(f"  adapt: transferring {sorted(old)} ({cfg.transfer})")

    geom = attrs.get("geometry_space", "dg0")
    Q_g = FunctionSpace(mesh_new, "DG" if geom == "dg0" else "CG", 0 if geom == "dg0" else 1)
    Qc = FunctionSpace(mesh_new, "CG", 1)
    rr = Constant(RHO_I / RHO_W)

    new = {}
    H_old = old["thickness"]
    b = bed_sampler(Q_g, Qc)
    if rebuild_aref and thickness_sampler is not None:
        # Úa, first run-step of a time-dependent run: ALL geometry comes from
        # DefineGeometry on the new mesh, nothing is interpolated.
        H = thickness_sampler(Q_g, Qc)
        route = "data"
    elif cfg.geometry == "bh-FROM-sBS":
        # Úa's default: move the surface, re-sample the bed, derive h.
        # Úa's s is nodal, so its interpolation is smooth; the DG0 analogue is
        # the volume-preserving CG1 lift of the cell surface, point-evaluated
        # at the new cell centroids. Then h = s - b where that grounds, else
        # the flotation thickness s * rho_w / (rho_w - rho_i).
        s_old = old.get("surface")
        if s_old is None:
            s_old = Function(H_old.function_space()).interpolate(
                max_value(old["bed"] + H_old, (Constant(1.0) - rr) * H_old))
        s_lift = cg1_lift(s_old) if geom == "dg0" else s_old
        s_new = _xfer(s_lift, mesh_new, "interpolate", default=0.0)
        s_c = Function(Q_g).interpolate(s_new)
        h_ground = s_c.dat.data_ro - b.dat.data_ro
        h_float = s_c.dat.data_ro * RHO_W / (RHO_W - RHO_I)
        H = Function(Q_g)
        H.dat.data[:] = np.where(s_c.dat.data_ro > 0.0,
                                 np.maximum(np.minimum(h_ground, h_float), 0.0), 0.0)
        route = "bh-FROM-sBS"
    else:
        # "bs-FROM-hBS": move the thickness, derive the surface.
        how_h = cfg.transfer if geom == "dg0" else "interpolate"
        H = _xfer(H_old, mesh_new, how_h, default=cfg.thick_min)
        H.dat.data[:] = np.maximum(H.dat.data_ro, 0.0)
        route = "bs-FROM-hBS"
    if geom == "dg0" and route != "data" and cfg.front_preserve:
        # Úa's nodal transfer keeps boundary values because boundary nodes
        # interpolate from the old boundary edge alone. The DG0 analogue:
        # every new cell on the exterior boundary takes the thickness of the
        # nearest OLD boundary cell (by facet midpoint), so the thin front is
        # neither pulled up by interior values (the CG1 lift's boundary bias)
        # nor by a centroid landing in an old interior cell.
        n_fixed = _preserve_front(mesh_old, H_old, mesh_new, H, cfg)
        log(f"  adapt: front-preserving transfer set {n_fixed} boundary cells")
    s = Function(Q_g).interpolate(max_value(b + H, (Constant(1.0) - rr) * H))
    new["thickness"], new["bed"], new["surface"] = H, b, s
    log(f"  adapt: geometry route {route}")
    for name, f in old.items():
        if name in ("thickness", "bed", "surface"):
            continue
        if name == "a_ref_mb" and rebuild_aref:
            continue
        if name == "H_init" and rebuild_aref:
            # the t=0 extent anchor IS the (re-sampled) initial thickness
            new[name] = Function(Q_g).assign(H)
            continue
        how = cfg.transfer if (name in ("a_ref_mb", "thickness_dg", "H_init") and geom == "dg0") else "interpolate"
        dflt = cfg.thick_min if name in ("H_init", "thickness_dg") else 0.0
        new[name] = _xfer(f, mesh_new, how, default=dflt)

    # The apparent-mass-balance reference cancels the OLD mesh's discrete flux
    # divergence spike by spike (~1000 m/yr locally at the PIG grounding zone),
    # so moved to another mesh it becomes a field of misplaced sources and the
    # run blows up within a year (measured: outflux 782 -> 53,185 Gt/yr). What
    # transfers is the PHYSICAL divergence the balanced run had,
    #     P = flux/area - a_ref   (= SMB - melt at t=0 in balance mode),
    # from which the forward rebuilds a_ref = flux_new/area - P with its own
    # operator on the new mesh, exactly as it built the t=0 one.
    if "a_ref_mb" in old and not rebuild_aref and "velocity" in old:
        P_old = physical_divergence(mesh_old, old["thickness"], old["velocity"], old["a_ref_mb"])
        # P is the SMOOTH physical field (the spikes live in a_ref, which it
        # nets out), so point evaluation transfers it well and, unlike the
        # supermesh projection, works between independently partitioned
        # meshes in parallel.
        new["phys_div"] = _xfer(P_old, mesh_new, "interpolate")
        new.pop("a_ref_mb", None)
        comm = mesh_new.comm
        pmin = comm.allreduce(float(P_old.dat.data_ro.min()) if P_old.dat.data_ro.size else np.inf, op=MPI.MIN)
        pmax = comm.allreduce(float(P_old.dat.data_ro.max()) if P_old.dat.data_ro.size else -np.inf, op=MPI.MAX)
        amin = comm.allreduce(float(old["a_ref_mb"].dat.data_ro.min()) if P_old.dat.data_ro.size else np.inf, op=MPI.MIN)
        amax = comm.allreduce(float(old["a_ref_mb"].dat.data_ro.max()) if P_old.dat.data_ro.size else -np.inf, op=MPI.MAX)
        log(f"  adapt: a_ref_old in [{amin:.1f}, {amax:.1f}] m/yr; physical divergence P in [{pmin:.2f}, {pmax:.2f}] m/yr, "
            f"net {assemble(P_old * dx) / 1e9:+.1f} Gt/yr")
        log("  adapt: a_ref_mb replaced by the transferred physical divergence (phys_div); "
            "the forward rebuilds a_ref on this mesh")

    # Budd's N_ref is the effective pressure of the ORIGINAL geometry frozen at
    # t=0; the law reads N_hat = N_eff/N_ref (cap 3, 0 afloat). Re-sampling the
    # bed changes N_eff cell by cell while an interpolated N_ref does not
    # follow, so friction is scrambled: measured on a same-resolution null
    # remesh the grounded-cell ratio spread went from 0.91-1.48 (5th-95th
    # pct) to 0.36-2.97 with 10% of cells losing their reference (=> capped
    # 3x friction). What is physical is the RATIO the run had, so transfer
    # that and rebuild N_ref = N_eff_new / ratio on the new geometry.
    # Regularized Coulomb has no such reference (its cap uses the live N).
    if "N_ref" in old and "N_ref" in new:
        # dual_friction.effective_pressure, replicated: importing that module
        # here pulls icepack2 -> irksome, which refuses to load after any UFL
        # form has been assembled (IrksomeImportOrderException). Same
        # constants (icepack2.constants: rho*g in MPa/m, the year factors
        # cancel), so N_ref stays in the units the law divides by.
        def effective_pressure(H_, s_):
            Hs = max_value(H_, Constant(1.0))
            p_I = Constant(917.0 * 9.81e-6) * Hs
            p_W = Constant(1024.0 * 9.81e-6) * max_value(Constant(0.0), H_ - s_)
            return max_value(p_I - p_W, Constant(0.0))
        nhat_cap = 3.0
        s_old_f = old.get("surface")
        if s_old_f is None:
            s_old_f = Function(H_old.function_space()).interpolate(
                max_value(old["bed"] + H_old, (Constant(1.0) - rr) * H_old))
        N_old = Function(H_old.function_space()).interpolate(
            max_value(effective_pressure(H_old, s_old_f), Constant(0.0)))
        ratio_old = Function(H_old.function_space(), name="nhat")
        nr = old["N_ref"].dat.data_ro
        ratio_old.dat.data[:] = np.where(nr > 1e-6, np.minimum(N_old.dat.data_ro / np.maximum(nr, 1e-6), nhat_cap), nhat_cap)
        ratio_new = _xfer(ratio_old, mesh_new, "interpolate", default=1.0)
        N_new = Function(Q_g).interpolate(max_value(effective_pressure(H, s), Constant(0.0)))
        rn = np.clip(ratio_new.dat.data_ro, 1e-3, nhat_cap)
        N_ref_new = Function(Q_g, name="N_ref")
        N_ref_new.dat.data[:] = np.where(N_new.dat.data_ro > 0.0, N_new.dat.data_ro / rn, 0.0)
        new["N_ref"] = N_ref_new
        log("  adapt: N_ref rebuilt from the transferred N_eff/N_ref ratio on the new geometry")

    m_old = assemble(H_old * dx) / 1e9
    m_new = assemble(H * dx) / 1e9
    # Mean thickness of the cells on the exterior boundary where ice exists:
    # the number the calving-terminus traction (~h^2) and the outflux see.
    def _front(Hf):
        ice = conditional(Hf > 1.0, 1.0, 0.0)
        return assemble(Hf * ds(domain=Hf.function_space().mesh())) / max(assemble(ice * ds(domain=Hf.function_space().mesh())), 1.0)
    audit = {"volume_old_km3": m_old, "volume_new_km3": m_new,
             "volume_change_pct": 100.0 * (m_new - m_old) / max(m_old, 1e-30),
             "front_h_old": _front(H_old), "front_h_new": _front(H)}
    log(f"  adapt: front <h> {audit['front_h_old']:.1f} -> {audit['front_h_new']:.1f} m")
    if "a_ref_mb" in old:
        audit["aref_old"] = assemble(old["a_ref_mb"] * dx) / 1e9
        if "a_ref_mb" in new:
            audit["aref_new"] = assemble(new["a_ref_mb"] * dx) / 1e9
    log(f"  adapt: ice volume {m_old:.6e} -> {m_new:.6e} km3 ({audit['volume_change_pct']:+.4f}%)")

    with fd.CheckpointFile(chk_out, "w") as chk:
        chk.save_mesh(mesh_new)
        for name, f in new.items():
            f.rename(name)
            chk.save_function(f, name=name)
        for k, v in attrs.items():
            chk.set_attr("/", k, v)
        chk.set_attr("/", "mesh_basename", new_msh_basename)
        chk.set_attr("/", "adapted_from", os.path.basename(chk_in))
        chk.set_attr("/", "adapt_count", int(attrs.get("adapt_count", 0)) + 1)
        if rebuild_aref:
            chk.set_attr("/", "adapted_initial", 1)
    return audit
