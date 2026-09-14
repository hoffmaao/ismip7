#!/usr/bin/env python
r"""ISMIP7 output, part two: the model's yearly fields onto the 8 km AIS grid
in the request's files.

    python antarctica/scripts/write_ismip7_output.py ANNUAL.h5 --out-dir DIR
        --esm CESM2-WACCM --scenario ssp585 --exp C007
        [--source-id RICE] [--ism-id icepack2] [--set-id CORE] [--scalars CSV]

Serial. ANNUAL.h5 is the STEM ``<experiment>_<lc>_ismip7_annual.h5`` that
``icepack2_tools.ismip7_output`` names its output after; the years themselves
are the files ``<experiment>_<lc>_ismip7_annual_<year>.h5`` beside it (one per
year, each written atomically, so an interrupted run cannot damage the years
already banked) and this globs them in year order. It writes one NetCDF per
variable under
``DIR/AIS/<source_id>/<ism_id>/<set_id>/<exp>/``, named
``<var>_AIS_<source_id>_<ism_id>_m001_<ESM>_f001_<scenario>_<exp>_<y0>-<y1>.nc``.

Regridding is conservative: a supermesh mixed mass matrix between the model's
DG0 cells and a triangulated copy of the 8 km grid gives the exact area of
every (cell, pixel) overlap; a pixel's value is the area-weighted mean of the
cells under it. Which area the mean is taken over follows the request's fill
policy (``isschecker/data/ISMIP7_variable_request.csv``): ``forbidden``
(thickness, fluxes, fractions) means over the whole pixel with the uncovered
part counting as zero, so sums over the grid are the model's sums;
``outside_domain`` (elevations) means over the covered part and fills pixels
the model does not cover; ``no_ice`` and friends mean over the ice part and
fill pixels without it. The overlap operator is cached next to the input
(``<annual>.overlap.npz``) because it depends on the mesh only.

Model-to-SI conversions use icepack's year, 365.25 days (31557600 s), which
is the model's own time unit; the time axis in the files is the standard
calendar regardless.

``acabf`` is written as the forcing surface mass balance, always. The
apparent-mass-balance reference stays where the forward put it, as
``acabf_correction`` in the annual file: it is not a request variable, and
folded into the SMB it would sit two orders of magnitude outside the
request's range. A reader who wants a grid budget that closes adds the two
from the annual file; the submission files never carry the sum.

Time follows ismip/ismip7-time-encoding: ``days since 1850-01-01`` on the
standard calendar; state variables are stamped 1 January of the following
year, fluxes 1 July with bounds over the year; the initial state is not
written.
"""
import argparse
import csv
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(_ROOT)))

import icepack2_tools.dual_friction  # noqa: F401,E402  (icepack2 -> irksome import order)
from icepack2_tools.ismip7_output import (AnnualOutput, RHO_I, SCALARS, SECONDS_PER_YEAR,  # noqa: E402
                                          VARIABLES_2D, VARIABLES_BANKED)
from icepack2_tools.regrid import ISMIP7_DX, ISMIP7_NX, ISMIP7_NY, ISMIP7_X0, ISMIP7_Y0  # noqa: E402
import firedrake as fd  # noqa: E402

import netCDF4 as _nc4
FILL = _nc4.default_fillvals["f4"]          # the checker wants the netCDF4 default fill for the dtype
TIME_UNITS = "days since 1850-01-01"        # the checker's exact spelling
PIXEL_AREA = ISMIP7_DX * ISMIP7_DX
REQUEST = os.path.join(os.path.dirname(os.path.dirname(_ROOT)), "icepack2_tools", "ismip7_variable_request.csv")   # not under a data/ dir: .gitignore ignores those


def standard_name(meta):
    r"""The CF standard name to write, or "" to omit the attribute.

    A blank standard_name column is the request stating that no standard name
    is assigned to that variable, so the attribute is omitted rather than
    written empty (an empty standard_name is not valid CF). The Comment column
    proposes names for some of them ("Should get standard name: ..."), but a
    proposal is not a registered CF name and a CF checker rejects it, so it is
    not used.
    """
    return (meta.get("standard_name") or "").strip()


def request_table():
    with open(REQUEST) as f:
        rows = list(csv.DictReader(f))
    return {r["Variable Name"]: r for r in rows}


# model units -> request units
CONVERT = {
    "m": 1.0, "1": 1.0,
    "m s-1": 1.0 / SECONDS_PER_YEAR,                      # m/yr -> m/s
    "kg m-2 s-1": RHO_I / SECONDS_PER_YEAR,               # m ice/yr -> kg m-2 s-1
    "Pa": 1.0e6,                                          # MPa -> Pa
}


def grid_mesh():
    r"""A triangulated copy of the 8 km grid whose pixel centres are the
    ISMIP7 points, and the pixel index of every triangle."""
    half = ISMIP7_DX / 2
    Lx = ISMIP7_X0 + (ISMIP7_NX - 1) * ISMIP7_DX + half
    Ly = ISMIP7_Y0 + (ISMIP7_NY - 1) * ISMIP7_DX + half
    mesh = fd.RectangleMesh(ISMIP7_NX, ISMIP7_NY, Lx, Ly, originX=ISMIP7_X0 - half, originY=ISMIP7_Y0 - half,
                            diagonal="left", reorder=False)
    Q = fd.FunctionSpace(mesh, "DG", 0)
    xc = fd.Function(fd.VectorFunctionSpace(mesh, "DG", 0)).interpolate(fd.SpatialCoordinate(mesh)).dat.data_ro
    i = np.floor((xc[:, 0] - (ISMIP7_X0 - half)) / ISMIP7_DX).astype(int)
    j = np.floor((xc[:, 1] - (ISMIP7_Y0 - half)) / ISMIP7_DX).astype(int)
    return mesh, Q, j * ISMIP7_NX + i                     # pixel = row-major (y, x)


def overlap_operator(mesh_src, cache):
    r"""Sparse (npixels x ncells) matrix of overlap areas, cached."""
    import scipy.sparse as sp
    if cache and os.path.exists(cache):
        z = np.load(cache)
        return sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    from firedrake.supermeshing import assemble_mixed_mass_matrix
    Q_src = fd.FunctionSpace(mesh_src, "DG", 0)
    gmesh, Q_g, pix = grid_mesh()
    M = assemble_mixed_mass_matrix(Q_src, Q_g)            # rows: grid triangles, cols: model cells
    indptr, indices, data = M.getValuesCSR()
    Mt = sp.csr_matrix((data, indices, indptr), shape=(Q_g.dim(), Q_src.dim()))
    P = sp.csr_matrix((np.ones(len(pix)), (pix, np.arange(len(pix)))), shape=(ISMIP7_NX * ISMIP7_NY, len(pix)))
    W = (P @ Mt).tocsr()
    if cache:
        np.savez(cache, data=W.data, indices=W.indices, indptr=W.indptr, shape=np.array(W.shape))
    return W


def regrid(W, values, policy, mask=None):
    r"""Pixel values under the request's fill policy; NaN where filled."""
    num = W @ values
    if policy == "forbidden":
        return num / PIXEL_AREA
    if policy == "outside_domain":
        # mean over the covered part of the pixel; any coverage counts, so
        # every pixel that carries ice (sftgif > 0) also carries elevations
        cov = W @ np.ones_like(values)
        out = np.full(num.shape, np.nan)
        ok = cov > 0.0
        out[ok] = num[ok] / cov[ok]
        return out
    # no_ice / no_grounded_ice / no_floating_ice: mean over the masked part
    m = mask.astype(float)
    num = W @ (values * m); den = W @ m
    out = np.full(num.shape, np.nan)
    ok = den > 0.0
    out[ok] = num[ok] / den[ok]
    return out


GLOBAL = {}


def _global_attrs(ds, meta):
    ds.Conventions = "CF-1.8"
    for k, v in GLOBAL.items():
        setattr(ds, k, v)


def days_since_1850(year, month, day):
    import cftime
    return cftime.date2num(cftime.datetime(year, month, day, calendar="standard"),
                           "days since 1850-01-01 00:00:00", calendar="standard")


def create_2d(path, var, meta, years, is_flux):
    r"""Create the file with its axes, crs and attributes, and return the
    open dataset and the empty variable to write time slices into. The years
    are written one at a time so a 286-year submission never holds more than
    one year of grid in memory."""
    import netCDF4
    x = ISMIP7_X0 + ISMIP7_DX * np.arange(ISMIP7_NX)
    y = ISMIP7_Y0 + ISMIP7_DX * np.arange(ISMIP7_NY)
    ds = netCDF4.Dataset(path, "w", format="NETCDF4")
    ds.createDimension("time", None); ds.createDimension("y", ISMIP7_NY); ds.createDimension("x", ISMIP7_NX)
    if is_flux:
        ds.createDimension("nv", 2)
    vx = ds.createVariable("x", "f8", ("x",)); vx[:] = x
    vx.units = "m"; vx.standard_name = "projection_x_coordinate"; vx.long_name = "x coordinate of projection"; vx.axis = "X"
    vy = ds.createVariable("y", "f8", ("y",)); vy[:] = y
    vy.units = "m"; vy.standard_name = "projection_y_coordinate"; vy.long_name = "y coordinate of projection"; vy.axis = "Y"
    vt = ds.createVariable("time", "f4", ("time",))
    vt.units = TIME_UNITS; vt.calendar = "standard"; vt.standard_name = "time"; vt.long_name = "time"; vt.axis = "T"
    if is_flux:
        vt.bounds = "time_bnds"
        vb = ds.createVariable("time_bnds", "f4", ("time", "nv"))
        vt[:] = [days_since_1850(yr, 7, 1) for yr in years]
        vb[:] = [[days_since_1850(yr, 1, 1), days_since_1850(yr + 1, 1, 1)] for yr in years]
    else:
        vt[:] = [days_since_1850(yr + 1, 1, 1) for yr in years]
    crs = ds.createVariable("crs", "i4")
    crs.grid_mapping_name = "polar_stereographic"; crs.latitude_of_projection_origin = -90.0
    crs.standard_parallel = -71.0; crs.straight_vertical_longitude_from_pole = 0.0
    crs.false_easting = 0.0; crs.false_northing = 0.0; crs.semi_major_axis = 6378137.0
    crs.inverse_flattening = 298.257223563; crs.epsg_code = "EPSG:3031"
    v = ds.createVariable(var, "f4", ("time", "y", "x"), zlib=True, complevel=4, fill_value=np.float32(FILL))
    _sn = standard_name(meta)
    if _sn:
        v.standard_name = _sn
    v.long_name = meta["long_name"]; v.units = meta["units"]
    v.grid_mapping = "crs"; v.coordinates = "y x"
    if is_flux:
        v.cell_methods = "time: mean"
    else:
        v.cell_methods = "time: point"
    _global_attrs(ds, meta)
    return ds, v


def write_slice(v, k, plane):
    r"""One year of one gridded variable, filled where the policy left NaN."""
    arr = np.asarray(plane, dtype="f4").copy()
    arr[~np.isfinite(arr)] = FILL
    v[k] = arr


def write_scalar(path, var, meta, years, values, is_flux):
    import netCDF4
    with netCDF4.Dataset(path, "w", format="NETCDF4") as ds:
        ds.createDimension("time", None)
        vt = ds.createVariable("time", "f4", ("time",))
        vt.units = TIME_UNITS; vt.calendar = "standard"; vt.standard_name = "time"; vt.long_name = "time"; vt.axis = "T"
        if is_flux:
            ds.createDimension("nv", 2); vt.bounds = "time_bnds"
            vb = ds.createVariable("time_bnds", "f4", ("time", "nv"))
            vt[:] = [days_since_1850(yr, 7, 1) for yr in years]
            vb[:] = [[days_since_1850(yr, 1, 1), days_since_1850(yr + 1, 1, 1)] for yr in years]
        else:
            vt[:] = [days_since_1850(yr + 1, 1, 1) for yr in years]
        v = ds.createVariable(var, "f4", ("time",), fill_value=FILL)
        arr = np.array(values, dtype="f4"); arr[~np.isfinite(arr)] = FILL
        v[:] = arr
        _sn = standard_name(meta)
        if _sn:
            v.standard_name = _sn
        v.long_name = meta["long_name"]; v.units = meta["units"]
        v.cell_methods = "time: mean" if is_flux else "time: point"
        _global_attrs(ds, meta)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("annual"); ap.add_argument("--out-dir", required=True)
    ap.add_argument("--esm", required=True); ap.add_argument("--scenario", required=True); ap.add_argument("--exp", required=True)
    ap.add_argument("--source-id", default=os.environ.get("ISMIP7_SOURCE_ID", "RICE"))
    ap.add_argument("--ism-id", default=os.environ.get("ISMIP7_ISM_ID", "icepack2"))
    ap.add_argument("--set-id", default="CORE")
    ap.add_argument("--contact-name", default=os.environ.get("ISMIP7_CONTACT_NAME", "Andrew Hoffman"))
    ap.add_argument("--contact-email", default=os.environ.get("ISMIP7_CONTACT_EMAIL", "ah301@rice.edu"))
    ap.add_argument("--scalars", default=None, help="the *_ismip7_scalars.csv (default: next to the annual file)")
    a = ap.parse_args()
    req = request_table()
    years = AnnualOutput.years_on_disk(a.annual)
    if not years:
        raise FileNotFoundError(
            f"no yearly checkpoints {os.path.basename(a.annual)[:-3]}_<year>.h5 "
            f"beside {a.annual}; pass the stem the run was named after."
        )
    if years != list(range(years[0], years[-1] + 1)):
        missing = sorted(set(range(years[0], years[-1] + 1)) - set(years))
        raise FileNotFoundError(
            f"the yearly checkpoints beside {a.annual} run {years[0]}-{years[-1]} "
            f"with {len(missing)} missing ({missing[0]}..{missing[-1]}); a "
            f"submission needs a contiguous series."
        )
    with fd.CheckpointFile(AnnualOutput.year_path(a.annual, years[0]), "r") as chk:
        mesh = chk.load_mesh()
    print(f"{os.path.basename(a.annual)}: years {years[0]}-{years[-1]}, "
          f"{len(VARIABLES_2D)} variables", flush=True)
    W = overlap_operator(mesh, a.annual + ".overlap.npz")
    print(f"  overlap operator {W.shape}, {W.nnz} entries, pixels covered {int((W.sum(axis=1) > 0).sum())}", flush=True)

    scal = a.scalars or a.annual.replace("_ismip7_annual.h5", "_ismip7_scalars.csv")
    if not os.path.exists(scal):
        raise FileNotFoundError(
            f"no scalars csv at {scal}; all ten scalar variables are "
            f"mandatory, so the submission cannot be written without it. "
            f"Pass --scalars if the run's csv is named differently."
        )
    with open(scal) as f:
        rows = {int(r["year"]): r for r in csv.DictReader(f)}
    missing_rows = [yr for yr in years if yr not in rows]
    if missing_rows:
        raise FileNotFoundError(
            f"{scal} has no row for year(s) {missing_rows[0]}"
            + (f"..{missing_rows[-1]}" if len(missing_rows) > 1 else "")
            + " although the gridded output has them; the scalar series "
              "would be submitted with a hole in it."
        )

    outdir = os.path.join(a.out_dir, "AIS", a.source_id, a.ism_id, a.set_id, a.exp)
    os.makedirs(outdir, exist_ok=True)
    tag = f"AIS_{a.source_id}_{a.ism_id}_m001_{a.esm}_f001_{a.scenario}_{a.exp}_{years[0]}-{years[-1]}"
    GLOBAL.update({"model": a.ism_id, "group": a.source_id, "crs": "EPSG:3031",
                   "contact_name": a.contact_name,
                   "contact_email": a.contact_email,
                   "source_id": a.source_id, "ism_id": a.ism_id, "experiment_id": a.exp,
                   "forcing": f"{a.esm} {a.scenario}",
                   "title": f"ISMIP7 AIS {a.exp} {a.esm} {a.scenario}, {a.source_id} {a.ism_id}"})

    # One open file per gridded variable, one year of grid in memory at a
    # time: the cubes for a 286-year submission on the production mesh would
    # otherwise be tens of GB in a serial script. Everything is written to
    # <name>.nc.tmp and renamed into place only once the last year and the
    # scalars are through, so a failure part-way leaves no half-filled file
    # that looks like a finished submission, and any previous files stand.
    paths = {var: os.path.join(outdir, f"{var}_{tag}.nc") for var in VARIABLES_2D}
    for var in SCALARS:
        paths[var] = os.path.join(outdir, f"{var}_{tag}.nc")
    tmps = {var: path + ".tmp" for var, path in paths.items()}
    handles = {}
    try:
        for var in VARIABLES_2D:
            handles[var] = create_2d(tmps[var], var, req[var], years,
                                     req[var]["Type"] == "FL")
        stats = {var: [np.inf, -np.inf, 0] for var in VARIABLES_2D}
        for k, yr in enumerate(years):
            with fd.CheckpointFile(AnnualOutput.year_path(a.annual, yr), "r") as chk:
                ymesh = chk.load_mesh()
                cells = {var: chk.load_function(ymesh, name=var).dat.data_ro.copy()
                         for var in VARIABLES_BANKED}
            # the overlap operator is built once from the first year's mesh;
            # a series whose links remeshed cannot be regridded through it
            if len(cells["lithk"]) != W.shape[1]:
                raise ValueError(
                    f"year {yr} has {len(cells['lithk'])} cells but the overlap "
                    f"operator was built from year {years[0]}'s {W.shape[1]}: "
                    f"the mesh changed inside the series, so one conservative "
                    f"operator cannot cover it."
                )
            masks = {"no_ice": cells["sftgif"] > 0.5,
                     "no_grounded_ice": cells["sftgrf"] > 0.5,
                     "no_floating_ice": cells["sftflf"] > 0.5}

            def plane(var):
                meta = req[var]; policy = meta["fill_policy"]
                m = masks[policy] if policy in masks else None
                return (regrid(W, cells[var], policy, m)
                        * CONVERT[meta["units"]]).reshape(ISMIP7_NY, ISMIP7_NX).astype("f4")

            # the checker requires orog == base + lithk pixel by pixel and
            # orog >= 0; the elevations are covered-part means while lithk is a
            # whole-pixel mean (diluted where the pixel is partly covered), so
            # the base absorbs the dilution: base := orog - lithk on the grid.
            # Rebuilding orog instead put negative surfaces on partly covered
            # floating pixels (base < 0 plus a diluted thickness).
            done = {"orog": plane("orog"), "lithk": plane("lithk")}
            done["base"] = done["orog"] - done["lithk"]
            for var in VARIABLES_2D:
                p_yr = done.get(var)
                if p_yr is None:
                    p_yr = plane(var)
                write_slice(handles[var][1], k, p_yr)
                finite = np.isfinite(p_yr)
                st = stats[var]
                if finite.any():
                    st[0] = min(st[0], float(np.min(p_yr[finite])))
                    st[1] = max(st[1], float(np.max(p_yr[finite])))
                st[2] += int(finite.sum())
        for var in VARIABLES_2D:
            handles.pop(var)[0].close()
            lo, hi, nfin = stats[var]
            pct = 100.0 * nfin / (len(years) * ISMIP7_NY * ISMIP7_NX)
            rng = f"[{lo:.3e}, {hi:.3e}]" if nfin else "[all fill]"
            print(f"  {var:12s} {req[var]['units']:11s} {rng} "
                  f"{pct:5.1f}% valid  -> {os.path.basename(paths[var])}", flush=True)

        for var in SCALARS:
            meta = req[var]
            values = [float(rows[yr][var]) for yr in years]
            write_scalar(tmps[var], var, meta, years, values, meta["Type"] == "FL")
        print(f"  {len(SCALARS)} scalars from {os.path.basename(scal)}")

        for var, tmp in tmps.items():
            os.replace(tmp, paths[var])
    except BaseException:
        for ds, _ in handles.values():
            ds.close()
        for tmp in tmps.values():
            if os.path.exists(tmp):
                os.remove(tmp)
        raise
    print(f"wrote {outdir}")


if __name__ == "__main__":
    main()
