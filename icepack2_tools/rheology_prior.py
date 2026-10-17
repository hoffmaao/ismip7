r"""rheology_prior.py - a fluidity prior mean from a depth-averaged ice
temperature raster.

Recinos et al. (2023) take the prior mean of their rheology control from a
three-dimensional thermomechanical temperature field (Pattyn 2010), not from
a thermal model of their own. This module does the same with a gridded
depth-averaged temperature (``antarctica/data/temp/Pattyn_2013.tif``, EPSG:3031,
5 km, kelvin): the prior mean is ``A_prior = rate_factor(T)`` at every node of
the control space, with no strain heating, no geothermal flux and no water
content of ours on top. Andrew (26 Sep 2026): "use the depth averaged
temperature field from Pattyn".

Where the raster is undefined. The product is not NaN outside the modelled
ice: it carries an extrapolation halo (in-range values along the coast, then
-2400 K to 4500 K far out), so a range test alone does not find the holes.
A pixel is DEFINED when it is finite, within [200 K, 273.15 K] and, when a
BedMachine file is given, under BedMachine ice (mask 2, 3 or 4). Nodes whose
nearest pixel is undefined take a temperature supplied by the caller
(:func:`column_fill_temperature`: the mean-annual surface temperature on
grounded ice, its mean with the in-situ freezing point on floating ice), or
failing that the nearest defined pixel. Andrew (26 Sep): "fill the pattyn_2013
prior with temperatures where the mesh has undefined holes". The count of
ice-covered nodes that took a fill is reported. None of this caps the prior:
it replaces values that were never the model's.
"""
import numpy as np


T_MIN = 200.0          # K, colder than any ice on Earth
T_MAX = 273.15         # K, the melting point (no pressure dependence here)
ICE_MASKS = (2, 3, 4)  # BedMachine: grounded, floating, Lake Vostok


def clean_temperature(T, t_min=T_MIN, t_max=T_MAX, defined=None):
    r"""Replace undefined pixels by the nearest defined one.

    A pixel is undefined when non-finite, outside ``[t_min, t_max]`` or, with
    ``defined`` (a boolean array of the raster's shape), marked False there.
    Returns ``(T_clean, replaced)`` with ``replaced`` marking the filled
    pixels. Raises if no pixel is defined."""
    from scipy.ndimage import distance_transform_edt
    T = np.asarray(T, dtype=float)
    valid = np.isfinite(T) & (T >= t_min) & (T <= t_max)
    if defined is not None:
        defined = np.asarray(defined, dtype=bool)
        if defined.shape != T.shape:
            raise ValueError(f"defined mask shape {defined.shape} is not the "
                             f"raster's {T.shape}")
        valid &= defined
    if not valid.any():
        raise ValueError("temperature raster has no defined pixel in "
                         f"[{t_min:g}, {t_max:g}] K")
    replaced = ~valid
    if not replaced.any():
        return T.copy(), replaced
    _, (iy, ix) = distance_transform_edt(replaced, return_indices=True)
    return T[iy, ix], replaced


def read_temperature_raster(path):
    r"""``(x, y, T)`` of a single-band GeoTIFF: ``x`` ascending, ``y``
    ascending, ``T[j, i]`` at ``(y[j], x[i])`` (pixel centres)."""
    import rasterio
    with rasterio.open(path) as r:
        T = r.read(1).astype(float)
        if r.nodata is not None and np.isfinite(r.nodata):
            T[T == r.nodata] = np.nan
        dx, dy = r.res
        x = r.bounds.left + dx * (np.arange(r.width) + 0.5)
        y = r.bounds.top - dy * (np.arange(r.height) + 0.5)
    # rasterio rows run north to south; flip to ascending y
    return x, y[::-1], T[::-1, :]


def ice_extent_on_grid(bm_fn, x, y, stride=10):
    r"""Boolean ``(len(y), len(x))`` array: is the pixel centre under
    BedMachine ice (mask in :data:`ICE_MASKS`)? The BedMachine mask is read
    every ``stride`` pixels (500 m x 10 = 5 km, the temperature raster's
    resolution) and sampled nearest."""
    import netCDF4
    from scipy.interpolate import RegularGridInterpolator
    with netCDF4.Dataset(bm_fn, "r") as d:
        bx = np.asarray(d["x"][::stride], dtype=float)
        by = np.asarray(d["y"][::stride], dtype=float)
        m = np.asarray(d["mask"][::stride, ::stride])
    ice = np.isin(m, ICE_MASKS).astype(float)
    if by[0] > by[-1]:
        by, ice = by[::-1], ice[::-1, :]
    if bx[0] > bx[-1]:
        bx, ice = bx[::-1], ice[:, ::-1]
    # nearest sample everywhere, including the up-to-(stride-1) pixels the
    # strided read leaves past its last sample at the BedMachine edge
    f = RegularGridInterpolator((by, bx), ice, method="nearest",
                                bounds_error=False, fill_value=None)
    X, Y = np.meshgrid(x, y)
    return f(np.c_[Y.ravel(), X.ravel()]).reshape(len(y), len(x)) > 0.5


def column_fill_temperature(T_srf, H, b, Q, rho_i=917.0, rho_w=1024.0):
    r"""A depth-averaged temperature for nodes the raster does not define:
    the mean-annual surface temperature ``T_srf`` [K] where the ice is
    grounded (thin fringe ice is close to isothermal with its surface) and the
    mean of ``T_srf`` and the in-situ seawater freezing point at the draft
    where it floats (a linear column between the two). ``H``, ``b`` on ``Q``
    or liftable to it. Clipped to ``[T_MIN, T_MAX]``."""
    from firedrake import Constant, Function, conditional, max_value, min_value
    from icepack2_tools.thermo_model import ocean_melting_point
    haf = H - Constant(rho_w / rho_i) * max_value(-b, Constant(0.0))
    floating = conditional(haf <= 0.0, 1.0, 0.0)
    T_float = 0.5 * (T_srf + ocean_melting_point(H))
    T = (1.0 - floating) * T_srf + floating * T_float
    T = min_value(max_value(T, Constant(T_MIN)), Constant(T_MAX))
    return Function(Q, name="fill_temperature").interpolate(T)


def temperature_on_space(path, Q, t_min=T_MIN, t_max=T_MAX, defined=None,
                         fill=None):
    r"""Depth-averaged temperature [K] from ``path`` on the CG space ``Q``.

    Defined pixels are interpolated bilinearly (undefined pixels in a stencil
    hold their nearest defined value). A node whose bilinear stencil is more
    than half undefined takes ``fill`` (a Function on ``Q``) when given, else the
    interpolant of the nearest-filled raster.

    Returns ``(T, info)``; ``info`` has ``pixels_replaced`` / ``pixels_total``
    (raster), ``nodes_filled`` / ``nodes_total`` (this rank's owned nodes; the
    caller reduces), ``filled`` (the boolean node mask), ``fill_source`` and
    ``T_range`` over this rank's nodes."""
    from scipy.interpolate import RegularGridInterpolator
    from firedrake import Function, SpatialCoordinate, VectorFunctionSpace
    x, y, T_raw = read_temperature_raster(path)
    T_clean, replaced = clean_temperature(T_raw, t_min=t_min, t_max=t_max,
                                          defined=defined)
    mesh = Q.mesh()
    xy = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        SpatialCoordinate(mesh)).dat.data_ro
    pts = np.c_[xy[:, 1], xy[:, 0]]
    lin = RegularGridInterpolator((y, x), T_clean, method="linear",
                                  bounds_error=False, fill_value=None)
    # a node is undefined when the undefined pixels carry more than half of
    # its bilinear stencil, so a node beside the ice edge whose stencil is
    # mostly defined keeps the interpolant
    undefined = RegularGridInterpolator((y, x), replaced.astype(float),
                                        method="linear", bounds_error=False,
                                        fill_value=1.0)
    vals = np.clip(lin(pts), t_min, t_max)   # a node past the raster edge extrapolates
    filled = undefined(pts) > 0.5
    if fill is not None:
        vals = np.where(filled, np.asarray(fill.dat.data_ro, dtype=float), vals)
    T = Function(Q, name="temperature")
    T.dat.data[:] = vals
    info = {
        "pixels_replaced": int(replaced.sum()),
        "pixels_total": int(replaced.size),
        "nodes_filled": int(filled.sum()),
        "nodes_total": int(filled.size),
        "filled": filled,
        "fill_source": ("caller temperature" if fill is not None
                        else "nearest defined pixel"),
        "T_range": (float(vals.min()) if vals.size else np.inf,
                    float(vals.max()) if vals.size else -np.inf),
    }
    return T, info


def fluidity_prior_from_temperature_raster(path, Q, t_min=T_MIN, t_max=T_MAX,
                                           bm_fn=None, fill=None):
    r"""``A_prior = rate_factor(T)`` on ``Q`` from the temperature raster at
    ``path`` (icepack units, MPa^-3 yr^-1). With ``bm_fn`` (BedMachine), only
    pixels under its ice are defined; ``fill`` supplies the temperature of
    undefined nodes (:func:`column_fill_temperature`). Returns ``(A, info)``
    as :func:`temperature_on_space`, with ``info["temperature"]`` the field."""
    from firedrake import Function
    from icepack.models.viscosity import rate_factor
    defined = None
    if bm_fn is not None:
        x, y, _ = read_temperature_raster(path)
        defined = ice_extent_on_grid(bm_fn, x, y)
    T, info = temperature_on_space(path, Q, t_min=t_min, t_max=t_max,
                                   defined=defined, fill=fill)
    info["temperature"] = T
    A = Function(Q, name="fluidity_prior").interpolate(rate_factor(T))
    return A, info
