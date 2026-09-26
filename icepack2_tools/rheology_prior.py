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

Raster cleaning, which is not a cap on the prior: outside the ice the product
holds fill values (-2400 K to 4500 K over the ocean), and a handful of ice
pixels (about 1 in 10,000) fall outside [200 K, 273.15 K]. Those pixels are
replaced by their nearest valid neighbour before interpolation, and the
number of nodes that took a replaced value is reported so a run can see how
much of the ice it touched.
"""
import numpy as np


T_MIN = 200.0          # K, colder than any ice on Earth
T_MAX = 273.15         # K, the melting point (no pressure dependence here)


def clean_temperature(T, t_min=T_MIN, t_max=T_MAX):
    r"""Replace non-finite and out-of-range pixels by the nearest valid one.

    Returns ``(T_clean, replaced)`` with ``replaced`` a boolean array marking
    the pixels that were filled. Raises if no pixel is valid."""
    from scipy.ndimage import distance_transform_edt
    T = np.asarray(T, dtype=float)
    valid = np.isfinite(T) & (T >= t_min) & (T <= t_max)
    if not valid.any():
        raise ValueError("temperature raster has no pixel in "
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


def temperature_on_space(path, Q, t_min=T_MIN, t_max=T_MAX):
    r"""Depth-averaged temperature [K] from ``path`` on the CG space ``Q``
    (bilinear in the raster), with the cleaning above.

    Returns ``(T, info)``; ``info`` has ``pixels_replaced`` (in the raster),
    ``nodes_from_replaced`` (nodes of ``Q`` on this rank whose value came from
    a replaced pixel, nearest-pixel test), and ``T_range`` over this rank's
    nodes. Callers reduce the counts over ranks."""
    from scipy.interpolate import RegularGridInterpolator
    from firedrake import Function, SpatialCoordinate, VectorFunctionSpace
    x, y, T_raw = read_temperature_raster(path)
    T_clean, replaced = clean_temperature(T_raw, t_min=t_min, t_max=t_max)
    mesh = Q.mesh()
    xy = Function(VectorFunctionSpace(mesh, "CG", 1)).interpolate(
        SpatialCoordinate(mesh)).dat.data_ro
    pts = np.c_[xy[:, 1], xy[:, 0]]
    lin = RegularGridInterpolator((y, x), T_clean, method="linear",
                                  bounds_error=False, fill_value=None)
    near = RegularGridInterpolator((y, x), replaced.astype(float),
                                   method="nearest", bounds_error=False,
                                   fill_value=1.0)
    T = Function(Q, name="temperature")
    vals = lin(pts)
    # a node just outside the raster extrapolates; keep it physical
    T.dat.data[:] = np.clip(vals, t_min, t_max)
    from_replaced = near(pts) > 0.5
    info = {
        "pixels_replaced": int(replaced.sum()),
        "pixels_total": int(replaced.size),
        "nodes_from_replaced": int(from_replaced.sum()),
        "nodes_total": int(from_replaced.size),
        "T_range": (float(vals.min()) if vals.size else np.inf,
                    float(vals.max()) if vals.size else -np.inf),
    }
    return T, info


def fluidity_prior_from_temperature_raster(path, Q, t_min=T_MIN, t_max=T_MAX):
    r"""``A_prior = rate_factor(T)`` on ``Q`` from the temperature raster at
    ``path`` (icepack units, MPa^-3 yr^-1). Returns ``(A_prior, info)`` as
    :func:`temperature_on_space` does."""
    from firedrake import Function
    from icepack.models.viscosity import rate_factor
    T, info = temperature_on_space(path, Q, t_min=t_min, t_max=t_max)
    A = Function(Q, name="fluidity_prior").interpolate(rate_factor(T))
    return A, info
