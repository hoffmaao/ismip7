r"""Ice masks for the marine front of a buffered mesh (issue #167), and
the pixel classes the front is drawn from.

The _frontbm meshes and the ``vertex_front`` raster sampling take the front
from BedMachine itself (:func:`bedmachine_classes`): its ice edge is where
its thickness ends, so a mesh on that edge gives the front cells BedMachine's
own thickness with nothing to fill. The ISMIP7 observations MIPkit
(``AntarcticaObsISMIP7-v*.nc``) also carries ``icemask_greene``, a binary
ice / no-ice mask on a 500 m grid for 1997 and 2000 to 2021, dated 15 March;
:func:`icemask_tif` cuts one year into a GeoTIFF beside the dH/dt cache
(``<ISMIP7_OBS_DATA_ROOT>/icemask_cache``) for comparison
(front_mask_census.py, bm_front_flux.py --greene-tif). Its 2015 front lies
seaward of BedMachine's ice over 51,373 km2 where BedMachine holds no
thickness, up to 20 km out, which is why the front is BedMachine's.

:func:`classify` splits the pixels three ways: ice, marine (no ice, over
BedMachine ocean, floating ice or a bed below sea level) and land (the
rest). A marine front is an edge between ice and marine pixels.
"""

import datetime
import os
import re

import numpy as np

from .obs_dhdt import _obs_kit_path

MASK_VARIABLE = "icemask_greene"
TIME_VARIABLE = "greene_mask_time"

# Pixel classes of classify(); int8 so a continent-wide array is 150 MB.
LAND, ICE, MARINE = 0, 1, 2

# BedMachine mask values.
_BM_OCEAN, _BM_LAND, _BM_FLOATING = 0, 1, 3


def _cache_dir(cache_dir=None):
    if cache_dir is None:
        from .runconfig import obs_data_root
        cache_dir = os.path.join(obs_data_root(), "icemask_cache")
    return cache_dir


def _version(name):
    m = re.search(r"v(\d+)\.(\d+)", os.path.basename(name))
    return (int(m.group(1)), int(m.group(2))) if m else (0, 0)


def icemask_source(year, data_root=None, cache_dir=None):
    r"""``(kit, tif)``: the MIPkit the mask of ``year`` is cut from (None on a
    machine without one) and the GeoTIFF it is cached in. Without a kit the
    tif is the newest cached one for that year. Raises ``FileNotFoundError``
    when there is neither, or when ``ISMIP7_OBS_KIT`` names a missing file."""
    year = int(year)
    cache_dir = _cache_dir(cache_dir)
    try:
        kit = _obs_kit_path(data_root)
    except FileNotFoundError as no_kit:
        if os.environ.get("ISMIP7_OBS_KIT"):
            raise
        import glob
        cands = sorted(glob.glob(os.path.join(
            cache_dir, f"{MASK_VARIABLE}_{year}_AntarcticaObsISMIP7-v*.tif")))
        if not cands:
            raise FileNotFoundError(
                f"{no_kit} The cached mask would do in its place, and "
                f"{cache_dir} holds no {MASK_VARIABLE}_{year}_*.tif.") from no_kit
        return None, max(cands, key=_version)
    stem = os.path.basename(kit)[:-len(".nc")]
    return kit, os.path.join(cache_dir, f"{MASK_VARIABLE}_{year}_{stem}.tif")


def epoch_index(days_since_1900, year):
    r"""The index of the mask dated in ``year``: the 15 March one where a year
    has two (2000 also has 1 October)."""
    dates = [datetime.date(1900, 1, 1) + datetime.timedelta(days=int(d))
             for d in np.asarray(days_since_1900).ravel()]
    hits = [i for i, d in enumerate(dates) if d.year == int(year)]
    march = [i for i in hits if (dates[i].month, dates[i].day) == (3, 15)]
    if march:
        return march[0], dates[march[0]]
    if len(hits) == 1:
        return hits[0], dates[hits[0]]
    raise ValueError(
        f"{MASK_VARIABLE} has {len(hits)} masks dated {year}: "
        f"{[str(dates[i]) for i in hits]}")


def icemask_tif(year=2015, data_root=None, cache_dir=None):
    r"""Path to the int8 GeoTIFF (1 ice, 0 no ice, EPSG:3031) of the Greene
    mask dated ``year``, cutting it from the kit on first use. The write is
    atomic, as the dH/dt cache's is, so concurrent ranks or sites at worst
    write the same file twice."""
    import uuid

    kit, tif = icemask_source(year, data_root, cache_dir)
    if kit is None or os.path.exists(tif):
        return tif
    import netCDF4
    import rasterio
    from rasterio.transform import from_origin

    with netCDF4.Dataset(kit) as d:
        idx, date = epoch_index(d.variables[TIME_VARIABLE][:], year)
        x = np.asarray(d.variables["x"][:], dtype="f8")
        y = np.asarray(d.variables["y"][:], dtype="f8")
        var = d.variables[MASK_VARIABLE]
        var.set_auto_mask(False)
        arr = np.asarray(var[idx], dtype="i1")
    arr = np.where(arr == 1, 1, 0).astype("i1")
    if y[1] > y[0]:
        arr = arr[::-1, :]
        y = y[::-1]
    dx = abs(float(x[1] - x[0]))
    dy = abs(float(y[0] - y[1]))
    prof = dict(driver="GTiff", height=arr.shape[0], width=arr.shape[1],
                count=1, dtype="int8", crs="EPSG:3031",
                transform=from_origin(float(x[0]) - dx / 2.0,
                                      float(y[0]) + dy / 2.0, dx, dy),
                compress="deflate", tiled=True)
    os.makedirs(os.path.dirname(tif), exist_ok=True)
    tmp = f"{tif}.{os.getpid()}.{uuid.uuid4().hex}.part"
    try:
        with rasterio.open(tmp, "w", **prof) as dst:
            dst.write(arr, 1)
            dst.update_tags(source=os.path.basename(kit), variable=MASK_VARIABLE,
                            date=str(date), time_index=str(idx))
        os.replace(tmp, tif)
    except BaseException:
        try:
            os.remove(tmp)
        except OSError:
            pass
        raise
    return tif


def classify(ice, bm_mask, bed):
    r"""Pixel classes (:data:`LAND`, :data:`ICE`, :data:`MARINE`) from the
    Greene ice flag and BedMachine's mask and bed on the same pixels. A no-ice
    pixel is marine over BedMachine ocean or floating ice, or over any ice
    BedMachine grounds on a bed below sea level; over ice-free land, or ice
    grounded above sea level, it is land."""
    ice = np.asarray(ice).astype(bool)
    bm_mask = np.asarray(bm_mask)
    bed = np.asarray(bed)
    out = np.full(ice.shape, LAND, dtype="i1")
    out[ice] = ICE
    marine = (~ice & (bm_mask != _BM_LAND)
              & ((bm_mask == _BM_OCEAN) | (bm_mask == _BM_FLOATING) | (bed < 0.0)))
    out[marine] = MARINE
    return out


def _aligned_window(src, bounds):
    r"""The window of ``src`` covering ``bounds`` exactly, refusing a raster
    whose pixels do not line up with the mask's."""
    from rasterio.windows import from_bounds
    win = from_bounds(*bounds, transform=src.transform)
    off = (win.col_off, win.row_off, win.width, win.height)
    if any(abs(v - round(v)) > 1e-6 for v in off):
        raise ValueError(
            f"{src.name} is not on the ice mask's pixel grid (window {off})")
    return win.round_offsets().round_lengths()


def read_classes(mask_tif, bm_fn, bounds=None):
    r"""``(classes, transform)`` on the mask's grid, over ``bounds`` (left,
    bottom, right, top; the whole mask when None, which must then be whole
    pixels). BedMachine is read on the same pixels."""
    import rasterio
    from rasterio.windows import from_bounds

    with rasterio.open(mask_tif) as src:
        if bounds is None:
            bounds = tuple(src.bounds)
        win = from_bounds(*bounds, transform=src.transform)
        win = win.round_offsets(op="floor").round_lengths(op="ceil")
        ice = src.read(1, window=win) == 1
        transform = src.window_transform(win)
        snapped = rasterio.windows.bounds(win, src.transform)
    with rasterio.open(f"netcdf:{bm_fn}:mask") as src:
        bm_mask = src.read(1, window=_aligned_window(src, snapped))
    with rasterio.open(f"netcdf:{bm_fn}:bed") as src:
        bed = src.read(1, window=_aligned_window(src, snapped))
    if bm_mask.shape != ice.shape or bed.shape != ice.shape:
        raise ValueError(
            f"BedMachine window {bm_mask.shape} does not match the mask's "
            f"{ice.shape} over {snapped}")
    return classify(ice, bm_mask, bed), transform


# BedMachine mask values that hold ice: grounded, floating, Lake Vostok.
_BM_ICE = (2, 3, 4)


def bedmachine_classes(bm_fn, bounds=None):
    r"""``(classes, transform)`` of BedMachine's own ice (mask grounded,
    floating or lake) on its grid, over ``bounds`` (left, bottom, right, top;
    the whole raster when None)."""
    import rasterio
    from rasterio.windows import from_bounds

    with rasterio.open(f"netcdf:{bm_fn}:mask") as src:
        if bounds is None:
            win = rasterio.windows.Window(0, 0, src.width, src.height)
        else:
            win = from_bounds(*bounds, transform=src.transform)
            win = win.round_offsets(op="floor").round_lengths(op="ceil")
        bm_mask = src.read(1, window=win)
        transform = src.window_transform(win)
        snapped = rasterio.windows.bounds(win, src.transform)
    with rasterio.open(f"netcdf:{bm_fn}:bed") as src:
        bed = src.read(1, window=_aligned_window(src, snapped))
    return classify(np.isin(bm_mask, _BM_ICE), bm_mask, bed), transform

