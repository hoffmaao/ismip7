r"""ISMIP7 forcing data reader for Antarctic simulations."""

import os
import numpy as np

_SEC_PER_YEAR = 31556926.0
_RHO_ICE = 917.0
_RHO_WATER = 1000.0

_DEFAULT_DATA_ROOT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "ISMIP7", "AIS",
)


def smb_kgm2s_to_myr(smb_kgm2s):
    r"""Convert SMB from kg/m^2/s to m/yr ice equivalent."""
    return smb_kgm2s * _SEC_PER_YEAR / _RHO_ICE * (_RHO_WATER / _RHO_ICE)


def load_racmo_smb_climatology(Q, clim_start=2000, clim_end=2029, data_dir=None,
                               target_res=8000.0, rho_ice=_RHO_ICE):
    r"""RACMO2.4p1 mean-annual SMB (m/yr ice equiv) as a Function on Q's mesh.

    The RACMO ANT11 grid is rotated-pole, so this reprojects the climatology to
    an intermediate EPSG:3031 raster with rasterio, then samples it onto the mesh
    with icepack.interpolate -- the same path used for BedMachine -- avoiding any
    scattered-point interpolation. ``smbgl`` is a monthly mass sum (kg/m^2), so the
    annual SMB is the sum of the 12 months, averaged over the climatology window.
    """
    import xarray as xr
    import pyproj
    import icepack
    from affine import Affine
    from rasterio.io import MemoryFile
    from rasterio.warp import reproject, Resampling

    if data_dir is None:
        base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        data_dir = os.path.join(base, "antarctica", "data", "racmo")
    fn = os.path.join(
        data_dir, "smbgl_monthlyS_ANT11_RACMO2.4p1_ERA5_197901_202312.nc"
    )
    ds = xr.open_dataset(fn)

    smb = ds["smbgl"].squeeze("height")  # (time, rlat, rlon), kg/m^2 per month
    yrs = smb["time"].dt.year.values
    sel = (yrs >= clim_start) & (yrs <= clim_end)
    n_years = len(np.unique(yrs[sel]))
    src = (smb.isel(time=sel).sum("time").values / n_years / rho_ice).astype("float64")

    rlon, rlat = ds["rlon"].values, ds["rlat"].values
    d = float(rlon[1] - rlon[0])
    src_crs = pyproj.CRS.from_cf(ds["rotated_pole"].attrs)
    src_transform = Affine.translation(rlon[0] - d / 2, rlat[0] - d / 2) * Affine.scale(d, d)
    ds.close()

    # Reproject onto the standard ISMIP AIS grid (EPSG:3031, centers +/- 3040 km).
    half = 3040000.0
    N = int(round(2 * half / target_res)) + 1
    dst_transform = (
        Affine.translation(-half - target_res / 2, half + target_res / 2)
        * Affine.scale(target_res, -target_res)
    )
    dst = np.full((N, N), np.nan, dtype="float64")
    reproject(
        src, dst,
        src_transform=src_transform, src_crs=src_crs,
        dst_transform=dst_transform, dst_crs="EPSG:3031",
        src_nodata=np.nan, dst_nodata=np.nan, resampling=Resampling.bilinear,
    )
    dst[~np.isfinite(dst)] = 0.0  # zero SMB over ocean / outside RACMO coverage

    with MemoryFile() as mf:
        with mf.open(
            driver="GTiff", height=N, width=N, count=1, dtype="float64",
            crs="EPSG:3031", transform=dst_transform,
        ) as out:
            out.write(dst, 1)
        with mf.open() as raster:
            return icepack.interpolate(raster, Q)


def _find_ismip7_data(data_root=None):
    if data_root is not None:
        return data_root
    env = os.environ.get("ISMIP7_DATA_ROOT")
    if env and os.path.isdir(env):
        return env
    if os.path.isdir(_DEFAULT_DATA_ROOT):
        return _DEFAULT_DATA_ROOT
    return None


def atmosphere_path(scenario, esm="CESM2-WACCM", variable="acabf-anomaly",
                    resolution="8000m", version="v2", data_root=None):
    root = _find_ismip7_data(data_root)
    if root is None:
        return None
    return os.path.join(root, esm, scenario, f"SDBN1-{resolution}",
                        variable, version)


def ocean_path(scenario, esm="CESM2-WACCM", variable="tf",
               version="v3", data_root=None):
    root = _find_ismip7_data(data_root)
    if root is None:
        return None
    return os.path.join(root, esm, scenario, "ocean", variable, version)


class ISMIP7Atmosphere:
    r"""Read ISMIP7 downscaled atmosphere forcing for Antarctica."""

    def __init__(self, data_root=None, esm="CESM2-WACCM", scenario="ssp585",
                 resolution="8000m", version="v2"):
        self.data_root = _find_ismip7_data(data_root)
        self.esm = esm
        self.scenario = scenario
        self.resolution = resolution
        self.version = version
        self._cache = {}
        self._grid_x = None
        self._grid_y = None

    def _var_dir(self, variable):
        return atmosphere_path(
            self.scenario, self.esm, variable,
            self.resolution, self.version, self.data_root,
        )

    def _load_year(self, variable, year):
        import xarray as xr

        vdir = self._var_dir(variable)
        if vdir is None or not os.path.isdir(vdir):
            return None

        pattern = f"{variable}_AIS_{self.esm}_{self.scenario}_SDBN1-{self.resolution}_{self.version}_{int(year)}.nc"
        path = os.path.join(vdir, pattern)

        if not os.path.exists(path):
            return None

        ds = xr.open_dataset(path)

        if self._grid_x is None:
            for xname in ["x", "X", "lon"]:
                if xname in ds.coords or xname in ds.dims:
                    self._grid_x = ds[xname].values
                    break
            for yname in ["y", "Y", "lat"]:
                if yname in ds.coords or yname in ds.dims:
                    self._grid_y = ds[yname].values
                    break

        data_vars = [v for v in ds.data_vars if v not in ("x", "y", "time")]
        if data_vars:
            da = ds[data_vars[0]]
            if "time" in da.dims:
                da = da.isel(time=0)
            result = da.load()
            ds.close()
            return result
        ds.close()
        return None

    def available_years(self, variable="acabf-anomaly"):
        r"""List available years for a variable."""
        vdir = self._var_dir(variable)
        if vdir is None or not os.path.isdir(vdir):
            return []
        years = []
        for fn in sorted(os.listdir(vdir)):
            if fn.endswith(".nc"):
                try:
                    yr = int(fn.rstrip(".nc").split("_")[-1])
                    years.append(yr)
                except ValueError:
                    pass
        return years

    def get_field(self, variable, year, mesh_x, mesh_y):
        r"""Get a forcing field interpolated to mesh coordinates."""
        import xarray as xr

        yr = int(round(year))
        da = self._load_year(variable, yr)
        if da is None:
            return np.zeros(len(mesh_x))

        mx = xr.DataArray(np.asarray(mesh_x), dims="node")
        my = xr.DataArray(np.asarray(mesh_y), dims="node")

        xdim = [d for d in da.dims if d.lower() == "x"]
        ydim = [d for d in da.dims if d.lower() == "y"]
        if xdim and ydim:
            vals = da.interp({xdim[0]: mx, ydim[0]: my}, method="nearest")
        else:
            dims = [d for d in da.dims if d not in ("time",)]
            if len(dims) >= 2:
                vals = da.interp({dims[-1]: mx, dims[-2]: my}, method="nearest")
            else:
                return np.zeros(len(mesh_x))

        return np.nan_to_num(vals.values.flatten(), nan=0.0)

    def get_smb(self, year, mesh_x, mesh_y, anomaly=True):
        r"""Get SMB field in m/yr ice equivalent."""
        var = "acabf-anomaly" if anomaly else "acabf"
        raw = self.get_field(var, year, mesh_x, mesh_y)
        return smb_kgm2s_to_myr(raw)

    def get_smb_gradient(self, year, mesh_x, mesh_y):
        r"""Get SMB elevation gradient (dacabfdz) for ice-elevation feedback."""
        return self.get_field("dacabfdz", year, mesh_x, mesh_y)

    def get_temperature(self, year, mesh_x, mesh_y, anomaly=True):
        r"""Get surface temperature (K or K anomaly)."""
        var = "ts-anomaly" if anomaly else "ts"
        return self.get_field(var, year, mesh_x, mesh_y)


class ISMIP7Ocean:
    r"""Read ISMIP7 ocean forcing for Antarctica."""

    def __init__(self, data_root=None, esm="CESM2-WACCM", scenario="ssp585",
                 version="v3"):
        self.data_root = _find_ismip7_data(data_root)
        self.esm = esm
        self.scenario = scenario
        self.version = version
        self._ds_cache = {}

    def _var_dir(self, variable):
        return ocean_path(
            self.scenario, self.esm, variable,
            self.version, self.data_root,
        )

    def _load_variable(self, variable):
        import xarray as xr

        if variable in self._ds_cache:
            return self._ds_cache[variable]

        vdir = self._var_dir(variable)
        if vdir is None or not os.path.isdir(vdir):
            return None

        nc_files = sorted(
            os.path.join(vdir, f) for f in os.listdir(vdir)
            if f.endswith(".nc")
        )
        if not nc_files:
            return None

        ds = xr.open_mfdataset(nc_files, combine="by_coords")
        self._ds_cache[variable] = ds
        return ds

    def get_thermal_forcing(self, year, mesh_x, mesh_y, draft=None):
        r"""Get thermal forcing at ice shelf base, interpolated to mesh."""
        import xarray as xr

        ds = self._load_variable("tf")
        if ds is None:
            return np.zeros(len(mesh_x))

        tf_var = None
        for name in ds.data_vars:
            if name.lower() in ("tf", "thermal_forcing", "thermalforcing"):
                tf_var = ds[name]
                break
        if tf_var is None:
            tf_var = ds[list(ds.data_vars)[0]]

        if "time" in tf_var.dims:
            tf_slice = tf_var.sel(time=year, method="nearest")
        else:
            tf_slice = tf_var

        zdim = None
        for d in tf_slice.dims:
            if d.lower() in ("z", "depth", "lev"):
                zdim = d
                break

        if zdim is not None:
            if draft is not None:
                draft_arr = xr.DataArray(np.abs(np.asarray(draft)), dims="node")
                tf_slice = tf_slice.interp({zdim: draft_arr}, method="nearest")
            else:
                tf_slice = tf_slice.isel({zdim: 0})

        mx = xr.DataArray(np.asarray(mesh_x), dims="node")
        my = xr.DataArray(np.asarray(mesh_y), dims="node")

        xdim = [d for d in tf_slice.dims if d.lower() == "x"]
        ydim = [d for d in tf_slice.dims if d.lower() == "y"]
        if xdim and ydim:
            vals = tf_slice.interp({xdim[0]: mx, ydim[0]: my}, method="nearest")
        else:
            return np.zeros(len(mesh_x))

        return np.nan_to_num(vals.values.flatten(), nan=0.0)

    def close(self):
        for ds in self._ds_cache.values():
            ds.close()
        self._ds_cache.clear()


class ISMIP7Fracture:
    r"""Read ISMIP7 fracture / ice shelf collapse forcing."""

    def __init__(self, data_root=None, esm="CESM2-WACCM", scenario="ssp585"):
        self.data_root = _find_ismip7_data(data_root)
        self.esm = esm
        self.scenario = scenario
        self._collapse_mask = None
        self._excess_melt = None

    def _fracture_dir(self):
        root = _find_ismip7_data(self.data_root)
        if root is None:
            return None
        return os.path.join(root, self.esm, self.scenario, "fracture")

    def load(self):
        import xarray as xr

        fdir = self._fracture_dir()
        if fdir is None or not os.path.isdir(fdir):
            return self

        for fn in os.listdir(fdir):
            path = os.path.join(fdir, fn)
            if "collapse_mask" in fn:
                self._collapse_mask = xr.open_dataset(path)
            elif "excess_melt" in fn:
                self._excess_melt = xr.open_dataset(path)

        return self

    def get_collapse_mask(self, year, mesh_x, mesh_y):
        r"""Get ice shelf collapse mask (0/1) at given year."""
        import xarray as xr

        if self._collapse_mask is None:
            return np.zeros(len(mesh_x))

        ds = self._collapse_mask
        var = list(ds.data_vars)[0]
        da = ds[var]

        if "time" in da.dims:
            da = da.sel(time=year, method="nearest")

        mx = xr.DataArray(np.asarray(mesh_x), dims="node")
        my = xr.DataArray(np.asarray(mesh_y), dims="node")

        xdim = [d for d in da.dims if d.lower() == "x"]
        ydim = [d for d in da.dims if d.lower() == "y"]
        if xdim and ydim:
            vals = da.interp({xdim[0]: mx, ydim[0]: my}, method="nearest")
        else:
            return np.zeros(len(mesh_x))

        return np.nan_to_num(vals.values.flatten(), nan=0.0)

    def close(self):
        if self._collapse_mask is not None:
            self._collapse_mask.close()
        if self._excess_melt is not None:
            self._excess_melt.close()


def make_forcing_callback(atm=None, ocean=None, fracture=None):
    r"""Build a forcing callback for use with simulation.run_simulation()."""
    def callback(ctx, t_yr):
        mesh_x = ctx["mesh"].coordinates.dat.data_ro[:, 0]
        mesh_y = ctx["mesh"].coordinates.dat.data_ro[:, 1]

        if atm is not None:
            smb = atm.get_smb(t_yr, mesh_x, mesh_y, anomaly=True)
            ctx["accum"].dat.data[:] = smb

    return callback
