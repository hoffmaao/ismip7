r"""ISMIP7 forcing data reader for Antarctic simulations."""

import os
import numpy as np

_SEC_PER_YEAR = 31556926.0
_RHO_ICE = 917.0
_RHO_WATER = 1000.0

# Constants for the Burgard et al. 2022 quadratic-mixed-slope melt
# parameterization, taken verbatim from multimelt.constants
# (https://github.com/ClimateClara/multimelt). The ISMIP7 ocean-forcing
# pipeline calibrates K under exactly this decomposition.
_RHO_SW = 1028.0       # seawater, kg/m^3
_RHO_I = 917.0         # ice,      kg/m^3
_C_PO = 3974.0         # seawater specific heat, J/(kg K)
_L_I = 3.34e5          # latent heat of fusion of ice, J/kg
_BETA_S = 7.86e-4      # haline contraction coefficient (Lazeroms), 1/PSU
_G = 9.81              # gravity, m/s^2
_F_CORIOLIS = 1.4e-4   # representative Antarctic Coriolis parameter, 1/s

# melt_factor = (rho_sw * c_po) / (rho_i * L_i)   [1/K]
_MELT_FACTOR = (_RHO_SW * _C_PO) / (_RHO_I * _L_I)

# K50 median from Burgard 2022 calibration
# (parameter_selection_quadratic_example.ipynb).
_K_DEFAULT = 11.5e-5

_DEFAULT_DATA_ROOT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "ISMIP7", "AIS",
)


def smb_kgm2s_to_myr(smb_kgm2s):
    r"""Convert SMB from kg/m^2/s to m/yr ice equivalent."""
    return smb_kgm2s * _SEC_PER_YEAR / _RHO_ICE * (_RHO_WATER / _RHO_ICE)


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
            times = ds["time"].values
            # Handle cftime calendars by matching on year
            if hasattr(times[0], "year"):
                yr = int(round(year))
                idx = min(range(len(times)), key=lambda i: abs(times[i].year - yr))
                tf_slice = tf_var.isel(time=idx)
            else:
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
                # z coords are negative (depth below sea level), draft is also negative
                draft_arr = xr.DataArray(np.asarray(draft), dims="node")
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

    def get_salinity(self, year, mesh_x, mesh_y, draft=None, fill=34.5):
        r"""Get ambient salinity at ice draft depth, interpolated to mesh."""
        import xarray as xr

        ds = self._load_variable("so")
        if ds is None:
            return np.full(len(mesh_x), fill)

        so_var = None
        for name in ds.data_vars:
            if name.lower() in ("so", "salinity"):
                so_var = ds[name]
                break
        if so_var is None:
            so_var = ds[list(ds.data_vars)[0]]

        if "time" in so_var.dims:
            times = ds["time"].values
            if hasattr(times[0], "year"):
                yr = int(round(year))
                idx = min(range(len(times)),
                          key=lambda i: abs(times[i].year - yr))
                so_slice = so_var.isel(time=idx)
            else:
                so_slice = so_var.sel(time=year, method="nearest")
        else:
            so_slice = so_var

        zdim = None
        for d in so_slice.dims:
            if d.lower() in ("z", "depth", "lev"):
                zdim = d
                break

        if zdim is not None:
            if draft is not None:
                draft_arr = xr.DataArray(np.asarray(draft), dims="node")
                so_slice = so_slice.interp({zdim: draft_arr}, method="nearest")
            else:
                so_slice = so_slice.isel({zdim: 0})

        mx = xr.DataArray(np.asarray(mesh_x), dims="node")
        my = xr.DataArray(np.asarray(mesh_y), dims="node")

        xdim = [d for d in so_slice.dims if d.lower() == "x"]
        ydim = [d for d in so_slice.dims if d.lower() == "y"]
        if xdim and ydim:
            vals = so_slice.interp({xdim[0]: mx, ydim[0]: my}, method="nearest")
        else:
            return np.full(len(mesh_x), fill)

        return np.nan_to_num(vals.values.flatten(), nan=fill)

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


def quadratic_mixed_slope(tf, salinity, sin_alpha, K=_K_DEFAULT):
    r"""ISMIP7 / Burgard et al. 2022 quadratic-mixed-slope local melt.

    Matches multimelt.melt_functions.quadratic_mixed_slope with TF_avg = TF
    (local-quadratic variant):

        m = K * melt_factor * U_factor * TF * |TF| * sin(alpha)

    where

        melt_factor = (rho_sw * c_po) / (rho_i * L_i)             [1/K]
        U_factor    = (c_po / L_i) * beta_S * g/(2|f|) * S0       [m/s/K]

    Inputs:
        tf        : thermal forcing T - T_f at ice base, K (numpy array)
        salinity  : ambient salinity at ice draft, PSU (numpy array)
        sin_alpha : sin of local ice-draft slope, dimensionless (numpy array)
        K         : dimensionless tuning factor (scalar or per-node array).
                    Burgard K50 = 1.15e-4.

    Returns melt rate in m/yr ice equivalent (positive = melting).
    """
    U_factor = (_C_PO / _L_I) * _BETA_S * (_G / (2.0 * abs(_F_CORIOLIS))) \
        * salinity
    melt = K * _MELT_FACTOR * U_factor * tf * np.abs(tf) * sin_alpha
    return melt * _SEC_PER_YEAR


def load_K_per_basin(npz_path, mesh_x, mesh_y, fill=0.0):
    r"""Load a per-basin calibrated K and return a per-node array.

    `npz_path` should be the output of `antarctica/scripts/calibrate_melt.py`
    and is expected to contain `basin_ids` (int) and `K_basin` (float) plus
    the IMBIE2 basin file path (the IMBIE2 8 km grid is re-read here so
    that the K-field can be remapped to *any* mesh, not just the one used
    during calibration).

    Returns an array of shape (len(mesh_x),) of per-node K values, with
    `fill` outside the calibrated basin set or where K_basin is NaN.
    """
    import xarray as xr
    from scipy.interpolate import RegularGridInterpolator

    data = np.load(npz_path)
    bids = np.asarray(data["basin_ids"]).astype(int)
    Kbas = np.asarray(data["K_basin"]).astype(float)

    imbie2 = os.path.join(
        os.environ.get("ISMIP7_DATA_ROOT",
                       os.path.join(os.path.dirname(
                           os.path.dirname(os.path.abspath(__file__))),
                           "ISMIP7", "AIS")),
        "parameterisations", "ocean", "imbie2",
        "basin_numbers_ismip8km_v2.nc",
    )
    ds = xr.open_dataset(imbie2)
    xa = ds["x"].values; ya = ds["y"].values
    bn = ds["basinNumber"].values
    if ya[0] > ya[-1]:
        ya = ya[::-1]; bn = bn[::-1, :]
    if xa[0] > xa[-1]:
        xa = xa[::-1]; bn = bn[:, ::-1]
    interp = RegularGridInterpolator(
        (ya, xa), bn.astype(np.float32),
        method="nearest", bounds_error=False, fill_value=-1.0,
    )
    pts = np.column_stack([np.asarray(mesh_y), np.asarray(mesh_x)])
    basin_node = np.round(interp(pts)).astype(int)
    ds.close()

    K_field = np.full(len(mesh_x), fill, dtype=float)
    for bid, kb in zip(bids, Kbas):
        if np.isfinite(kb):
            K_field[basin_node == bid] = kb
    return K_field


def compute_sin_alpha(ctx):
    r"""Return sin(alpha) of local ice-draft slope at CG1 nodes.

    Computes draft = s - h, projects grad(draft) into the vector CG1
    space, and returns sin(arctan(|grad|)) = |grad|/sqrt(1 + |grad|^2).
    """
    import firedrake as fd
    Q = ctx["Q"]
    V = ctx["V"]
    h = ctx["h"]
    s = ctx["s"]
    draft = fd.Function(Q).interpolate(s - h)
    grad_draft = fd.project(fd.grad(draft), V)
    g = grad_draft.dat.data_ro
    gmag = np.sqrt(g[:, 0] ** 2 + g[:, 1] ** 2)
    return gmag / np.sqrt(1.0 + gmag * gmag)


def make_forcing_callback(atm=None, ocean=None, fracture=None,
                          K=_K_DEFAULT, K_per_basin_npz=None,
                          smb_anomaly=True):
    r"""Build a forcing callback for use with simulation.run_simulation().

    Ocean melt uses the ISMIP7 Burgard quadratic_mixed_slope formula
    (local-quadratic variant, TF_avg = TF).

    K can be:
      - a scalar (single dimensionless K applied everywhere), or
      - a per-node numpy array (same length as mesh CG1 dofs), or
      - left at default while `K_per_basin_npz` points to the output of
        `calibrate_melt.py`; the per-basin K is then looked up on the mesh
        on first call and reused for subsequent steps.

    smb_anomaly selects between `acabf-anomaly` (True) and the full
    `acabf` field (False). The anomaly files are referenced to the ESM's
    1960-1989 climatology, so anomaly-only SMB is NOT a usable total —
    pass False (full field) unless you are adding the matching baseline
    to ctx["accum"] yourself.
    """
    K_field_cache = {"arr": None}

    def callback(ctx, t_yr):
        mesh_x = ctx["mesh"].coordinates.dat.data_ro[:, 0]
        mesh_y = ctx["mesh"].coordinates.dat.data_ro[:, 1]

        if atm is not None:
            smb = atm.get_smb(t_yr, mesh_x, mesh_y, anomaly=smb_anomaly)
            ctx["accum"].dat.data[:] = smb

        if ocean is not None and "ocean_melt" in ctx:
            h = ctx["h"].dat.data_ro
            b = ctx["b"].dat.data_ro
            s = ctx["s"].dat.data_ro
            # Ice shelf draft (negative depth below sea level)
            draft = np.minimum(s - h, 0.0)

            tf = ocean.get_thermal_forcing(t_yr, mesh_x, mesh_y, draft=draft)
            sal = ocean.get_salinity(t_yr, mesh_x, mesh_y, draft=draft)
            sin_alpha = compute_sin_alpha(ctx)

            # Resolve K: per-basin npz takes precedence if supplied.
            if K_per_basin_npz is not None:
                if K_field_cache["arr"] is None:
                    K_field_cache["arr"] = load_K_per_basin(
                        K_per_basin_npz, mesh_x, mesh_y, fill=0.0
                    )
                K_use = K_field_cache["arr"]
            else:
                K_use = K

            melt = quadratic_mixed_slope(tf, sal, sin_alpha, K=K_use)

            # Only apply melt where ice is floating (haf <= 0)
            haf = s - (b + (_RHO_WATER / _RHO_ICE) * np.maximum(-b, 0.0))
            floating = haf <= 0
            ctx["ocean_melt"].dat.data[:] = np.where(floating, melt, 0.0)

    return callback
