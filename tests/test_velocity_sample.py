r"""The area-mean velocity target (``ISMIP7_VELOCITY_SAMPLE=cell_mean``).

``geometry.raster_patch_mean`` gives each CG1 node the area-weighted mean of
the raster's cell means over the cells around it, counting only cells that
hold a valid pixel. Serial, a 20 km square of 1 km cells on a 500 m raster.
"""

# dual_friction first: it pulls icepack2 -> irksome, which must be imported
# before any UFL form is assembled.
import icepack2_tools.dual_friction  # noqa: F401,E402

import numpy as np                                              # noqa: E402
import netCDF4                                                  # noqa: E402
import pytest                                                   # noqa: E402
import rasterio                                                 # noqa: E402

from firedrake import FunctionSpace, RectangleMesh              # noqa: E402

from icepack2_tools import geometry as G                        # noqa: E402
from icepack2_tools.handoff import objective_mismatches         # noqa: E402
from icepack2_tools.relaxation import relaxation_attrs          # noqa: E402
from icepack2_tools.runconfig import velocity_sample            # noqa: E402

PX = 500.0
X = np.arange(-2000.0, 22000.0 + 1.0, PX)
Y = X[::-1].copy()


def _raster(path, values):
    with netCDF4.Dataset(path, "w") as d:
        d.createDimension("x", len(X))
        d.createDimension("y", len(Y))
        for name, vals in (("x", X), ("y", Y)):
            v = d.createVariable(name, "f8", (name,))
            v[:] = vals
            v.standard_name = f"projection_{name}_coordinate"
        v = d.createVariable("v", "f4", ("y", "x"), fill_value=np.float32(-9999.0))
        v[:] = np.ma.masked_invalid(values.astype("f4"))
    return rasterio.open(f"netcdf:{path}:v")


def _spaces():
    mesh = RectangleMesh(20, 20, 20e3, 20e3)
    return FunctionSpace(mesh, "CG", 1), FunctionSpace(mesh, "DG", 0)


def test_the_patch_mean_of_a_ramp_sits_near_the_node(tmp_path):
    XX, _ = np.meshgrid(X, Y)
    Q, Q0 = _spaces()
    out, observed = G.raster_patch_mean(_raster(tmp_path / "ramp.nc", XX), Q0, Q)
    x = Q.mesh().coordinates.dat.data_ro[:, 0]
    assert observed.all()
    # interior nodes: the patch of cells is symmetric about the node
    interior = (x > 1500.0) & (x < 18500.0)
    assert np.abs(out.dat.data_ro[interior] - x[interior]).max() < PX


def test_nodes_with_no_valid_pixel_around_them_are_unobserved(tmp_path):
    XX, _ = np.meshgrid(X, Y)
    vals = np.where(XX < 10000.0, 100.0, np.nan)
    Q, Q0 = _spaces()
    out, observed = G.raster_patch_mean(_raster(tmp_path / "half.nc", vals), Q0, Q)
    x = Q.mesh().coordinates.dat.data_ro[:, 0]
    assert observed[x < 9500.0].all()
    assert not observed[x > 11500.0].any()
    assert np.all(out.dat.data_ro[~observed] == 0.0)
    # a node on the data edge averages only its observed cells, so no blend
    # with the zero that the unobserved ones would otherwise contribute
    assert np.allclose(out.dat.data_ro[observed], 100.0)


def test_the_velocity_sample_defaults_to_points_and_refuses_others(monkeypatch):
    monkeypatch.delenv("ISMIP7_VELOCITY_SAMPLE", raising=False)
    assert velocity_sample() == "point"
    monkeypatch.setenv("ISMIP7_VELOCITY_SAMPLE", "Cell_Mean")
    assert velocity_sample() == "cell_mean"
    monkeypatch.setenv("ISMIP7_VELOCITY_SAMPLE", "bilinear")
    with pytest.raises(ValueError):
        velocity_sample()


def _relaxed(**source):
    return relaxation_attrs(
        source_map="/maps/inversion_icepack2_rc_n3_dg0_2000.h5", source_sha256="ab12",
        source_attrs={"misfit_norm": "logvel", **source}, t_start=2014.0,
        t_end=2015.0, dt=0.0125, backdate_years=1.0, forcing="ocx protocol")


def test_a_relaxed_end_state_keeps_the_maps_velocity_sample():
    end = _relaxed(velocity_sample="cell_mean")
    assert end["velocity_sample"] == "cell_mean"
    # a re-inversion from it at the default sampling is a different objective
    assert objective_mismatches(end, {"misfit_norm": "logvel",
                                      "velocity_sample": "point"}) == \
        ["velocity_sample: 'cell_mean' -> 'point'"]
    assert objective_mismatches(end, {"misfit_norm": "logvel",
                                      "velocity_sample": "cell_mean"}) == []


def test_a_map_without_the_record_still_warm_starts():
    end = _relaxed()
    assert "velocity_sample" not in end
    for run in ("point", "cell_mean"):
        assert objective_mismatches(end, {"misfit_norm": "logvel",
                                          "velocity_sample": run}) == []
