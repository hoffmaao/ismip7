r"""What the forward records every year: the annual file survives a chained
resume, and ``limnsw`` is the mass above flotation.

``nots_projection.sbatch`` self-chains, so ``run_simulation`` is entered once
per Slurm link and rebuilds :class:`AnnualOutput` against the file the
previous link left. The years attribute of that file is what the writer turns
into the submitted time axis, so a link that truncated it would silently
submit only its own segment.

``limnsw`` is the request's "mass above floatation ... volume times density":
the integrand is the THICKNESS above flotation, which on a marine bed is
larger than the height above flotation by 1/rho_ratio and on a dry bed is
just the thickness.
"""
import numpy as np
import pytest

fd = pytest.importorskip("firedrake")

from icepack2_tools.ismip7_output import AnnualOutput, RHO_I  # noqa: E402

RHO_RATIO = 917.0 / 1024.0


@pytest.fixture
def two_cells():
    r"""Two DG0 cells of area 0.5 each on the unit square."""
    mesh = fd.UnitSquareMesh(1, 1)
    Q_dg = fd.FunctionSpace(mesh, "DG", 0)
    V = fd.VectorFunctionSpace(mesh, "CG", 1)
    return mesh, Q_dg, V


def _dg(Q, values):
    f = fd.Function(Q)
    f.dat.data[:] = values
    return f


def _write_year(annual, Q, V, h, bed):
    r"""One year end with a still ice sheet: h and bed per cell, no flow."""
    h_dg = _dg(Q, h)
    b = _dg(Q, bed)
    s = _dg(Q, np.asarray(bed) + np.asarray(h))
    u = fd.Function(V)
    tau = fd.Function(V)
    annual.year_end(h_dg, s, b, u, tau,
                    np.ones(len(h), dtype=bool), np.ones(len(h), dtype=bool))


def test_limnsw_is_the_thickness_above_flotation(two_cells, tmp_path):
    mesh, Q, V = two_cells
    h = [2000.0, 2000.0]
    bed = [-1000.0, 1000.0]          # one marine cell, one above sea level
    annual = AnnualOutput(mesh, Q, V, str(tmp_path / "out" / "annual.h5"),
                          str(tmp_path / "out" / "scalars.csv"),
                          first_year=2015, rho_ratio=RHO_RATIO)
    annual.start_year(_dg(Q, h))
    _write_year(annual, Q, V, h, bed)
    annual.close()

    import csv
    with open(tmp_path / "out" / "scalars.csv") as f:
        row = next(iter(csv.DictReader(f)))
    # marine cell: 2000 - 1000 / rho_ratio; dry cell: the whole thickness
    expected = RHO_I * 0.5 * ((2000.0 - 1000.0 / RHO_RATIO) + 2000.0)
    assert float(row["limnsw"]) == pytest.approx(expected, rel=1e-9)
    # height above flotation would have given rho_ratio times less on the
    # marine cell and bed + rho_ratio * h on the dry one
    haf_height = RHO_I * 0.5 * ((2000.0 - (1.0 - RHO_RATIO) * 2000.0 - 1000.0)
                                + (1000.0 + RHO_RATIO * 2000.0))
    assert float(row["limnsw"]) != pytest.approx(haf_height, rel=1e-3)


def test_a_resume_appends_to_the_annual_file(two_cells, tmp_path):
    mesh, Q, V = two_cells
    h = [1500.0, 1500.0]
    bed = [-500.0, -500.0]
    out = str(tmp_path / "out" / "annual.h5")
    scalars = str(tmp_path / "out" / "scalars.csv")

    first = AnnualOutput(mesh, Q, V, out, scalars, first_year=2015, rho_ratio=RHO_RATIO)
    first.start_year(_dg(Q, h))
    _write_year(first, Q, V, h, bed)     # year 2015
    _write_year(first, Q, V, h, bed)     # year 2016
    first.close()

    # the next link of the chain resumes from the checkpoint at 2017
    second = AnnualOutput(mesh, Q, V, out, scalars, first_year=2017, rho_ratio=RHO_RATIO)
    assert second.year == 2017
    second.start_year(_dg(Q, h))
    _write_year(second, Q, V, h, bed)    # year 2017
    second.close()

    with fd.CheckpointFile(out, "r") as chk:
        years = [int(y) for y in chk.get_attr("/", "years").split(",")]
        loaded = chk.load_mesh()
        thickness = [chk.load_function(loaded, name="lithk", idx=k).dat.data_ro.copy()
                     for k in range(len(years))]
    assert years == [2015, 2016, 2017]
    assert all(np.allclose(t, h) for t in thickness)

    import csv
    with open(scalars) as f:
        assert [int(r["year"]) for r in csv.DictReader(f)] == [2015, 2016, 2017]


def test_a_resume_off_a_year_boundary_is_refused(two_cells, tmp_path):
    mesh, Q, V = two_cells
    with pytest.raises(ValueError, match="year boundary"):
        AnnualOutput(mesh, Q, V, str(tmp_path / "out" / "annual.h5"),
                     str(tmp_path / "out" / "scalars.csv"),
                     first_year=2047.6, rho_ratio=RHO_RATIO)


def test_a_resume_that_would_leave_a_gap_is_refused(two_cells, tmp_path):
    mesh, Q, V = two_cells
    h = [1500.0, 1500.0]
    out = str(tmp_path / "out" / "annual.h5")
    scalars = str(tmp_path / "out" / "scalars.csv")
    first = AnnualOutput(mesh, Q, V, out, scalars, first_year=2015, rho_ratio=RHO_RATIO)
    first.start_year(_dg(Q, h))
    _write_year(first, Q, V, h, [-500.0, -500.0])
    first.close()
    with pytest.raises(ValueError, match="gap in the submitted series"):
        AnnualOutput(mesh, Q, V, out, scalars, first_year=2030, rho_ratio=RHO_RATIO)
