r"""What the forward records every year: the annual series survives a chained
resume, and ``limnsw`` is the mass above flotation.

``projection.sbatch`` self-chains, so ``run_simulation`` is entered once
per Slurm link and rebuilds :class:`AnnualOutput` against the years the
previous link left on disk. Those per-year files are what the writer turns
into the submitted time axis, so a link that truncated or renamed them would
silently submit the wrong series.

``limnsw`` is the request's "mass above floatation ... volume times density":
the integrand is the THICKNESS above flotation, which on a marine bed is
larger than the height above flotation by 1/rho_ratio and on a dry bed is
just the thickness.
"""
import numpy as np
import pytest

fd = pytest.importorskip("firedrake")

from icepack2_tools.ismip7_output import (AnnualOutput, FRONT_MELT, FRONT_MELT_ATTR,  # noqa: E402
                                          RHO_I, SECONDS_PER_YEAR, net_reference,
                                          net_snowfall, split_melt)

RHO_RATIO = 917.0 / 1024.0
# The scalars integrate over true area (issue #97). The unit meshes below sit
# within metres of the pole, where af2 = (1/k)^2 of EPSG:3031 is 1/k0^2.
AF2_POLE = 1.0567701662990163


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
    annual = AnnualOutput(mesh, Q, str(tmp_path / "out" / "annual.h5"),
                          str(tmp_path / "out" / "scalars.csv"),
                          first_year=2015, rho_ratio=RHO_RATIO)
    annual.start_year(_dg(Q, h))
    _write_year(annual, Q, V, h, bed)
    annual.close()

    import csv
    with open(tmp_path / "out" / "scalars.csv") as f:
        row = next(iter(csv.DictReader(f)))
    # marine cell: 2000 - 1000 / rho_ratio; dry cell: the whole thickness
    expected = RHO_I * 0.5 * AF2_POLE * ((2000.0 - 1000.0 / RHO_RATIO) + 2000.0)
    assert float(row["limnsw"]) == pytest.approx(expected, rel=1e-6)   # the csv carries 7 digits
    # height above flotation would have given rho_ratio times less on the
    # marine cell and bed + rho_ratio * h on the dry one
    haf_height = RHO_I * 0.5 * AF2_POLE * ((2000.0 - (1.0 - RHO_RATIO) * 2000.0 - 1000.0)
                                           + (1000.0 + RHO_RATIO * 2000.0))
    assert float(row["limnsw"]) != pytest.approx(haf_height, rel=1e-3)


def test_the_scalars_integrate_over_true_area(tmp_path):
    r"""ismip7-scalars weights every pixel by af2 = (1/k)^2 of EPSG:3031, so
    the native scalars weight every cell by it too (issue #97): the cell's
    map-plane area times af2 at its centroid. Two cells of 1e12 m2 whose
    centroids sit 687 and 1344 km from the pole, where af2 is 1.050 and
    1.033; the mass budget's own cell areas stay map-plane."""
    from icepack2_tools.regrid import area_factor
    mesh = fd.RectangleMesh(1, 1, 2.0e6, 0.5e6, originX=0.0, originY=-0.5e6)
    Q = fd.FunctionSpace(mesh, "DG", 0)
    V = fd.VectorFunctionSpace(mesh, "CG", 1)
    h, bed = [1000.0, 3000.0], [500.0, 500.0]        # both grounded on dry beds
    annual = AnnualOutput(mesh, Q, str(tmp_path / "out" / "annual.h5"),
                          str(tmp_path / "out" / "scalars.csv"),
                          first_year=2015, rho_ratio=RHO_RATIO)
    annual.start_year(_dg(Q, h))
    _write_year(annual, Q, V, h, bed)
    annual.close()
    import csv
    with open(tmp_path / "out" / "scalars.csv") as f:
        row = next(iter(csv.DictReader(f)))
    area = fd.assemble(fd.TestFunction(Q) * fd.dx).dat.data_ro
    X = fd.SpatialCoordinate(mesh)
    af2 = area_factor(*(fd.Function(Q).interpolate(X[i]).dat.data_ro for i in (0, 1)))
    assert np.allclose(annual.cell_area, area)
    assert float(row["iareagr"]) == pytest.approx(np.sum(area * af2), rel=1e-6)
    assert float(row["lim"]) == pytest.approx(RHO_I * np.sum(np.array(h) * area * af2), rel=1e-6)
    assert float(row["limnsw"]) == pytest.approx(float(row["lim"]), rel=1e-6)
    assert float(row["iareagr"]) / np.sum(area) == pytest.approx(1.0414669, rel=1e-6)


def test_a_resume_appends_to_the_annual_file(two_cells, tmp_path):
    mesh, Q, V = two_cells
    h = [1500.0, 1500.0]
    bed = [-500.0, -500.0]
    out = str(tmp_path / "out" / "annual.h5")
    scalars = str(tmp_path / "out" / "scalars.csv")

    first = AnnualOutput(mesh, Q, out, scalars, first_year=2015, rho_ratio=RHO_RATIO)
    first.start_year(_dg(Q, h))
    _write_year(first, Q, V, h, bed)     # year 2015
    _write_year(first, Q, V, h, bed)     # year 2016
    first.close()

    # the next link of the chain resumes from the checkpoint at 2017
    second = AnnualOutput(mesh, Q, out, scalars, first_year=2017, rho_ratio=RHO_RATIO)
    assert second.year == 2017
    second.start_year(_dg(Q, h))
    _write_year(second, Q, V, h, bed)    # year 2017
    second.close()

    assert AnnualOutput.years_on_disk(out) == [2015, 2016, 2017]
    for yr in (2015, 2016, 2017):
        with fd.CheckpointFile(AnnualOutput.year_path(out, yr), "r") as chk:
            loaded = chk.load_mesh()
            assert np.allclose(chk.load_function(loaded, name="lithk").dat.data_ro, h)

    import csv
    with open(scalars) as f:
        assert [int(r["year"]) for r in csv.DictReader(f)] == [2015, 2016, 2017]


def test_a_midyear_resume_carries_the_partial_year(two_cells, tmp_path):
    r"""The wall-clock budget stops a link between steps, so with the default
    dt=0.1 a chained resume lands mid-year. The flux means of that year must
    survive the link boundary and stay labelled with their own year."""
    mesh, Q, V = two_cells
    h = [1500.0, 1500.0]
    bed = [-500.0, -500.0]
    out = str(tmp_path / "out" / "annual.h5")
    scalars = str(tmp_path / "out" / "scalars.csv")

    first = AnnualOutput(mesh, Q, out, scalars, first_year=2015, rho_ratio=RHO_RATIO)
    first.start_year(_dg(Q, h))
    _write_year(first, Q, V, h, bed)            # year 2015 written, now inside 2016
    first.year_acc["acabf"][:] = [0.3, 0.3]     # 0.4 yr of accumulated sources
    first.year_time = 0.4
    state = (first.state_fields(), first.state_attrs())
    first.close()
    assert state[1]["ismip7_year"] == 2016
    assert state[1]["ismip7_year_time"] == pytest.approx(0.4)

    resume = {
        "year": state[1]["ismip7_year"],
        "year_time": state[1]["ismip7_year_time"],
        "series": state[1][AnnualOutput.STATE_SERIES],
        "acc": {k: state[0][AnnualOutput.STATE_PREFIX + k].dat.data_ro.copy()
                for k in AnnualOutput.ACCUMULATORS},
        "h_year_start": state[0][AnnualOutput.STATE_THICKNESS].dat.data_ro.copy(),
    }
    # the next link picks up at t=2016.4, inside the year already in progress
    second = AnnualOutput(mesh, Q, out, scalars, first_year=2016.4,
                          rho_ratio=RHO_RATIO, resume=resume)
    assert second.year == 2016
    assert second.year_time == pytest.approx(0.4)
    assert np.allclose(second.year_acc["acabf"], [0.3, 0.3])
    assert second.h_year_start is not None       # the year is in progress, not restarted

    second.begin_step()
    second.step_acc["acabf"][:] = [0.3, 0.3]
    second.step_time = 0.6
    second.commit_step()
    _write_year(second, Q, V, h, bed)             # year 2016, a whole year of it
    second.close()

    assert AnnualOutput.years_on_disk(out) == [2015, 2016]
    with fd.CheckpointFile(AnnualOutput.year_path(out, 2016), "r") as chk:
        acabf = chk.load_function(chk.load_mesh(), name="acabf").dat.data_ro.copy()
    # mean over the whole year, both links pooled: 0.6 m over 1.0 yr
    assert np.allclose(acabf, 0.6)


def test_the_year_in_progress_round_trips_through_a_checkpoint(two_cells, tmp_path):
    r"""state_fields/state_attrs and read_state are the two halves of what
    _save_state writes into the run's own checkpoint."""
    mesh, Q, V = two_cells
    h = [1500.0, 1500.0]
    annual = AnnualOutput(mesh, Q, str(tmp_path / "out" / "annual.h5"),
                          str(tmp_path / "out" / "scalars.csv"),
                          first_year=2015, rho_ratio=RHO_RATIO)
    annual.start_year(_dg(Q, h))
    annual.year_acc["licalvf"][:] = [-2.0, -5.0]
    annual.year_time = 0.7
    state_path = str(tmp_path / "state.h5")
    with fd.CheckpointFile(state_path, "w") as chk:
        chk.save_mesh(mesh)
        for name, f in annual.state_fields().items():
            chk.save_function(f, name=name)
        for key, val in annual.state_attrs().items():
            chk.set_attr("/", key, val)
    annual.close()

    with fd.CheckpointFile(state_path, "r") as chk:
        loaded = chk.load_mesh()
        resume = AnnualOutput.read_state(chk, loaded)
    assert resume["year"] == 2015
    assert resume["year_time"] == pytest.approx(0.7)
    assert resume["series"] == "annual.h5"
    assert resume["skipping"] is False
    assert np.allclose(resume["acc"]["licalvf"], [-2.0, -5.0])
    assert np.allclose(resume["h_year_start"], h)


def test_a_checkpoint_without_output_has_no_state(two_cells, tmp_path):
    r"""A run that had ISMIP7_OUTPUT off writes no such attrs, and the reader
    says so rather than raising."""
    mesh, Q, V = two_cells
    path = str(tmp_path / "plain.h5")
    with fd.CheckpointFile(path, "w") as chk:
        chk.save_mesh(mesh)
        chk.set_attr("/", "t_yr", 2020.0)
    with fd.CheckpointFile(path, "r") as chk:
        assert AnnualOutput.read_state(chk, chk.load_mesh()) is None


def test_a_year_killed_mid_write_leaves_the_banked_years_intact(two_cells, tmp_path):
    r"""Each year is written to <name>.tmp and renamed, so a kill during one
    year cannot damage the years already on disk, and the next link simply
    rewrites the missing one."""
    mesh, Q, V = two_cells
    h = [1500.0, 1500.0]
    bed = [-500.0, -500.0]
    out = str(tmp_path / "out" / "annual.h5")
    scalars = str(tmp_path / "out" / "scalars.csv")

    first = AnnualOutput(mesh, Q, out, scalars, first_year=2015, rho_ratio=RHO_RATIO)
    first.start_year(_dg(Q, h))
    _write_year(first, Q, V, h, bed)          # 2015 banked
    first.close()

    # a job killed part-way through writing 2016 leaves only the temp file
    killed = AnnualOutput.year_path(out, 2016) + ".tmp"
    with open(killed, "wb") as f:
        f.write(b"\x89HDF\r\n\x1a\n truncated")

    assert AnnualOutput.years_on_disk(out) == [2015]
    with fd.CheckpointFile(AnnualOutput.year_path(out, 2015), "r") as chk:
        assert np.allclose(chk.load_function(chk.load_mesh(), name="lithk").dat.data_ro, h)

    # the next link picks up at 2016 and rewrites it
    second = AnnualOutput(mesh, Q, out, scalars, first_year=2016, rho_ratio=RHO_RATIO)
    assert second.year == 2016
    second.start_year(_dg(Q, h))
    _write_year(second, Q, V, h, bed)
    second.close()
    assert AnnualOutput.years_on_disk(out) == [2015, 2016]


def test_a_resume_state_from_another_year_is_refused(two_cells, tmp_path):
    mesh, Q, V = two_cells
    resume = {
        "year": 2016, "year_time": 0.4, "series": "annual.h5",
        "acc": {k: np.zeros(2) for k in AnnualOutput.ACCUMULATORS},
        "h_year_start": np.zeros(2),
    }
    with pytest.raises(ValueError, match="does not belong to it"):
        AnnualOutput(mesh, Q, str(tmp_path / "out" / "annual.h5"),
                     str(tmp_path / "out" / "scalars.csv"),
                     first_year=2020.2, rho_ratio=RHO_RATIO, resume=resume)


def test_an_unclean_kill_past_the_checkpoint_discards_the_stale_years(two_cells, tmp_path):
    r"""State checkpoints are written every 5 model years while years are
    banked every 1, so a node failure leaves years on disk that the restart's
    trajectory never simulated. They belong to an abandoned run and are
    discarded and re-simulated rather than spliced into the series."""
    mesh, Q, V = two_cells
    h = [1500.0, 1500.0]
    bed = [-500.0, -500.0]
    out = str(tmp_path / "out" / "annual.h5")
    scalars = str(tmp_path / "out" / "scalars.csv")

    first = AnnualOutput(mesh, Q, out, scalars, first_year=2015, rho_ratio=RHO_RATIO)
    first.start_year(_dg(Q, h))
    for _ in range(4):
        _write_year(first, Q, V, h, bed)              # 2015-2018 banked
    first.close()
    assert AnnualOutput.years_on_disk(out) == [2015, 2016, 2017, 2018]

    # the newest state checkpoint is at t=2017.0, i.e. inside year 2017, and
    # carries that year's (empty) accumulation state: 2017 and 2018 were
    # simulated after it and are stale
    resume = {
        "year": 2017, "year_time": 0.0, "series": "annual.h5",
        "acc": {k: np.zeros(2) for k in AnnualOutput.ACCUMULATORS},
        "h_year_start": np.asarray(h),
    }
    second = AnnualOutput(mesh, Q, out, scalars, first_year=2017,
                          rho_ratio=RHO_RATIO, resume=resume)
    assert AnnualOutput.years_on_disk(out) == [2015, 2016]
    assert second.year == 2017
    import csv
    with open(scalars) as f:
        assert [int(r["year"]) for r in csv.DictReader(f)] == [2015, 2016]

    _write_year(second, Q, V, h, bed)                 # 2017 re-simulated
    second.close()
    assert AnnualOutput.years_on_disk(out) == [2015, 2016, 2017]
    with open(scalars) as f:
        assert [int(r["year"]) for r in csv.DictReader(f)] == [2015, 2016, 2017]


def test_a_resume_that_would_leave_a_hole_is_refused(two_cells, tmp_path):
    r"""Discarding stale years is recovery; skipping forward past unwritten
    years would submit a series with a gap, so it stays an error."""
    mesh, Q, V = two_cells
    h = [1500.0, 1500.0]
    out = str(tmp_path / "out" / "annual.h5")
    scalars = str(tmp_path / "out" / "scalars.csv")
    first = AnnualOutput(mesh, Q, out, scalars, first_year=2015, rho_ratio=RHO_RATIO)
    first.start_year(_dg(Q, h))
    _write_year(first, Q, V, h, [-500.0, -500.0])     # year 2015 on disk
    first.close()
    resume = {
        "year": 2030, "year_time": 0.0, "series": "annual.h5",
        "acc": {k: np.zeros(2) for k in AnnualOutput.ACCUMULATORS},
        "h_year_start": np.asarray(h),
    }
    with pytest.raises(ValueError, match="hole in it"):
        AnnualOutput(mesh, Q, out, scalars, first_year=2030,
                     rho_ratio=RHO_RATIO, resume=resume)


def test_a_foreign_resume_refuses_to_rewrite_a_banked_series(two_cells, tmp_path):
    r"""A projection restarting from the historical endpoint carries the
    HISTORICAL run's accumulation state, which belongs to another submitted
    series. That is a cold start for this output, not its resume, so the
    banked years stand."""
    mesh, Q, V = two_cells
    h = [1500.0, 1500.0]
    bed = [-500.0, -500.0]
    proj = str(tmp_path / "out" / "ssp126_2500_ismip7_annual.h5")
    proj_scalars = str(tmp_path / "out" / "ssp126_2500_ismip7_scalars.csv")

    banked = AnnualOutput(mesh, Q, proj, proj_scalars, first_year=2015, rho_ratio=RHO_RATIO)
    banked.start_year(_dg(Q, h))
    for _ in range(3):
        _write_year(banked, Q, V, h, bed)             # 2015-2017 banked
    banked.close()

    # the historical run's own state, stamped with ITS series
    hist = AnnualOutput(mesh, Q, str(tmp_path / "out" / "hist_2500_ismip7_annual.h5"),
                        str(tmp_path / "out" / "hist_2500_ismip7_scalars.csv"),
                        first_year=2015, rho_ratio=RHO_RATIO)
    hist.start_year(_dg(Q, h))
    hist_attrs = hist.state_attrs()
    hist_fields = hist.state_fields()
    hist.close()
    assert hist_attrs[AnnualOutput.STATE_SERIES] == "hist_2500_ismip7_annual.h5"

    foreign = {
        "year": hist_attrs[AnnualOutput.STATE_YEAR],
        "year_time": hist_attrs[AnnualOutput.STATE_YEAR_TIME],
        "series": hist_attrs[AnnualOutput.STATE_SERIES],
        "acc": {k: hist_fields[AnnualOutput.STATE_PREFIX + k].dat.data_ro.copy()
                for k in AnnualOutput.ACCUMULATORS},
        "h_year_start": hist_fields[AnnualOutput.STATE_THICKNESS].dat.data_ro.copy(),
    }
    with pytest.raises(ValueError, match="move the series aside"):
        AnnualOutput(mesh, Q, proj, proj_scalars, first_year=2015,
                     rho_ratio=RHO_RATIO, resume=foreign)
    assert AnnualOutput.years_on_disk(proj) == [2015, 2016, 2017]
    import csv
    with open(proj_scalars) as f:
        assert [int(r["year"]) for r in csv.DictReader(f)] == [2015, 2016, 2017]


def test_a_foreign_resume_starts_a_new_series_where_none_is_banked(two_cells, tmp_path):
    r"""The ordinary first link of a projection: it restarts from the
    historical endpoint, and its own series does not exist yet, so it starts
    clean rather than carrying the historical accumulators into year 2015."""
    mesh, Q, V = two_cells
    h = [1500.0, 1500.0]
    proj = str(tmp_path / "out" / "ssp126_2500_ismip7_annual.h5")
    foreign = {
        "year": 2015, "year_time": 0.4,
        "series": "hist_2500_ismip7_annual.h5",
        "acc": {k: np.full(2, 7.0) for k in AnnualOutput.ACCUMULATORS},
        "h_year_start": np.full(2, 999.0),
    }
    started = AnnualOutput(mesh, Q, proj,
                           str(tmp_path / "out" / "ssp126_2500_ismip7_scalars.csv"),
                           first_year=2015, rho_ratio=RHO_RATIO, resume=foreign)
    assert started.year == 2015
    assert started.h_year_start is None              # not mid-year: a fresh series
    assert started.year_time == 0.0
    assert np.allclose(started.year_acc["acabf"], 0.0)
    started.start_year(_dg(Q, h))
    _write_year(started, Q, V, h, [-500.0, -500.0])
    started.close()
    assert AnnualOutput.years_on_disk(proj) == [2015]


def test_a_cold_start_refuses_to_rewrite_a_banked_series(two_cells, tmp_path):
    r"""Without the ISMIP7 accumulation state there is no evidence the years
    on disk are wrong, and they are the only copy of what gets submitted, so
    a run that would rewrite them stops instead of deleting them, and names
    the files so the operator can move them aside deliberately."""
    mesh, Q, V = two_cells
    h = [1500.0, 1500.0]
    bed = [-500.0, -500.0]
    out = str(tmp_path / "out" / "annual.h5")
    scalars = str(tmp_path / "out" / "scalars.csv")

    first = AnnualOutput(mesh, Q, out, scalars, first_year=2015, rho_ratio=RHO_RATIO)
    first.start_year(_dg(Q, h))
    for _ in range(3):
        _write_year(first, Q, V, h, bed)              # 2015-2017 banked
    first.close()

    with pytest.raises(ValueError, match=r"2015-2017.*move the series aside"):
        AnnualOutput(mesh, Q, out, scalars, first_year=2015, rho_ratio=RHO_RATIO)
    assert AnnualOutput.years_on_disk(out) == [2015, 2016, 2017]
    import csv
    with open(scalars) as f:
        assert [int(r["year"]) for r in csv.DictReader(f)] == [2015, 2016, 2017]

    # the operator moves them aside; the run then starts clean
    import os
    for yr in (2015, 2016, 2017):
        os.remove(AnnualOutput.year_path(out, yr))
    os.remove(scalars)
    again = AnnualOutput(mesh, Q, out, scalars, first_year=2015, rho_ratio=RHO_RATIO)
    assert AnnualOutput.years_on_disk(out) == []
    assert again.year == 2015
    again.close()


def test_a_cold_start_beside_earlier_years_is_untouched(two_cells, tmp_path):
    r"""Only years the run would rewrite are in question: a cold start that
    appends after what is on disk is an ordinary continuation."""
    mesh, Q, V = two_cells
    h = [1500.0, 1500.0]
    out = str(tmp_path / "out" / "annual.h5")
    scalars = str(tmp_path / "out" / "scalars.csv")
    first = AnnualOutput(mesh, Q, out, scalars, first_year=2015, rho_ratio=RHO_RATIO)
    first.start_year(_dg(Q, h))
    _write_year(first, Q, V, h, [-500.0, -500.0])     # 2015 banked
    first.close()

    second = AnnualOutput(mesh, Q, out, scalars, first_year=2016, rho_ratio=RHO_RATIO)
    assert AnnualOutput.years_on_disk(out) == [2015]
    assert second.year == 2016
    second.close()


def test_a_midyear_cold_start_does_not_bank_the_partial_year(two_cells, tmp_path):
    r"""Enabling ISMIP7_OUTPUT on a link whose predecessor ran without it lands
    part-way through a year with no record of its earlier months. Banking that
    fraction as the year's mean would submit a wrong annual value, so the
    partial year is dropped and accumulation begins at the next 1 January."""
    mesh, Q, V = two_cells
    h = [1500.0, 1500.0]
    bed = [-500.0, -500.0]
    out = str(tmp_path / "out" / "annual.h5")
    scalars = str(tmp_path / "out" / "scalars.csv")

    annual = AnnualOutput(mesh, Q, out, scalars, first_year=2050.4, rho_ratio=RHO_RATIO)
    assert annual.year == 2051                  # the next whole year, not 2050
    annual.start_year(_dg(Q, h))

    # the 0.6 yr remaining in 2050
    annual.begin_step()
    annual.step_acc["acabf"][:] = [6.0, 6.0]
    annual.step_time = 0.6
    annual.commit_step()
    _write_year(annual, Q, V, h, bed)           # the year end at t=2051.0
    assert AnnualOutput.years_on_disk(out) == []        # 2050 is not banked
    assert annual.year == 2051
    assert np.allclose(annual.year_acc["acabf"], 0.0)   # the partial year is dropped

    # a whole year of 2051 follows and IS banked, as a full-year mean
    annual.begin_step()
    annual.step_acc["acabf"][:] = [2.0, 2.0]
    annual.step_time = 1.0
    annual.commit_step()
    _write_year(annual, Q, V, h, bed)
    annual.close()

    assert AnnualOutput.years_on_disk(out) == [2051]
    with fd.CheckpointFile(AnnualOutput.year_path(out, 2051), "r") as chk:
        acabf = chk.load_function(chk.load_mesh(), name="acabf").dat.data_ro.copy()
    assert np.allclose(acabf, 2.0)


def test_a_whole_year_cold_start_banks_its_first_year(two_cells, tmp_path):
    r"""The ordinary cold start begins on 1 January, so nothing is skipped."""
    mesh, Q, V = two_cells
    h = [1500.0, 1500.0]
    out = str(tmp_path / "out" / "annual.h5")
    annual = AnnualOutput(mesh, Q, out, str(tmp_path / "out" / "scalars.csv"),
                          first_year=2015.0, rho_ratio=RHO_RATIO)
    assert annual.year == 2015
    annual.start_year(_dg(Q, h))
    _write_year(annual, Q, V, h, [-500.0, -500.0])
    annual.close()
    assert AnnualOutput.years_on_disk(out) == [2015]


def test_a_midyear_resume_still_carries_its_months(two_cells, tmp_path):
    r"""The skip is for a COLD start only: a resume carries the months already
    accumulated, so its year must not be pushed forward."""
    mesh, Q, V = two_cells
    resume = {
        "year": 2050, "year_time": 0.4, "series": "annual.h5",
        "acc": {k: np.zeros(2) for k in AnnualOutput.ACCUMULATORS},
        "h_year_start": np.full(2, 1500.0),
    }
    annual = AnnualOutput(mesh, Q, str(tmp_path / "out" / "annual.h5"),
                          str(tmp_path / "out" / "scalars.csv"),
                          first_year=2050.4, rho_ratio=RHO_RATIO, resume=resume)
    assert annual.year == 2050
    assert annual.year_time == pytest.approx(0.4)
    annual.close()


def _round_trip_state(annual, mesh, path):
    r"""Write the run's ISMIP7 state the way _save_state does, read it back the
    way setup_model does."""
    with fd.CheckpointFile(path, "w") as chk:
        chk.save_mesh(mesh)
        for name, f in annual.state_fields().items():
            chk.save_function(f, name=name)
        for key, val in annual.state_attrs().items():
            chk.set_attr("/", key, val)
    with fd.CheckpointFile(path, "r") as chk:
        return AnnualOutput.read_state(chk, chk.load_mesh())


def test_a_checkpoint_inside_a_skipped_year_resumes_the_skip(two_cells, tmp_path):
    r"""A run can stop (wall clock, stall, crash) inside the partial year it
    is not banking. The checkpoint records the year the MODEL is in plus the
    skip flag, so the resume re-enters the skip; stamping the first bankable
    year instead would not match the timeline written beside it and the resume
    would refuse its own state."""
    mesh, Q, V = two_cells
    h = [1500.0, 1500.0]
    bed = [-500.0, -500.0]
    out = str(tmp_path / "out" / "annual.h5")
    scalars = str(tmp_path / "out" / "scalars.csv")

    first = AnnualOutput(mesh, Q, out, scalars, first_year=2050.4, rho_ratio=RHO_RATIO)
    assert first.year == 2051                      # the first year it will bank
    first.start_year(_dg(Q, h))
    first.begin_step()
    first.step_acc["acabf"][:] = [9.0, 9.0]
    first.step_time = 0.3
    first.commit_step()

    # the run stops at t=2050.7, still inside the skipped year
    state = _round_trip_state(first, mesh, str(tmp_path / "state.h5"))
    first.close()

    # the resume must accept its own state: stamping the first bankable year
    # made this raise "does not belong to it" on every retry
    second = AnnualOutput(mesh, Q, out, scalars, first_year=2050.7,
                          rho_ratio=RHO_RATIO, resume=state)
    assert state["year"] == 2050                   # the model's own year
    assert state["skipping"] is True
    assert second.year == 2051                     # still the first full year
    assert np.allclose(second.year_acc["acabf"], 0.0)   # the partial window is dropped

    # the year end at t=2051.0 is still skipped, and 2051 is banked in full
    second.start_year(_dg(Q, h))
    _write_year(second, Q, V, h, bed)
    assert AnnualOutput.years_on_disk(out) == []
    second.begin_step()
    second.step_acc["acabf"][:] = [4.0, 4.0]
    second.step_time = 1.0
    second.commit_step()
    _write_year(second, Q, V, h, bed)
    second.close()

    assert AnnualOutput.years_on_disk(out) == [2051]
    with fd.CheckpointFile(AnnualOutput.year_path(out, 2051), "r") as chk:
        acabf = chk.load_function(chk.load_mesh(), name="acabf").dat.data_ro.copy()
    assert np.allclose(acabf, 4.0)


def test_a_checkpoint_past_the_skipped_year_resumes_normally(two_cells, tmp_path):
    r"""Once the skip is over the state is an ordinary accumulation again: the
    flag is clear and the year is the one being banked."""
    mesh, Q, V = two_cells
    h = [1500.0, 1500.0]
    out = str(tmp_path / "out" / "annual.h5")
    scalars = str(tmp_path / "out" / "scalars.csv")

    annual = AnnualOutput(mesh, Q, out, scalars, first_year=2050.4, rho_ratio=RHO_RATIO)
    annual.start_year(_dg(Q, h))
    _write_year(annual, Q, V, h, [-500.0, -500.0])     # the skipped year end
    annual.begin_step()
    annual.step_acc["acabf"][:] = [1.0, 1.0]
    annual.step_time = 0.5
    annual.commit_step()

    state = _round_trip_state(annual, mesh, str(tmp_path / "state.h5"))
    annual.close()
    assert state["year"] == 2051 and state["skipping"] is False

    resumed = AnnualOutput(mesh, Q, out, scalars, first_year=2051.5,
                           rho_ratio=RHO_RATIO, resume=state)
    assert resumed.year == 2051
    assert resumed.year_time == pytest.approx(0.5)
    assert np.allclose(resumed.year_acc["acabf"], [1.0, 1.0])
    resumed.close()


# --- the grounding-line flux carries a sign ----------------------------------
#
# ``ligroundf`` takes the grounded sheet as its reference (discussion #22,
# settled 22 September 2026): positive for grounded ice going afloat, negative
# where floating ice flows onto grounded ice, booked in the floating cell.

def _strip(n):
    r"""``n`` unit quadrilateral cells in a row, plus the dof order that
    reads them left to right."""
    mesh = fd.RectangleMesh(n, 1, float(n), 1.0, quadrilateral=True)
    Q = fd.FunctionSpace(mesh, "DG", 0)
    V = fd.VectorFunctionSpace(mesh, "CG", 1)
    x = fd.Function(Q).interpolate(fd.SpatialCoordinate(mesh)[0]).dat.data_ro
    return mesh, Q, V, np.argsort(x)


def _one_advance(mesh, Q, V, grounded, ux, tmp_path, bed=None):
    r"""One unit advance of a unit-thick sheet at velocity (ux, 0), no sources.
    Returns the writer and the ``ligroundf`` booked per cell in dof order,
    which with unit cells is the facet flux u h L. ``bed`` (per cell, dof
    order) goes to the writer."""
    annual = AnnualOutput(mesh, Q, str(tmp_path / "out" / "annual.h5"),
                          str(tmp_path / "out" / "scalars.csv"),
                          first_year=2015, rho_ratio=RHO_RATIO, bed=bed)
    h = _dg(Q, np.ones(len(grounded)))
    u = fd.Function(V)
    u.dat.data[:, 0] = ux
    annual.start_year(h)
    annual.begin_step()
    annual.book_advance(1.0, fd.Constant(0.0), fd.Constant(0.0), None, h, u, grounded)
    annual.commit_step()
    return annual, annual.year_acc["ligroundf"].copy(), h, u


def test_ligroundf_is_positive_for_grounded_ice_going_afloat(tmp_path):
    mesh, Q, V, order = _strip(4)
    grounded = np.zeros(4, dtype=bool)
    grounded[order[:2]] = True                     # x < 2 grounded, x > 2 afloat
    _, booked, _, _ = _one_advance(mesh, Q, V, grounded, 1.0, tmp_path)
    assert np.allclose(booked[order], [0.0, 0.0, 1.0, 0.0])


def test_ligroundf_is_negative_where_floating_ice_grounds(tmp_path):
    r"""The same facet with the flow reversed books the same magnitude with
    the opposite sign, in the same floating cell, and the year's scalar
    carries it."""
    mesh, Q, V, order = _strip(4)
    grounded = np.zeros(4, dtype=bool)
    grounded[order[:2]] = True
    annual, booked, h, u = _one_advance(mesh, Q, V, grounded, -1.0, tmp_path)
    assert np.allclose(booked[order], [0.0, 0.0, -1.0, 0.0])

    from icepack2_tools.ismip7_output import SECONDS_PER_YEAR
    bed = _dg(Q, -1000.0 * np.ones(4))
    annual.year_end(h, h, bed, u, fd.Function(V), grounded, np.ones(4, dtype=bool))
    annual.close()
    import csv
    with open(tmp_path / "out" / "scalars.csv") as f:
        row = next(iter(csv.DictReader(f)))
    assert float(row["tendligroundf"]) == pytest.approx(-RHO_I / SECONDS_PER_YEAR * AF2_POLE,
                                                        rel=1e-6)   # the csv carries 7 digits


def test_an_ice_rumple_nets_to_zero(tmp_path):
    r"""Floating, grounded, floating under uniform flow: what flows onto the
    pinning point books negative upstream of it and positive downstream, so
    the discharge sums to zero. Booking the grounded cell's outflow alone
    gave +1 here, the rumple's throughput counted as loss."""
    mesh, Q, V, order = _strip(5)
    grounded = np.zeros(5, dtype=bool)
    grounded[order[2]] = True
    _, booked, _, _ = _one_advance(mesh, Q, V, grounded, 1.0, tmp_path)
    assert np.allclose(booked[order], [0.0, -1.0, 0.0, 1.0, 0.0])
    assert abs(booked.sum()) < 1e-12


def test_ice_pushed_onto_ice_free_land_is_no_grounding_line_flux(tmp_path):
    r"""A land-terminating margin pushing ice onto ice-free land beyond a
    pinned front. The flotation test reads the ice-free land cells as afloat
    (a height above flotation of exactly 0), and the booking counted the flow
    as grounding-line flux (issue #136, 8 to 35 Gt/yr in the 25 km
    rehearsal). Given the bed, a cell on a bed at or above sea level counts as
    grounded; grounded ice going afloat onto a marine bed still books."""
    mesh, Q, V, order = _strip(4)
    grounded = np.zeros(4, dtype=bool)
    grounded[order[:2]] = True                     # as the flotation test reads it
    _, booked, _, _ = _one_advance(mesh, Q, V, grounded, 1.0, tmp_path / "land",
                                   bed=np.full(4, 100.0))
    assert np.allclose(booked, 0.0)
    bed = np.full(4, 100.0)
    bed[order[2:]] = -500.0                        # the margin reaches the sea
    _, booked, _, _ = _one_advance(mesh, Q, V, grounded, 1.0, tmp_path / "coast", bed=bed)
    assert np.allclose(booked[order], [0.0, 0.0, 1.0, 0.0])


def test_the_writer_refuses_a_bed_of_the_wrong_size(tmp_path):
    mesh, Q, _, _ = _strip(4)
    with pytest.raises(ValueError, match="owns 4 cells"):
        AnnualOutput(mesh, Q, str(tmp_path / "annual.h5"), str(tmp_path / "scalars.csv"),
                     first_year=2015, rho_ratio=RHO_RATIO, bed=np.zeros(3))


# The fluxes book what the transport applied. The positivity limiter holds
# back the part of a net sink that would draw a cell below the floor; that
# part comes off the SMB, melt and reference sinks in proportion, so a cell
# held at the floor does not report melt of ice it did not have.

def _book(two_cells, tmp_path, smb, melt, withheld, a_ref=None):
    r"""One unit advance of two floating 100 m cells with no flow, the given
    per-cell sources (m/yr) and withheld sink; returns the booked year."""
    mesh, Q, V = two_cells
    annual = AnnualOutput(mesh, Q, str(tmp_path / "out" / "annual.h5"),
                          str(tmp_path / "out" / "scalars.csv"),
                          first_year=2015, rho_ratio=RHO_RATIO)
    h = _dg(Q, [100.0, 100.0])
    annual.start_year(h)
    annual.begin_step()
    annual.book_advance(1.0, _dg(Q, smb), _dg(Q, melt),
                        None if a_ref is None else _dg(Q, a_ref), h, fd.Function(V),
                        np.zeros(2, dtype=bool),
                        withheld=None if withheld is None else np.asarray(withheld, float))
    annual.commit_step()
    return annual.year_acc


def test_without_a_withheld_sink_the_requested_sources_are_booked(two_cells, tmp_path):
    acc = _book(two_cells, tmp_path, [-2.0, 3.0], [6.0, 5.0], None)
    assert np.allclose(acc["acabf"], [-2.0, 3.0])
    assert np.allclose(acc["libmassbffl"], [-6.0, -5.0])


def test_a_withheld_sink_comes_off_the_melt(two_cells, tmp_path):
    r"""Melt alone: the cell books the melt it had ice for."""
    acc = _book(two_cells, tmp_path, [0.0, 0.0], [10.0, 10.0], [4.0, 0.0])
    assert np.allclose(acc["libmassbffl"], [-6.0, -10.0])
    assert np.allclose(acc["acabf"], [0.0, 0.0])


def test_a_withheld_sink_splits_in_proportion_and_closes(two_cells, tmp_path):
    r"""Cell 0 sinks 2 of SMB, 6 of melt and 2 of reference, and the limiter
    holds back 5, half of them, so each keeps half. Cell 1 gains 3 of SMB,
    which the limiter never touches, against 5 of melt, 2 held back. Either
    way the booked sources sum to the requested source plus what was held
    back."""
    smb, melt, ref, held = [-2.0, 3.0], [6.0, 5.0], [-2.0, 0.0], [5.0, 2.0]
    acc = _book(two_cells, tmp_path, smb, melt, held, a_ref=ref)
    assert np.allclose(acc["acabf"], [-1.0, 3.0])
    assert np.allclose(acc["libmassbffl"], [-3.0, -3.0])
    assert np.allclose(acc["acabf_correction"], [-1.0, 0.0])
    applied = acc["acabf"] + acc["libmassbffl"] + acc["acabf_correction"]
    requested = np.array(smb) - np.array(melt) + np.array(ref)
    assert np.allclose(applied, requested + np.array(held))


# Melt where no ice floats at year end (issue #109, option 3 of the 25
# September 2026 meeting, and issue #136 of 28 September). The request fills
# libmassbffl wherever no ice floats at year end, so the melt of every cell the
# written masks show without floating ice at year end is written as lifmassbf,
# which the request never fills. Before that, on a cell holding no ice at
# either end of the year, the reference and the snowfall it received are
# booked against the melt and SMB that only cancelled them (issue #136).

def _front_fields(tmp_path, cases, year_time=1.0):
    r"""One year of a strip of unit cells, one per case, with the year's sums
    set as the transport books them, then ``year_end`` with the forward's own
    masks. A case gives the start and end thickness ``h0`` and ``h1`` (m),
    the bed (m, marine by default), the surface ``orog`` at year end (m, h1 by
    default), whether it is ``grounded`` at year end (afloat by default), and
    the year's sums in metres of ice: ``melt`` (booked ``libmassbffl``),
    ``smb``, ``ref`` (the apparent-MB reference) and ``calv`` (``licalvf``).
    Returns the year file's flux fields in case order (m/yr), its scalars row
    and its front-melt stamp."""
    import csv
    n = len(cases)
    mesh, Q, V, order = _strip(n)
    out = str(tmp_path / "out" / "annual.h5")
    annual = AnnualOutput(mesh, Q, out, str(tmp_path / "out" / "scalars.csv"),
                          first_year=2015, rho_ratio=RHO_RATIO)

    def column(key, default=0.0):
        v = np.zeros(n)
        v[order] = [c.get(key, default) for c in cases]
        return v

    annual.start_year(_dg(Q, column("h0")))
    annual.begin_step()
    for acc, key in (("libmassbffl", "melt"), ("acabf", "smb"),
                     ("acabf_correction", "ref"), ("licalvf", "calv")):
        annual.step_acc[acc][:] = column(key)
    annual.step_time = year_time
    annual.commit_step()
    h1 = _dg(Q, column("h1"))
    orog = np.zeros(n)
    orog[order] = [c.get("orog", c.get("h1", 0.0)) for c in cases]
    grounded = column("grounded", False) > 0.5
    annual.year_end(h1, _dg(Q, orog), _dg(Q, column("bed", -500.0)), fd.Function(V),
                    fd.Function(V), grounded, h1.dat.data_ro > AnnualOutput.ICE_THICKNESS)
    annual.close()
    names = ("libmassbffl", "lifmassbf", "acabf", "acabf_correction", "licalvf")
    with fd.CheckpointFile(AnnualOutput.year_path(out, 2015), "r") as chk:
        m = chk.load_mesh()
        got = {k: chk.load_function(m, name=k) for k in names}
        stamp = str(chk.get_attr("/", FRONT_MELT_ATTR)) if chk.has_attr("/", FRONT_MELT_ATTR) else None
    x = fd.Function(got["lifmassbf"].function_space()).interpolate(fd.SpatialCoordinate(m)[0])
    left_to_right = np.argsort(x.dat.data_ro)
    with open(tmp_path / "out" / "scalars.csv") as f:
        row = next(iter(csv.DictReader(f)))
    return {k: got[k].dat.data_ro[left_to_right].copy() for k in names}, row, stamp


def _front_year(tmp_path, cases, year_time=1.0):
    r""":func:`_front_fields`, returning the year file's ``libmassbffl`` and
    ``lifmassbf``, its scalars row and its stamp."""
    fields, row, stamp = _front_fields(tmp_path, cases, year_time)
    return fields["libmassbffl"], fields["lifmassbf"], row, stamp


# The cases, year sums in metres of ice. The inflow a cell kept is the
# remainder of its thickness budget, h1 - h0 - smb - ref - melt - calv.
INFLOW = dict(h0=0.0, h1=0.5, melt=-40.0, smb=-5.0)          # kept 45.5, all of the supply
MIXED = dict(h0=0.0, h1=0.0, melt=-40.0, ref=30.0)            # kept 10 of a supply of 40
SNOWFALL = dict(h0=0.0, h1=0.0, smb=3.0, melt=-3.0)           # snow on open water, melted
SNOW_AND_INFLOW = dict(h0=0.0, h1=0.0, smb=3.0, melt=-13.0)   # kept 10 of a supply of 13
CALVED_INFLOW = dict(h0=0.0, h1=0.0, smb=2.0, melt=-2.0, calv=-7.0)
# every cell with no floating ice at year end books its melt as lifmassbf
NO_FLOATING_ICE_AT_YEAR_END = {
    "shelf gone": dict(h0=100.0, h1=0.0, melt=-100.0),
    "grounded at year end": dict(h0=100.0, h1=100.0, melt=-5.0, grounded=True),
    "within 1 cm of the bed": dict(h0=100.0, h1=100.0, melt=-5.0, bed=-90.0, orog=10.005),
    "own residue": dict(h0=0.8, h1=0.0, melt=-0.8),
    "at the ice threshold": dict(h0=1.0, h1=1.0, melt=-30.0),
}
# and a cell floating at year end, or refreezing anywhere, keeps it basal
STAYS_BASAL = {
    "ice at year end": dict(h0=0.0, h1=5.0, melt=-20.0),
    "shelf": dict(h0=100.0, h1=95.0, melt=-10.0),
    "refreezing": dict(h0=0.0, h1=0.3, melt=0.3),
}


def test_the_melt_of_inflow_into_an_empty_cell_is_front_melt(tmp_path):
    r"""Grounded ice that goes afloat into an emptied shelf cell and melts
    there: all of the melt is front melt, and the year's scalar carries it
    with its sign."""
    lib, lif, row, stamp = _front_year(tmp_path, [INFLOW])
    assert lif == pytest.approx([-40.0]) and lib == pytest.approx([0.0], abs=1e-12)
    assert float(row["tendlifmassbf"]) == pytest.approx(
        -40.0 * AF2_POLE * RHO_I / SECONDS_PER_YEAR, rel=1e-12)
    assert float(row["tendlibmassbffl"]) == pytest.approx(0.0, abs=1e-12)
    assert stamp == FRONT_MELT


def test_the_reference_an_empty_cell_received_is_booked_against_its_melt(tmp_path):
    r"""An emptied shelf cell under a pinned front keeps receiving the frozen
    reference, and its melt removes it in the step it arrives (issue #136).
    The 30 m of melt that only cancelled the reference come off the melt and
    the reference; what stays is the melt of the 10 m of inflow, all of it
    front melt, and no reference is booked."""
    f, row, _ = _front_fields(tmp_path, [MIXED])
    assert f["lifmassbf"] == pytest.approx([-10.0])
    assert f["libmassbffl"] == pytest.approx([0.0], abs=1e-12)
    assert f["acabf_correction"] == pytest.approx([0.0], abs=1e-12)
    assert f["acabf"] == pytest.approx([0.0], abs=1e-12)
    assert float(row["tendlibmassbffl"]) == pytest.approx(0.0, abs=1e-12)


def test_a_reference_fed_land_cell_books_no_smb(tmp_path):
    r"""Emptied land cells: one the reference alone feeds, whose negative SMB
    removes it (the 25 km rehearsal's ssp585 had 9 to 13 such cells late in
    the run), and one where the reference outruns the SMB and the rest flows
    on. Neither books SMB of the reference, and land takes no front melt."""
    fed = dict(h0=0.0, h1=0.0, smb=-4.0, ref=4.0, bed=50.0)
    outrun = dict(h0=0.0, h1=0.0, smb=-1.0, ref=4.0, bed=50.0)     # 3 m flows on
    f, _, _ = _front_fields(tmp_path, [fed, outrun])
    assert f["acabf"] == pytest.approx([0.0, 0.0], abs=1e-12)
    assert f["acabf_correction"] == pytest.approx([0.0, 3.0])
    assert f["libmassbffl"] == pytest.approx([0.0, 0.0], abs=1e-12)
    assert f["lifmassbf"] == pytest.approx([0.0, 0.0], abs=1e-12)


def test_snowfall_on_open_water_is_booked_against_its_melt(tmp_path):
    r"""Snow that falls on a cell holding no ice at either end and melts
    there was never ice of the ice sheet (issue #136): neither acabf nor the
    melt fields carry it. Where inflow fed the same melt, its share stays as
    front melt, and ice the front removed stays calving."""
    f, row, _ = _front_fields(tmp_path, [SNOWFALL, SNOW_AND_INFLOW, CALVED_INFLOW])
    assert f["acabf"] == pytest.approx([0.0, 0.0, 0.0], abs=1e-12)
    assert f["libmassbffl"] == pytest.approx([0.0, 0.0, 0.0], abs=1e-12)
    assert f["lifmassbf"] == pytest.approx([0.0, -10.0, 0.0], abs=1e-12)
    assert f["licalvf"] == pytest.approx([0.0, 0.0, -7.0], abs=1e-12)
    assert float(row["tendacabf"]) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("case", sorted(NO_FLOATING_ICE_AT_YEAR_END))
def test_melt_where_no_ice_floats_at_year_end_is_front_melt(tmp_path, case):
    r"""Shelf ice gone within the year, a cell that floated and ends the year
    grounded, one the writer writes as grounded for sitting within 1 cm of
    the bed, a cell's own residue and a cell at the ice threshold: the fill
    would leave each one's melt out of libmassbffl, so all of it is
    lifmassbf (issue #136)."""
    c = NO_FLOATING_ICE_AT_YEAR_END[case]
    lib, lif, _, _ = _front_year(tmp_path, [c])
    assert lif == pytest.approx([c["melt"]]) and lib == pytest.approx([0.0], abs=1e-12)


@pytest.mark.parametrize("case", sorted(STAYS_BASAL))
def test_melt_of_ice_afloat_at_year_end_stays_basal(tmp_path, case):
    r"""A cell afloat at year end keeps its melt in libmassbffl, where the fill
    keeps it, and refreezing never moves: lifmassbf's range excludes it."""
    c = STAYS_BASAL[case]
    lib, lif, _, _ = _front_year(tmp_path, [c])
    assert lif[0] == 0.0
    assert lib[0] == pytest.approx(c["melt"])


def test_the_year_s_sums_become_means_over_its_length(tmp_path):
    r"""A cell whose residue grew over the year, and a year of half length
    whose sums become means over it; the other cells ride along to check
    that each case keeps its own cell. In the first, 7 m of the 21 m of melt
    and negative SMB only cancelled the reference, taken off in proportion
    (20:1), and the melt left is front melt."""
    cases = [dict(h0=0.2, h1=0.9, melt=-20.0, smb=-1.0, ref=7.0), INFLOW, MIXED]
    lib, lif, _, _ = _front_year(tmp_path, cases, year_time=0.5)
    melt = -20.0 + 7.0 * 20.0 / 21.0
    assert lif == pytest.approx([melt / 0.5, -80.0, -20.0])
    assert lib == pytest.approx([0.0, 0.0, 0.0], abs=1e-12)


def test_the_booked_sources_sum_to_what_the_transport_applied(tmp_path):
    r"""The netting moves melt and SMB into the reference and back, and
    snowfall out of both, never the total; the melt fields sum to the melt
    applied less the part that only cancelled the reference (30 of MIXED's
    40 m) or the snowfall (3 m of each snowfall case, 2 of the calved one)."""
    cases = ([INFLOW, MIXED, SNOWFALL, SNOW_AND_INFLOW, CALVED_INFLOW]
             + [NO_FLOATING_ICE_AT_YEAR_END[k] for k in sorted(NO_FLOATING_ICE_AT_YEAR_END)]
             + [STAYS_BASAL[k] for k in sorted(STAYS_BASAL)])
    f, row, _ = _front_fields(tmp_path, cases)
    applied = [c.get("melt", 0.0) + c.get("smb", 0.0) + c.get("ref", 0.0) for c in cases]
    booked = f["libmassbffl"] + f["lifmassbf"] + f["acabf"] + f["acabf_correction"]
    assert booked == pytest.approx(applied)
    melt = [c["melt"] for c in cases]
    melt[1:5] = [-10.0, 0.0, -10.0, 0.0]
    assert f["libmassbffl"] + f["lifmassbf"] == pytest.approx(melt)
    assert np.all(f["lifmassbf"] <= 0.0)
    total = float(row["tendlibmassbffl"]) + float(row["tendlifmassbf"])
    assert total == pytest.approx(sum(melt) * AF2_POLE * RHO_I / SECONDS_PER_YEAR, rel=1e-12)


def test_a_midyear_resume_splits_the_year_as_one_link_would(tmp_path):
    r"""The netting reads the thickness the year began with, which the run's
    checkpoint carries, so a chained link that picks the year up at 0.4
    books the same front melt as one that runs it whole."""
    mesh, Q, V, order = _strip(2)
    first_h = np.zeros(2)
    sums = {"libmassbffl": [-40.0, -40.0], "acabf": [-5.0, 0.0],
            "acabf_correction": [0.0, 30.0]}
    h1 = np.zeros(2); h1[order[0]] = 0.5

    def book(annual, share, time):
        annual.begin_step()
        for k, v in sums.items():
            annual.step_acc[k][order] = np.asarray(v) * share
        annual.step_time = time
        annual.commit_step()

    whole = AnnualOutput(mesh, Q, str(tmp_path / "whole" / "annual.h5"),
                         str(tmp_path / "whole" / "scalars.csv"),
                         first_year=2015, rho_ratio=RHO_RATIO)
    whole.start_year(_dg(Q, first_h))
    book(whole, 1.0, 1.0)
    bed = _dg(Q, [-500.0, -500.0])
    whole.year_end(_dg(Q, h1), _dg(Q, h1), bed, fd.Function(V), fd.Function(V),
                   np.zeros(2, dtype=bool), h1 > AnnualOutput.ICE_THICKNESS)
    whole.close()

    out = str(tmp_path / "linked" / "annual.h5")
    scalars = str(tmp_path / "linked" / "scalars.csv")
    first = AnnualOutput(mesh, Q, out, scalars, first_year=2015, rho_ratio=RHO_RATIO)
    first.start_year(_dg(Q, first_h))
    book(first, 0.4, 0.4)
    resume = _round_trip_state(first, mesh, str(tmp_path / "state.h5"))
    first.close()
    second = AnnualOutput(mesh, Q, out, scalars, first_year=2015.4,
                          rho_ratio=RHO_RATIO, resume=resume)
    book(second, 0.6, 0.6)
    second.year_end(_dg(Q, h1), _dg(Q, h1), bed, fd.Function(V), fd.Function(V),
                    np.zeros(2, dtype=bool), h1 > AnnualOutput.ICE_THICKNESS)
    second.close()

    got = {}
    for tag, path in (("whole", str(tmp_path / "whole" / "annual.h5")), ("linked", out)):
        with fd.CheckpointFile(AnnualOutput.year_path(path, 2015), "r") as chk:
            m = chk.load_mesh()
            got[tag] = {k: chk.load_function(m, name=k).dat.data_ro.copy()
                        for k in ("libmassbffl", "lifmassbf")}
    for k in ("libmassbffl", "lifmassbf"):
        assert got["linked"][k] == pytest.approx(got["whole"][k])
    assert got["whole"]["lifmassbf"][order] == pytest.approx([-40.0, -40.0 * 10.0 / 40.0])


def test_split_melt_never_books_a_gain_and_keeps_the_total():
    r"""On random cells: every value finite, lifmassbf never above zero, zero
    on the floating cells and wherever the year refroze, and libmassbffl
    left with no melt off the floating cells; the two fields sum to the melt
    booked, bit for bit."""
    rng = np.random.default_rng(109)
    n = 20000
    melt = rng.normal(0.0, 30.0, n)
    melt[rng.random(n) < 0.1] = 0.0
    floating = rng.random(n) < 0.5
    lib, lif = split_melt(melt, floating)
    assert np.isfinite(lib).all() and np.isfinite(lif).all()
    assert (lif <= 0.0).all()
    assert (lif[floating] == 0.0).all() and (lif[melt >= 0.0] == 0.0).all()
    assert (lib[~floating] >= 0.0).all()
    assert (lib + lif == melt).all()


def _random_year(seed, n=20000):
    r"""Random year sums in metres of ice, with refreezing, bare cells, both
    signs of reference and SMB, and a cell set."""
    rng = np.random.default_rng(seed)
    melt = rng.normal(0.0, 30.0, n)
    melt[rng.random(n) < 0.1] = 0.0
    smb, ref = rng.normal(0.0, 5.0, n), rng.normal(0.0, 20.0, n)
    ref[rng.random(n) < 0.3] = 0.0
    cells = rng.random(n) < 0.7
    return melt, smb, ref, cells


def test_net_reference_never_books_a_gain_and_keeps_the_sum():
    r"""On random cells: the three keep their sum; the melt and the negative
    SMB only shrink and the reference only falls, to no less than zero; no
    positive reference is left where a sink could have taken it; a negative
    reference, refreezing, a positive SMB and every cell outside the set are
    untouched, bit for bit."""
    melt, smb, ref, cells = _random_year(136)
    m2, s2, r2 = net_reference(melt, smb, ref, cells)
    assert np.allclose(m2 + s2 + r2, melt + smb + ref, rtol=0.0, atol=1e-12)
    assert (m2 >= melt).all() and (s2 >= smb).all() and (r2 <= ref).all()
    assert (m2[melt <= 0.0] <= 1e-12).all() and (s2[smb <= 0.0] <= 1e-12).all()
    assert (r2[(ref > 0.0) & cells] >= -1e-12).all()
    sinks = np.maximum(-m2, 0.0) + np.maximum(-s2, 0.0)
    assert np.minimum(np.maximum(r2, 0.0), sinks)[cells].max() <= 1e-12
    untouched = ~cells | (ref <= 0.0)
    for new, old in ((m2, melt), (s2, smb), (r2, ref)):
        assert (new[untouched] == old[untouched]).all()
    assert (m2[melt >= 0.0] == melt[melt >= 0.0]).all()
    assert (s2[smb >= 0.0] == smb[smb >= 0.0]).all()


def test_net_reference_without_a_reference_returns_its_inputs():
    melt, smb, _, cells = _random_year(1360)
    for ref in (np.zeros_like(melt), -np.abs(melt)):
        m2, s2, r2 = net_reference(melt, smb, ref, cells)
        assert (m2 == melt).all() and (s2 == smb).all() and (r2 == ref).all()


def test_net_snowfall_never_books_a_gain_and_keeps_the_sum():
    r"""On random cells: the two keep their sum; the melt only shrinks and
    the snowfall only falls, to no less than zero; no snowfall is left where
    the melt could have taken it; a negative SMB, refreezing and every cell
    outside the set are untouched, bit for bit."""
    melt, smb, _, cells = _random_year(1362)
    m2, s2 = net_snowfall(melt, smb, cells)
    assert np.allclose(m2 + s2, melt + smb, rtol=0.0, atol=1e-12)
    assert (m2 >= melt).all() and (s2 <= smb).all()
    assert (m2[melt <= 0.0] <= 1e-12).all() and (s2[smb > 0.0] >= -1e-12).all()
    assert np.minimum(np.maximum(s2, 0.0), np.maximum(-m2, 0.0))[cells].max() <= 1e-12
    untouched = ~cells | (smb <= 0.0) | (melt >= 0.0)
    for new, old in ((m2, melt), (s2, smb)):
        assert (new[untouched] == old[untouched]).all()


def test_what_flowed_in_is_the_front_melt_of_a_cell_that_ends_as_it_began():
    r"""Cells holding no ice at either end whose thickness did not change:
    the melt and the negative SMB removed everything supplied, so once the
    reference and the snowfall are booked against them, the melt and the
    negative SMB left are the inflow, and the melt is all lifmassbf."""
    rng = np.random.default_rng(1361)
    n = 5000
    inflow = rng.uniform(0.0, 50.0, n)
    ref = np.where(rng.random(n) < 0.8, rng.uniform(0.0, 40.0, n), 0.0)
    smb = rng.normal(0.0, 5.0, n)
    smb = np.maximum(smb, -(inflow + ref) * 0.9)        # the melt stays non-negative
    melt = -(inflow + ref + smb)                         # everything supplied leaves
    cells = np.ones(n, dtype=bool)
    m2, s2, r2 = net_reference(melt, smb, ref, cells)
    m3, s3 = net_snowfall(m2, s2, cells)
    lib, lif = split_melt(m3, np.zeros(n, dtype=bool))
    assert np.allclose(r2, 0.0, atol=1e-12) and (s3 <= 1e-12).all()
    assert np.allclose(lif + s3, -inflow, rtol=1e-12, atol=1e-9)
    assert np.allclose(lib, 0.0, atol=1e-12)


def test_ice_leaving_across_the_mesh_boundary_is_calving(tmp_path):
    r"""Floating ice flowing out across one end of a strip (issue #136). The
    transport's own outflow term books into licalvf in the cell the ice
    leaves and nowhere else, its integral is the outflux the forward tallies,
    and a year's booked fluxes, the reference among them, close the change in
    lim on the scalars. Before, no field carried the outflux: in the 25 km
    rehearsal's 2015 it was 495 Gt/yr with the reference and 627 without."""
    import csv
    mesh, Q, V, order = _strip(4)
    out = tmp_path / "out"
    annual = AnnualOutput(mesh, Q, str(out / "annual.h5"), str(out / "scalars.csv"),
                          first_year=2015, rho_ratio=RHO_RATIO)
    u = fd.Function(V)
    u.dat.data[:, 0] = 0.1                                  # m/yr, out across x = 4
    smb, melt, ref = fd.Constant(0.5), fd.Constant(2.0), fd.Constant(0.3)
    h, h_old = _dg(Q, np.full(4, 100.0)), fd.Function(Q)
    bed = _dg(Q, np.full(4, -500.0))
    afloat = np.zeros(4, dtype=bool)
    un = fd.dot(u, fd.FacetNormal(mesh))
    un_plus = (un + abs(un)) / 2
    phi, trial = fd.TestFunction(Q), fd.TrialFunction(Q)
    dt = 1.0
    # the forward's transport: implicit upwind DG0 with the applied source
    F = ((trial - h_old) / dt * phi * fd.dx
         + (un_plus("+") * trial("+") - un_plus("-") * trial("-")) * fd.jump(phi) * fd.dS
         + un_plus * trial * phi * fd.ds - (smb - melt + ref) * phi * fd.dx)
    outflux = []
    annual.start_year(h)
    for _ in (2015, 2016):
        h_old.assign(h)
        fd.solve(fd.lhs(F) == fd.rhs(F), h,
                 solver_parameters={"ksp_type": "preonly", "pc_type": "lu"})
        outflux.append(float(fd.assemble(un_plus * h * fd.ds)) * dt)   # m3 of ice
        annual.begin_step()
        annual.book_advance(dt, smb, melt, ref, h, u, afloat)
        annual.commit_step()
        annual.year_end(h, _dg(Q, h.dat.data_ro * (1.0 - RHO_RATIO)), bed, u,
                        fd.Function(V), afloat, h.dat.data_ro > AnnualOutput.ICE_THICKNESS)
    annual.close()
    got = {}
    for yr in (2015, 2016):
        with fd.CheckpointFile(AnnualOutput.year_path(str(out / "annual.h5"), yr), "r") as chk:
            m = chk.load_mesh()
            fields = {k: chk.load_function(m, name=k) for k in ("licalvf", "acabf_correction")}
        x = fd.Function(fields["licalvf"].function_space()).interpolate(
            fd.SpatialCoordinate(m)[0]).dat.data_ro
        got[yr] = {k: f.dat.data_ro[np.argsort(x)].copy() for k, f in fields.items()}
    for k, yr in enumerate((2015, 2016)):
        c = got[yr]["licalvf"]                             # unit cells, left to right
        assert (c[:3] == 0.0).all() and c[3] < 0.0
        assert -c.sum() == pytest.approx(outflux[k], rel=1e-12)
    with open(out / "scalars.csv") as f:
        rows = {int(r["year"]): r for r in csv.DictReader(f)}
    booked = sum(float(rows[2016][k]) for k in ("tendacabf", "tendlibmassbffl",
                                                "tendlifmassbf", "tendlicalvf"))
    booked = booked * SECONDS_PER_YEAR + got[2016]["acabf_correction"].sum() * AF2_POLE * RHO_I
    change = float(rows[2016]["lim"]) - float(rows[2015]["lim"])   # kg over the year
    assert change == pytest.approx(booked, rel=1e-9)
    # nothing but the outflux leaves this strip, and tendlicalvf carries it
    assert float(rows[2016]["tendlicalvf"]) * SECONDS_PER_YEAR == pytest.approx(
        -outflux[1] * AF2_POLE * RHO_I, rel=1e-9)
