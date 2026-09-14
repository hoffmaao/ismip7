r"""What the forward records every year: the annual series survives a chained
resume, and ``limnsw`` is the mass above flotation.

``nots_projection.sbatch`` self-chains, so ``run_simulation`` is entered once
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
