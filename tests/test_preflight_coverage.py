r"""What the preflight gate calls a missing input.

The reader bridges exactly one year past the end of a series (a CESM2-WACCM
atmosphere fetched while the empty 2300 files were withdrawn stops at 2299,
discussion #8), so a tree whose LAST year is absent runs correctly. The gate
has to agree: calling that a missing input reports cores 5 and 7 as BLOCKED
for such a tree. Any other hole stays an error, because the reader raises on
it.

``atm_years`` reads only filenames, so the synthetic trees here are empty
files in the layout ``atmosphere_path`` resolves. The shared mesh/MAP checks
are stubbed out so the status reflects the forcing coverage alone; everything
else, including the printed status line, is the shipped code.

The SMB-elevation feedback is on by default, so every tree here carries the
``dacabfdz`` gradient its core reads unless a test says otherwise.
"""
import importlib
import os
import sys

import numpy as np
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "antarctica", "scripts"))

pytest.importorskip("firedrake")

ESM, SCENARIO = "CESM2-WACCM", "ssp585"
CORE_7 = "core  7"
# the knobs that move a core's period, refuse its start or locate its dH/dt;
# cleared so a developer's shell cannot change what the tests see
KNOBS = ("ISMIP7_T_START", "ISMIP7_T_END", "ISMIP7_GEOMETRY_BACKDATE",
         "ISMIP7_GEOMETRY_SPACE", "ISMIP7_OBS_KIT", "ISMIP7_DHDT_VAR")


def _knobs(monkeypatch, tmp_path, env=None):
    for knob in KNOBS:
        monkeypatch.delenv(knob, raising=False)
    # the dH/dt cache is looked for under the obs root, which by default is
    # the checkout's, where a real cache may sit
    monkeypatch.setenv("ISMIP7_OBS_DATA_ROOT", str(tmp_path / "obs"))
    for knob, value in (env or {}).items():
        monkeypatch.setenv(knob, value)


def _tree(root, years, esm=ESM, scenario=SCENARIO,
          variables=("acabf-anomaly", "acabf", "dacabfdz")):
    r"""Empty files of `variables` for `years`, in the layout the reader uses."""
    for var in variables:
        d = os.path.join(root, esm, scenario, "SDBN1-8000m", var, "v2")
        os.makedirs(d, exist_ok=True)
        for y in years:
            head = f"{var}_AIS_{esm}_{scenario}_SDBN1-8000m_v2_{y}.nc"
            open(os.path.join(d, head), "wb").close()


def _run(monkeypatch, tmp_path, years, ocean_cover=lambda e, s: (2015, 2300),
         scenario=SCENARIO, env=None, dhdt=False, gradient=True, ctrl_gradient=True):
    r"""preflight.main() over a tree covering `years`, returning its output.
    ``dhdt=True`` keeps the real dH/dt lookup, which has tests of its own.
    `gradient` and `ctrl_gradient` add the ``dacabfdz`` that the core and the
    CESM2-WACCM control read."""
    root = str(tmp_path / "ISMIP7" / "AIS")
    os.makedirs(root, exist_ok=True)
    _tree(root, years, scenario=scenario, variables=("acabf-anomaly", "acabf")
          + (("dacabfdz",) if gradient else ()))
    if ctrl_gradient:
        _tree(root, range(2015, 2301), scenario="ctrl", variables=("dacabfdz",))
    monkeypatch.setenv("ISMIP7_DATA_ROOT", root)
    _knobs(monkeypatch, tmp_path, env)
    sys.modules.pop("preflight", None)
    preflight = importlib.import_module("preflight")
    # isolate the forcing-coverage verdict from the mesh/MAP/ocean checks
    monkeypatch.setattr(preflight, "shared_missing", lambda warn: [])
    monkeypatch.setattr(preflight, "racmo_ok", lambda: True)
    monkeypatch.setattr(preflight, "oi_ok", lambda: True)
    monkeypatch.setattr(preflight, "ocean_cover", ocean_cover)
    monkeypatch.setattr(preflight, "pool_status",
                        lambda *a, **k: ("ok", ""))
    if not dhdt:
        monkeypatch.setattr(preflight, "dhdt_missing", lambda: [])
    preflight.main()
    return None


def _core_line(capsys, core=CORE_7):
    line = [ln for ln in capsys.readouterr().out.splitlines() if core in ln]
    assert len(line) == 1, line
    return line[0]


def test_an_absent_final_year_is_a_note_not_a_blocker(monkeypatch, tmp_path, capsys):
    r"""CESM2-WACCM ssp585 covers 2015-2300 and the tree stops at 2299: the
    reader persists 2299 for 2300, so the core is READY."""
    _run(monkeypatch, tmp_path, range(2015, 2300))       # 2015..2299
    line = _core_line(capsys)
    assert "READY" in line, line
    assert "BLOCKED" not in line
    assert "2300 is absent" in line and "persists 2299" in line


def test_a_hole_inside_the_series_still_blocks(monkeypatch, tmp_path, capsys):
    r"""The reader raises on an interior gap, so the gate must not clear it."""
    years = [y for y in range(2015, 2301) if y != 2100]
    _run(monkeypatch, tmp_path, years)
    line = _core_line(capsys)
    assert "BLOCKED" in line, line
    assert "2100..2100" in line


def test_a_short_tree_still_blocks(monkeypatch, tmp_path, capsys):
    r"""Two years short is not the one-year bridge; it is a short download."""
    _run(monkeypatch, tmp_path, range(2015, 2299))       # 2015..2298
    line = _core_line(capsys)
    assert "BLOCKED" in line, line


def test_a_complete_tree_is_ready_without_a_note(monkeypatch, tmp_path, capsys):
    _run(monkeypatch, tmp_path, range(2015, 2301))       # 2015..2300
    line = _core_line(capsys)
    assert "READY" in line, line
    assert "absent" not in line


# --- the years a core needs are its driver's own period ---------------------

def _preflight():
    sys.modules.pop("preflight", None)
    return importlib.import_module("preflight")


def test_each_core_is_checked_over_its_drivers_period():
    r"""Read from the drivers, with the knobs applied as each driver applies
    them: the control starts at 2015 whatever ISMIP7_T_START says."""
    preflight = _preflight()
    period = {core: preflight.driver_period(driver, core, {})
              for core, _, driver, _, _ in preflight.CORES}
    assert period[1] == period[2] == (2003.0, 2015.0)
    assert period[3] == period[4] == (2015.0, 2101.0)
    assert {period[c] for c in (5, 6, 7, 8, 9, 10)} == {(2015.0, 2301.0)}
    assert period[11] == (2003.0, 2026.0)
    moved = {"ISMIP7_T_START": "1990", "ISMIP7_T_END": "2101"}
    assert preflight.driver_period("historical/cesm_waccm.py", 1, moved) == (1990.0, 2101.0)
    assert preflight.driver_period("projections/ocx.py", 11, moved) == (1990.0, 2101.0)
    assert preflight.driver_period("control/run.py", 9, moved) == (2015.0, 2101.0)


def test_a_driver_listed_wrongly_or_written_unreadably_is_refused(tmp_path):
    r"""The guard is worth nothing if it cannot fail."""
    preflight = _preflight()
    with pytest.raises(ValueError, match="runs core 7"):
        preflight.driver_period("projections/ssp585_cesm_waccm.py", 5, {})
    driver = tmp_path / "driver.py"
    driver.write_text("T_START = start_year()\nT_END = 2026.0\n")
    with pytest.raises(ValueError, match="start_year"):
        preflight.driver_period(str(driver), None, {})


def test_a_historical_tree_from_2003_is_ready(monkeypatch, tmp_path, capsys):
    r"""The historicals start in 2003 (issue 117), so a tree holding
    2003-2014 is all they read."""
    _run(monkeypatch, tmp_path, range(2003, 2015), scenario="historical",
         ocean_cover=lambda e, s: (2003, 2014))
    line = _core_line(capsys, "core  1")
    assert "READY" in line and "2003-2014" in line, line


def test_ISMIP7_T_START_moves_the_years_a_historical_reads(monkeypatch, tmp_path, capsys):
    _run(monkeypatch, tmp_path, range(2003, 2015), scenario="historical",
         ocean_cover=lambda e, s: (2003, 2014),
         env={"ISMIP7_T_START": "1990", "ISMIP7_GEOMETRY_BACKDATE": "25"})
    line = _core_line(capsys, "core  1")
    assert "BLOCKED" in line and "13 of 1990-2014 missing (1990..2002)" in line, line


@pytest.mark.parametrize("backdate,status", [(None, "BLOCKED"), ("0", "READY")])
def test_an_1850_start_runs_only_with_the_backdating_named(
        monkeypatch, tmp_path, capsys, backdate, status):
    r"""The driver refuses a start before the Smith dH/dt window unless
    ISMIP7_GEOMETRY_BACKDATE says how many years to undo."""
    env = {"ISMIP7_T_START": "1850"}
    if backdate is not None:
        env["ISMIP7_GEOMETRY_BACKDATE"] = backdate
    _run(monkeypatch, tmp_path, range(1850, 2015), scenario="historical",
         ocean_cover=lambda e, s: (1850, 2014), env=env)
    line = _core_line(capsys, "core  1")
    assert status in line, line
    assert ("past the window" in line) == (backdate is None), line


def test_a_backdated_start_needs_the_dg0_geometry(monkeypatch, tmp_path, capsys):
    _run(monkeypatch, tmp_path, range(2003, 2015), scenario="historical",
         ocean_cover=lambda e, s: (2003, 2014),
         env={"ISMIP7_GEOMETRY_SPACE": "cg1"})
    line = _core_line(capsys, "core  1")
    assert "BLOCKED" in line and "ISMIP7_GEOMETRY_SPACE=dg0" in line, line


# --- cores 9 and 10 melt under the ESM's own ctrl ocean ---------------------

def test_the_control_is_blocked_without_its_esms_ctrl_ocean(monkeypatch, tmp_path, capsys):
    r"""The control reads ``<ESM>/ctrl/ocean`` (icepack/ismip7#107), so that
    tree is what gates cores 9 and 10."""
    asked = []

    def cover(esm, scenario):
        asked.append((esm, scenario))
        return None if scenario == "ctrl" else (2015, 2300)

    _run(monkeypatch, tmp_path, range(2015, 2301), ocean_cover=cover)
    line = _core_line(capsys, "core  9")
    assert "BLOCKED" in line and "CESM2-WACCM/ctrl ocean tf/so" in line, line
    assert ("MRI-ESM2-0", "ctrl") in asked


def test_the_control_is_ready_on_its_ctrl_ocean(monkeypatch, tmp_path, capsys):
    _run(monkeypatch, tmp_path, range(2015, 2301))
    line = _core_line(capsys, "core  9")
    assert "READY" in line, line


# --- the SMB-elevation feedback reads each core's own gradient -------------

def test_core_7_is_blocked_without_its_gradient(monkeypatch, tmp_path, capsys):
    r"""The reader turns an absent variable into zeros, so a run on this tree
    would carry no feedback and say nothing; the gate refuses it as the run
    does (forcing.SMBElevationFeedback.check)."""
    monkeypatch.delenv("ISMIP7_SMB_ELEVATION_FEEDBACK", raising=False)
    _run(monkeypatch, tmp_path, range(2015, 2301), gradient=False)
    out = capsys.readouterr().out
    line = [ln for ln in out.splitlines() if CORE_7 in ln][0]
    assert "BLOCKED" in line and "dacabfdz" in line and "no files" in line, line
    assert "ISMIP7_SMB_ELEVATION_FEEDBACK=0" in line
    assert "SMB-elevation feedback: on" in out


def test_core_7_needs_no_gradient_with_the_feedback_off(monkeypatch, tmp_path, capsys):
    monkeypatch.setenv("ISMIP7_SMB_ELEVATION_FEEDBACK", "0")
    _run(monkeypatch, tmp_path, range(2015, 2301), gradient=False)
    out = capsys.readouterr().out
    line = [ln for ln in out.splitlines() if CORE_7 in ln][0]
    assert "READY" in line and "dacabfdz" not in line, line
    assert "SMB-elevation feedback: off" in out


def test_the_control_reads_the_ctrl_gradient(monkeypatch, tmp_path, capsys):
    r"""The control's feedback reads the ESM's ``ctrl`` gradient, so that tree
    gates cores 9 and 10 even with the scenario gradients all present."""
    monkeypatch.delenv("ISMIP7_SMB_ELEVATION_FEEDBACK", raising=False)
    _run(monkeypatch, tmp_path, range(2015, 2301), ctrl_gradient=False)
    line = _core_line(capsys, "core  9")
    assert "BLOCKED" in line and "dacabfdz for CESM2-WACCM ctrl" in line, line


# --- core 11 asks for the OCX product the driver will insist on -------------

def _ocx_tree(root, years):
    from icepack2_tools.forcing import OCX_ATMOSPHERE_SOURCE as SRC
    d = os.path.join(root, "OCX", SRC, "SDBN1-8000m", "acabf", "v1")
    os.makedirs(d, exist_ok=True)
    for y in years:
        open(os.path.join(d, f"acabf_AIS_{SRC}_OCX_SDBN1-8000m_v1_{y}.nc"), "wb").close()
    d = os.path.join(root, "OCX", "ocean", "main", "v1")
    os.makedirs(d, exist_ok=True)
    for var in ("tf", "so"):
        open(os.path.join(d, f"{var}_AIS_OCX_ocean_main_v1_1950-2025.nc"), "wb").close()


def _ocx_gradient(root, years, version):
    from icepack2_tools.forcing import OCX_ATMOSPHERE_SOURCE as SRC
    d = os.path.join(root, "OCX", SRC, "SDBN1-8000m", "dacabfdz", version)
    os.makedirs(d, exist_ok=True)
    for y in years:
        open(os.path.join(d, f"dacabfdz_AIS_{SRC}_OCX_SDBN1-8000m_{version}_{y}.nc"),
             "wb").close()


def _run_core_11(monkeypatch, tmp_path, years, forcing="protocol", env=None,
                 dhdt=False, gradient_versions=("v2",)):
    r"""preflight.main() over an OCX tree, with the ``dacabfdz`` of the
    SMB-elevation feedback at each of `gradient_versions` over 1979-2025."""
    root = str(tmp_path / "ISMIP7" / "AIS")
    os.makedirs(root, exist_ok=True)
    if years is not None:
        _ocx_tree(root, years)
    for version in gradient_versions:
        _ocx_gradient(root, range(1979, 2026), version)
    monkeypatch.delenv("ISMIP7_SMB_ELEVATION_FEEDBACK", raising=False)
    monkeypatch.setenv("ISMIP7_DATA_ROOT", root)
    monkeypatch.setenv("ISMIP7_OCX_FORCING", forcing)
    monkeypatch.delenv("ISMIP7_OCX_OCEAN", raising=False)
    _knobs(monkeypatch, tmp_path, env)
    sys.modules.pop("preflight", None)
    preflight = importlib.import_module("preflight")
    monkeypatch.setattr(preflight, "shared_missing", lambda warn: [])
    monkeypatch.setattr(preflight, "racmo_ok", lambda: True)
    monkeypatch.setattr(preflight, "oi_ok", lambda: True)
    monkeypatch.setattr(preflight, "pool_status", lambda *a, **k: ("ok", ""))
    if not dhdt:
        monkeypatch.setattr(preflight, "dhdt_missing", lambda: [])
    preflight.main()


def test_core_11_is_blocked_until_the_ocx_product_is_on_disk(monkeypatch, tmp_path, capsys):
    r"""RACMO and the OI climatology being present used to be enough, which is
    how the core ran on the stopgap with the product sitting unread."""
    _run_core_11(monkeypatch, tmp_path, None)
    line = _core_line(capsys, "core 11")
    assert "BLOCKED" in line and "OCX atmosphere acabf" in line and "absent" in line


def test_core_11_is_ready_on_the_product_and_points_at_the_tripwire(monkeypatch, tmp_path, capsys):
    r"""OCX starts in 2003 (issue 117), so a product holding 2003-2025 is
    all it reads."""
    _run_core_11(monkeypatch, tmp_path, range(2003, 2026))
    line = _core_line(capsys, "core 11")
    assert "READY" in line and "check_melt_bound.py --ocx" in line and "#48" in line
    assert "2003-2025" in line


def test_core_11_on_the_stopgap_says_that_is_what_it_is(monkeypatch, tmp_path, capsys):
    _run_core_11(monkeypatch, tmp_path, None, forcing="stopgap")
    line = _core_line(capsys, "core 11")
    assert "READY" in line and "ISMIP7_OCX_FORCING=stopgap" in line


def test_core_11_refuses_the_shifted_ocx_gradient(monkeypatch, tmp_path, capsys):
    r"""The OCX ``dacabfdz`` v1 is the spatially shifted file of discussion
    #45. With no v2 beside it the reader would fall back to it, so the gate
    blocks it, on the stopgap too, whose feedback reads the same gradient."""
    for forcing in ("protocol", "stopgap"):
        _run_core_11(monkeypatch, tmp_path, range(1979, 2026), forcing=forcing,
                     gradient_versions=("v1",))
        line = _core_line(capsys, "core 11")
        assert "BLOCKED" in line and "dacabfdz" in line, line
        assert "is v1 on disk" in line and "v2 or newer" in line


def test_core_11_reads_the_v2_gradient_beside_the_v1(monkeypatch, tmp_path, capsys):
    _run_core_11(monkeypatch, tmp_path, range(1979, 2026), gradient_versions=("v1", "v2"))
    line = _core_line(capsys, "core 11")
    assert "READY" in line, line


# --- a backdated cold start reads the Smith dH/dt ---------------------------

def _stage_kit(tmp_path):
    kit = tmp_path / "ISMIP7" / "AIS" / "obs" / "mipkit"
    kit.mkdir(parents=True)
    (kit / "AntarcticaObsISMIP7-v1.2.nc").touch()


def _stage_cache(tmp_path, *rasters):
    cache = tmp_path / "obs" / "dhdt_cache"
    cache.mkdir(parents=True)
    for raster in rasters:
        (cache / f"dhdt_smith_AntarcticaObsISMIP7-v1.2_{raster}.tif").touch()


def _historical_line(monkeypatch, tmp_path, capsys):
    r"""Core 1 from 2003 on a complete tree, with the real dH/dt lookup."""
    _run(monkeypatch, tmp_path, range(2003, 2015), scenario="historical",
         ocean_cover=lambda e, s: (2003, 2014), dhdt=True)
    return _core_line(capsys, "core  1")


def test_a_backdated_start_without_the_dhdt_is_blocked(monkeypatch, tmp_path, capsys):
    r"""setup_model would stop there, after loading the mesh and the MAP."""
    line = _historical_line(monkeypatch, tmp_path, capsys)
    assert "BLOCKED" in line and "dH/dt for the geometry backdating" in line, line
    assert "MIPkit not found" in line and "_valid.tif" in line, line


def test_the_mipkit_is_a_dhdt_source(monkeypatch, tmp_path, capsys):
    _stage_kit(tmp_path)
    line = _historical_line(monkeypatch, tmp_path, capsys)
    assert "READY" in line, line


def test_the_cached_rasters_are_a_dhdt_source(monkeypatch, tmp_path, capsys):
    r"""A cluster may stage the two small rasters and no 9 GB kit."""
    _stage_cache(tmp_path, "value", "valid")
    line = _historical_line(monkeypatch, tmp_path, capsys)
    assert "READY" in line, line


def test_a_cached_value_raster_without_its_coverage_is_no_source(monkeypatch, tmp_path, capsys):
    _stage_cache(tmp_path, "value")
    line = _historical_line(monkeypatch, tmp_path, capsys)
    assert "BLOCKED" in line and "dH/dt" in line, line


def test_a_start_from_2015_reads_no_dhdt(monkeypatch, tmp_path, capsys):
    r"""Nothing is backdated from 2015 on, so OCX started there needs no
    source."""
    _run_core_11(monkeypatch, tmp_path, range(2015, 2026),
                 env={"ISMIP7_T_START": "2015"}, dhdt=True)
    line = _core_line(capsys, "core 11")
    assert "READY" in line and "2015-2025" in line, line


def test_the_dhdt_lookup_writes_nothing_and_honours_a_named_kit(monkeypatch, tmp_path):
    r"""The preflight must not build the cache it asks about, and a named
    ISMIP7_OBS_KIT that does not exist stops the loader even beside a cache."""
    from icepack2_tools.obs_dhdt import dhdt_source
    _knobs(monkeypatch, tmp_path)
    monkeypatch.setenv("ISMIP7_DATA_ROOT", str(tmp_path / "ISMIP7" / "AIS"))
    _stage_kit(tmp_path)
    kit, val_fn, cov_fn = dhdt_source()
    assert kit.endswith("AntarcticaObsISMIP7-v1.2.nc")
    assert val_fn.endswith("dhdt_smith_AntarcticaObsISMIP7-v1.2_value.tif")
    assert cov_fn.endswith("dhdt_smith_AntarcticaObsISMIP7-v1.2_valid.tif")
    assert not (tmp_path / "obs").exists()
    _stage_cache(tmp_path, "value", "valid")
    monkeypatch.setenv("ISMIP7_OBS_KIT", str(tmp_path / "moved.nc"))
    with pytest.raises(FileNotFoundError, match="ISMIP7_OBS_KIT"):
        dhdt_source()


def test_the_preflight_imports_without_firedrake():
    r"""The preflight imports in 0.1 s because nothing it imports loads
    Firedrake, which is why obs_dhdt imports Firedrake inside the loader."""
    import subprocess
    code = ("import sys; sys.path[:0] = sys.argv[1:]; import preflight; "
            "loaded = [m for m in sys.modules if m.split('.')[0] == 'firedrake']; "
            "assert not loaded, loaded")
    result = subprocess.run(
        [sys.executable, "-c", code, REPO, os.path.join(REPO, "antarctica", "scripts")],
        capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


# ── the melt calibration a run reads ────────────────────────────────────────
def _shared(monkeypatch, tmp_path):
    r"""preflight.shared_missing under the environment the test set."""
    sys.modules.pop("preflight", None)
    preflight = importlib.import_module("preflight")
    return preflight.shared_missing([])


def test_a_deltat_file_is_the_melt_calibration(monkeypatch, tmp_path):
    r"""With ISMIP7_DELTAT_PER_BASIN_NPZ naming a file the drivers read no
    per-basin K, so the gate must not ask for one."""
    for knob in ("ISMIP7_K_PER_BASIN_NPZ", "ISMIP7_K_SCALE"):
        monkeypatch.delenv(knob, raising=False)
    dT = tmp_path / "deltaT_per_basin_1000_K8.500e-05.npz"
    np.savez(dT, K=8.5e-5)
    monkeypatch.setenv("ISMIP7_DELTAT_PER_BASIN_NPZ", str(dT))
    assert not [m for m in _shared(monkeypatch, tmp_path) if "per-basin" in m]


def test_a_named_deltat_file_that_is_absent_is_reported(monkeypatch, tmp_path):
    monkeypatch.delenv("ISMIP7_K_SCALE", raising=False)
    monkeypatch.setenv("ISMIP7_DELTAT_PER_BASIN_NPZ", str(tmp_path / "nope.npz"))
    missing = _shared(monkeypatch, tmp_path)
    assert any("ISMIP7_DELTAT_PER_BASIN_NPZ" in m and "does not exist" in m
               for m in missing)
