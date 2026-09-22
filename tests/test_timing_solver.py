r"""The solver a timing campaign's lanes time is a parameter; the caches' is not.

The matrix was built around condensed MUMPS, and `scpc_mumps` was spelled into
the tag, the record check, the cache check and every export. A campaign under
`scpc_gamg` has to be a different campaign (its own tag, so no record of one
ever passes as the other's) that starts every lane from the very state the
MUMPS lane started from: the caches stay `scpc_mumps`'s, files and fingerprint
both, and only the lane's own solver changes.
"""
import os
import shlex
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "antarctica"
sys.path.insert(0, str(ROOT / "scripts"))

import timing_campaign as tc  # noqa: E402
from icepack2_tools.solverconfig import solver_provenance  # noqa: E402

GAMG_TAG = "scpc_gamg_10step_dt0p125at2500_dg0_logvelnet_cached_strict_v4"
MUMPS_TAG = "scpc_mumps_10step_dt0p125at2500_dg0_logvelnet_cached_strict_v4"


def test_the_module_s_own_selftest_passes():
    tc.selftest()


def test_the_lane_solver_leads_the_tag_and_the_default_campaign_keeps_its_name(monkeypatch):
    monkeypatch.delenv("ISMIP7_TIMING_SOLVER", raising=False)
    assert tc.campaign_tag(10, 0.125, "strict") == MUMPS_TAG
    assert tc.campaign_tag(10, 0.125, "strict", "scpc_gamg") == GAMG_TAG
    monkeypatch.setenv("ISMIP7_TIMING_SOLVER", "scpc_gamg")
    assert tc.campaign_tag(10, 0.125, "strict") == GAMG_TAG
    assert tc.parse_campaign_tag(GAMG_TAG + "_probe") == {
        "solver": "scpc_gamg", "steps": 10, "dt_2500": 0.125, "contract": "strict",
        "lane": "probe", "version": 4}
    # The v3 archive, written before the solver was a parameter, still parses.
    assert tc.parse_campaign_tag(MUMPS_TAG[:-1] + "3")["solver"] == "scpc_mumps"


@pytest.mark.parametrize("solver", ["full_mumps", "schur_gamg", "iterative", ""])
def test_only_a_condensed_solver_names_a_campaign(solver):
    r"""The caches are condensed-MUMPS states and the lanes share one record
    contract; the reference and the legacy selfp modes are qualified and timed
    through `make reference` / `make qualify`, never as a matrix."""
    with pytest.raises(ValueError, match="timing solver must be one of"):
        tc.campaign_tag(10, 0.125, "strict", solver)
    with pytest.raises(ValueError, match="does not name a cached-strict lane"):
        tc.parse_campaign_tag(MUMPS_TAG.replace("scpc_mumps", solver or "x", 1))


def test_a_record_passes_only_under_the_solver_its_own_tag_names():
    gamg = tc.synthetic_record(2500, 25000, 16, GAMG_TAG)
    mumps = tc.synthetic_record(2500, 25000, 16, MUMPS_TAG)
    assert gamg["diagnostic_solver_mode"] == "scpc_gamg"
    assert tc.validate_timing_record(gamg, timing_tag=GAMG_TAG)[0]
    assert tc.validate_timing_record(mumps, timing_tag=MUMPS_TAG)[0]
    # Neither campaign's report accepts the other's record...
    assert not tc.validate_timing_record(gamg, timing_tag=MUMPS_TAG)[0]
    assert not tc.validate_timing_record(mumps, timing_tag=GAMG_TAG)[0]
    # ...nor a lane that ran another solver than its tag says.
    valid, detail = tc.validate_timing_record(
        dict(gamg, diagnostic_solver_mode="scpc_mumps"), timing_tag=GAMG_TAG)
    assert not valid and "tag names 'scpc_gamg'" in detail


def test_the_caches_are_the_prepare_solver_s_under_every_campaign(monkeypatch):
    monkeypatch.setenv("ISMIP7_TIMING_SOLVER", "scpc_gamg")
    assert tc.CACHE_SOLVER_MODE == "scpc_mumps"
    assert tc.cache_stem(2500, 25000) == (
        "initial_state_scpc_mumps_dg0_logvelnet_v4_2500_25000_buffered20000")


def manifest(tmp_path, fingerprint):
    cache = tmp_path / f"{tc.cache_stem(2500, 25000)}.h5"
    cache.write_text("x")
    return cache, {
        "cache_schema_version": tc.CACHE_SCHEMA_VERSION,
        "cache_role": tc.CACHE_ROLE,
        "lc": 2500, "lc_coarse": 25000, "buffer_m": tc.BUFFER_M,
        "diagnostic_solver_mode": "scpc_mumps",
        "friction": "budd", "friction_gate": tc.BUDD_SHELF_GATE,
        "geometry_space": "dg0", "n_flow": 3.0, "a4_factor": 1.0,
        "t_yr": tc.MATRIX_T_START,
        "mesh_basename": tc.mesh_basename(2500, 25000),
        "geometry_source_method": tc.TARGET_MESH_GEOMETRY_METHOD,
        "source_inversion_basename": tc.SOURCE_INVERSION_BASENAME,
        "source_inversion_sha256": "a", "source_mesh_sha256": "b",
        "geometry_source": "/data/bed.nc", "geometry_source_basename": "bed.nc",
        "checkpoint_fields": list(tc.CACHE_REQUIRED_FIELDS),
        "cache_path": str(cache),
        "solver_configuration_fingerprint": fingerprint,
    }


def test_a_gamg_lane_accepts_the_cache_the_mumps_campaign_prepared(tmp_path, monkeypatch):
    r"""What blocked the campaign: a lane fingerprinted its OWN PETSc options
    and held them against the manifest, so every GAMG lane found the cache
    "stale". It is the solver the cache was prepared under that is compared, in
    the lane's environment."""
    monkeypatch.setenv("ISMIP7_DIAGNOSTIC_LINEAR_SOLVER", "scpc_mumps")
    prepared = tc.solver_configuration_fingerprint(solver_provenance())
    cache, doc = manifest(tmp_path, prepared)

    monkeypatch.setenv("ISMIP7_DIAGNOSTIC_LINEAR_SOLVER", "scpc_gamg")
    own = solver_provenance()
    assert own["diagnostic_mode"] == "scpc_gamg"
    assert own["diagnostic_petsc_options"]["condensed_field_pc_type"] == "gamg"
    as_prepared = tc.solver_configuration_fingerprint(
        solver_provenance(tc.CACHE_SOLVER_MODE))
    assert as_prepared == prepared != tc.solver_configuration_fingerprint(own)
    assert tc.validate_cache_manifest(
        doc, lc=2500, lc_coarse=25000, cache_path=cache,
        solver_fingerprint=as_prepared) == (True, "cache provenance matches")

    # A lane whose nonlinear tolerances differ from the prepare's is still stale.
    monkeypatch.setenv("ISMIP7_SNES_RTOL", "1e-6")
    loosened = tc.solver_configuration_fingerprint(solver_provenance(tc.CACHE_SOLVER_MODE))
    assert tc.validate_cache_manifest(
        doc, lc=2500, lc_coarse=25000, cache_path=cache,
        solver_fingerprint=loosened) == (False, "cache solver configuration is stale")
    # And a cache that was itself prepared under GAMG is not a campaign cache.
    assert not tc.validate_cache_manifest(
        dict(doc, diagnostic_solver_mode="scpc_gamg"), lc=2500, lc_coarse=25000)[0]


def dry_run(tmp_path, stage, *args):
    env = dict(os.environ, ISMIP7_SITE="iu_quartz", ISMIP7_LOCAL_ENV=os.devnull,
               ISMIP7_FIREDRAKE=os.devnull,
               # An ambient forward-run solver must not reach a campaign.
               ISMIP7_DIAGNOSTIC_LINEAR_SOLVER="full_mumps")
    env.pop("ISMIP7_TIMING_SOLVER", None)
    proc = subprocess.run(
        [sys.executable, str(ROOT / "scripts/manage_timing_campaign.py"), stage,
         "--dry-run", "--assume-valid-caches", "--only-mesh", "2500/25000",
         "--cache-dir", str(tmp_path / "cache"), "--timing-dir", str(tmp_path / "timing"),
         *args],
        env=env, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    line = next(line for line in proc.stdout.splitlines() if line.startswith("DRY RUN:"))
    # submit.sh prints the command shell-quoted: `env KEY=VALUE ... sbatch ...
    # --export=ALL script`, every value set in sbatch's own environment.
    words = shlex.split(line.split(": ", 2)[2])
    assert words[0] == "env" and "--export=ALL" in words, words
    return dict(word.split("=", 1) for word in words[1:words.index("sbatch")])


@pytest.mark.parametrize("stage, kind, suffix", [
    ("scout", "matrix", ""), ("probe", "cache_probe", "_probe")])
def test_the_manager_runs_a_gamg_lane_from_the_mumps_cache(tmp_path, stage, kind, suffix):
    pytest.importorskip("firedrake")
    exports = dry_run(tmp_path, stage, "--solver", "scpc_gamg")
    assert exports["ISMIP7_DIAGNOSTIC_LINEAR_SOLVER"] == "scpc_gamg"
    assert exports["ISMIP7_TIMING_SOLVER"] == "scpc_gamg"
    assert exports["ISMIP7_TIMING_KIND"] == kind
    assert exports["ISMIP7_TIMING_TAG"] == GAMG_TAG + suffix
    assert Path(exports["ISMIP7_TIMING_STATUS"]).name == (
        f"status_{GAMG_TAG}{suffix}_2500_25000_16.txt")
    assert Path(exports["ISMIP7_RESTART"]).name.startswith(
        "initial_state_scpc_mumps_dg0_logvelnet_v4_2500_25000_")
    assert exports["ISMIP7_RESCUE_ENABLED"] == "0" and exports["ISMIP7_SUBCYCLES"] == "1"


def test_the_default_campaign_is_unchanged(tmp_path):
    pytest.importorskip("firedrake")
    exports = dry_run(tmp_path, "scout")
    assert exports["ISMIP7_DIAGNOSTIC_LINEAR_SOLVER"] == "scpc_mumps"
    assert exports["ISMIP7_TIMING_TAG"] == MUMPS_TAG


def test_cache_work_stays_under_the_prepare_solver_whatever_the_campaign(tmp_path):
    r"""`make timing TIMING_SOLVER=scpc_gamg` still walks through prepare; a
    cache it had to rebuild must come out as the MUMPS campaign's would."""
    pytest.importorskip("firedrake")
    inputs = [ROOT / "mesh" / tc.mesh_basename(2000, 20000),
              ROOT / "mesh" / f"boundary_ids_antarctica_20000_2000_buffered{tc.BUFFER_M}.json"]
    if not all(path.is_file() for path in inputs):
        pytest.skip("prepare looks for the 2000/20000 mesh and its boundary sidecar")
    source = tmp_path / tc.SOURCE_INVERSION_BASENAME
    source.write_text("stands in for the campaign source MAP\n")
    env = dict(os.environ, ISMIP7_SITE="iu_quartz", ISMIP7_LOCAL_ENV=os.devnull,
               ISMIP7_FIREDRAKE=os.devnull, ISMIP7_TIMING_SOLVER="scpc_gamg")
    proc = subprocess.run(
        [sys.executable, str(ROOT / "scripts/manage_timing_campaign.py"), "prepare",
         "--dry-run", "--force", "--only-mesh", "2000/20000", "--inversion", str(source),
         "--cache-dir", str(tmp_path / "cache"), "--timing-dir", str(tmp_path / "timing")],
        env=env, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    submitted = [line for line in proc.stdout.splitlines() if line.startswith("DRY RUN:")]
    assert submitted, proc.stdout
    for line in submitted:
        assert "ISMIP7_DIAGNOSTIC_LINEAR_SOLVER=scpc_mumps" in line
        assert "scpc_gamg" not in line


def test_an_unnamed_solver_is_still_the_full_jacobian_reference(monkeypatch):
    r"""Production forwards run scpc_gamg because projection.sbatch says so, not
    because the default moved. The inversion resolves the same default to stamp
    `diagnostic_solver_mode` on its MAP, and `redistribute_checkpoint.py`
    fingerprints a published cache with it, so moving it would relabel every
    inversion run outside the campaign manager."""
    from icepack2_tools import solverconfig

    monkeypatch.delenv("ISMIP7_DIAGNOSTIC_LINEAR_SOLVER", raising=False)
    assert solverconfig.DIAGNOSTIC_SOLVER_DEFAULT == "full_mumps"
    assert solverconfig.diagnostic_solver_mode() == "full_mumps"
