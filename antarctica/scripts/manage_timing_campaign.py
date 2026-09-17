#!/usr/bin/env python3
"""Advance the cached timing campaign by one idempotent Slurm stage."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.fspath(Path(__file__).resolve().parents[2]))

from icepack2_tools.solverconfig import (
    SNES_DIVERGENCE_TOL_DEFAULT,
    solver_provenance,
)
from timing_campaign import (
    BUFFER_M,
    CACHE_TAG,
    CONTRACTS,
    LANE_INITIAL_STATE_DEFAULT,
    MATRIX_DT_MAX,
    TRIPWIRE_DEFAULTS,
    MEMORY_BY_LC,
    SOLVER_MODE,
    SOURCE_INVERSION_BASENAME,
    atomic_write_status,
    cache_paths,
    expected_dt,
    expected_t_end,
    inversion_cores,
    inversion_maxiter,
    inversion_memory,
    inversion_required,
    lane_tag,
    campaign_tag,
    contract_exports,
    contract_name,
    matrix_dt_2500,
    matrix_steps,
    probe_tag,
    mesh_basename,
    mesh_inversion_basename,
    mesh_inversion_source_mesh,
    mesh_inversion_map_path,
    mesh_inversion_status_path,
    mesh_inversion_timing_json_path,
    mesh_rows,
    pristine_cache_paths,
    scaling_lanes,
    scout_lanes,
    sha256_file,
    solver_configuration_fingerprint,
    timing_record_basename,
    timing_status_basename,
    validate_cache_manifest,
    validate_timing_record,
)

_ROOT = Path(__file__).resolve().parents[1]
_LIVE_SLURM_STATES = {
    "CONFIGURING",
    "COMPLETING",
    "PENDING",
    "REQUEUED",
    "REQUEUE_FED",
    "REQUEUE_HOLD",
    "RESIZING",
    "RUNNING",
    "STAGE_OUT",
    "SUSPENDED",
}
_TERMINAL_SLURM_STATES = {
    "BOOT_FAIL",
    "CANCELLED",
    "COMPLETED",
    "DEADLINE",
    "FAILED",
    "NODE_FAIL",
    "OUT_OF_MEMORY",
    "PREEMPTED",
    "REVOKED",
    "SPECIAL_EXIT",
    "STOPPED",
    "TIMEOUT",
}


def read_status(path):
    try:
        tokens = Path(path).read_text().strip().split()
    except OSError:
        return None
    if not tokens:
        return {"state": "invalid"}
    parsed = {"state": tokens[0]}
    for token in tokens[1:]:
        if "=" in token:
            key, value = token.split("=", 1)
            parsed[key] = value
    return parsed


# submit.sh's exit status for a request one node of the site cannot hold.
_SUBMIT_NOT_RUNNABLE = 3


def _not_runnable_reason(message):
    """The single-token ``reason=`` submit.sh printed, for the status file."""
    for token in message.split():
        if token.startswith("reason="):
            return token.split("=", 1)[1]
    return "exceeds_site_node"


def read_record(path):
    try:
        with open(path) as stream:
            return json.load(stream), None
    except FileNotFoundError:
        return None, None
    except (OSError, json.JSONDecodeError) as exc:
        return None, str(exc)


def _slurm_outcome_from_outputs(queue_output, accounting_output):
    """Classify one Slurm allocation from machine-readable command output."""
    if queue_output.strip():
        return "active", None, None
    rows = [line for line in accounting_output.splitlines() if line.strip()]
    if not rows:
        return "unknown", None, None
    fields = rows[0].split("|")
    raw_state = fields[0].strip()
    state = raw_state.split()[0].rstrip("+").upper() if raw_state else ""
    exit_code = fields[1].strip() if len(fields) > 1 else ""
    if state in _LIVE_SLURM_STATES:
        return "active", state, exit_code
    if state in _TERMINAL_SLURM_STATES:
        return "terminal", state, exit_code
    return "unknown", state or None, exit_code or None


def _slurm_outcome(status):
    """Return ``(kind, Slurm state, exit code)`` for an active-looking stamp.

    Unknown scheduler state is treated conservatively by callers: never submit
    a duplicate merely because accounting is delayed or unavailable.
    """
    if not status or status.get("state") not in {
        "submitting", "submitted", "running"
    }:
        return "not_active_stamp", None, None
    job_id = status.get("job_id")
    if not job_id or shutil.which("squeue") is None:
        return "unknown", None, None
    queue = subprocess.run(
        ["squeue", "-h", "-j", job_id, "-o", "%i"],
        check=False,
        capture_output=True,
        text=True,
    )
    if queue.returncode == 0 and queue.stdout.strip():
        return "active", None, None
    # A completed job may already be purged from squeue, which can be reported
    # as either an empty successful query or an "invalid job id" error.  In
    # both cases accounting, not the squeue return code, is authoritative.
    if shutil.which("sacct") is None:
        return "unknown", None, None
    accounting = subprocess.run(
        [
            "sacct", "-n", "-X", "-j", job_id,
            "-o", "State,ExitCode", "-P",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if accounting.returncode != 0:
        return "unknown", None, None
    return _slurm_outcome_from_outputs(queue.stdout, accounting.stdout)


def _timestamp():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def scaling_decision(scout_state):
    if scout_state == "passed":
        return "submit"
    if scout_state == "failed":
        return "block"
    return "wait"


def _parse_mesh(value):
    parts = value.replace(":", "/").split("/")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            "mesh must be written as LC/LC_coarse, for example 2500/25000"
        )
    try:
        pair = tuple(int(part) for part in parts)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("mesh values must be integers") from exc
    if pair not in mesh_rows():
        choices = ", ".join(f"{lc}/{lcc}" for lc, lcc in mesh_rows())
        raise argparse.ArgumentTypeError(
            f"mesh {value!r} is not in the timing policy; choose {choices}"
        )
    return pair


class CampaignManager:
    def __init__(self, args):
        self.root = Path(args.root).resolve()
        self.mesh_dir = self.root / "mesh"
        self.cache_dir = Path(args.cache_dir).resolve()
        self.timing_dir = Path(args.timing_dir).resolve()
        self.logs_dir = self.root / "results" / "logs"
        self.inversion = Path(args.inversion).resolve()
        self.queue = args.queue
        self.partition = args.partition
        self.constraint = args.constraint
        self.walltime = args.walltime
        self.maxiter = inversion_maxiter(args.maxiter)
        self.dry_run = args.dry_run
        self.force = args.force
        self.assume_valid_caches = args.assume_valid_caches
        self.follow_prepare = bool(getattr(args, "follow_prepare", False))
        self.follow_invert = bool(getattr(args, "follow_invert", False))
        self.only_mesh = args.only_mesh
        self.monitor = args.monitor
        # The interval (steps, 2.5 km dt) and the physics contract are
        # campaign parameters: they name the tag, so records from different
        # ladder rungs never mix, and they are pinned into this process's
        # environment so every path helper and validator agrees.
        self.matrix_steps = matrix_steps(getattr(args, "matrix_steps", None))
        self.matrix_dt_2500 = matrix_dt_2500(
            getattr(args, "matrix_dt_2500", None)
        )
        self.contract = contract_name(getattr(args, "contract", None))
        os.environ["ISMIP7_MATRIX_STEPS"] = str(self.matrix_steps)
        os.environ["ISMIP7_MATRIX_DT_2500"] = f"{self.matrix_dt_2500:g}"
        os.environ["ISMIP7_TIMING_CONTRACT"] = self.contract
        self.campaign_tag = campaign_tag(
            self.matrix_steps, self.matrix_dt_2500, self.contract
        )
        self.probe_tag = probe_tag(
            self.matrix_steps, self.matrix_dt_2500, self.contract
        )
        self.initial_state = getattr(
            args, "initial_state", LANE_INITIAL_STATE_DEFAULT
        )
        self.lane_tag = lane_tag(self.initial_state, self.campaign_tag)
        self.submit_failures = 0
        self.source_sha256 = None
        self.mesh_checksums = {}
        self.solver_configuration = solver_provenance()
        self.solver_fingerprint = solver_configuration_fingerprint(
            self.solver_configuration
        )
        # Keep path helpers and validate_cache_manifest in lockstep with the
        # CLI/Makefile maxiter for this process.
        os.environ["ISMIP7_TIMING_INVERSION_MAXITER"] = str(self.maxiter)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.timing_dir.mkdir(parents=True, exist_ok=True)
        self.logs_dir.mkdir(parents=True, exist_ok=True)

    def selected_rows(self):
        rows = mesh_rows()
        if self.only_mesh is None:
            return rows
        return tuple(row for row in rows if row == self.only_mesh)

    def selected_lanes(self, lanes):
        if self.only_mesh is None:
            return lanes
        return tuple(lane for lane in lanes if lane[:2] == self.only_mesh)

    def reconcile_status(self, status_path):
        """Turn a dead Slurm allocation's live-looking stamp into a failure."""
        status = read_status(status_path)
        kind, slurm_state, exit_code = _slurm_outcome(status)
        if kind in {"active", "unknown"}:
            return status, True
        if kind != "terminal":
            return status, False
        category = (
            "incomplete_output"
            if slurm_state == "COMPLETED"
            else "external_termination"
        )
        reconciled = {
            "state": "failed",
            "category": category,
            "phase": "scheduler",
            "slurm_state": slurm_state,
            "exit_code": exit_code,
            "job_id": status.get("job_id"),
            "timestamp": _timestamp(),
        }
        if not self.dry_run:
            atomic_write_status(
                status_path,
                reconciled.pop("state"),
                **reconciled,
            )
            reconciled["state"] = "failed"
        return reconciled, False

    def source_checksum(self):
        if self.source_sha256 is None:
            if not self.inversion.is_file():
                raise FileNotFoundError(self.inversion)
            self.source_sha256 = sha256_file(self.inversion)
        return self.source_sha256

    def mesh_path(self, lc, lc_coarse):
        return self.mesh_dir / mesh_basename(lc, lc_coarse)

    def mesh_checksum(self, lc, lc_coarse):
        key = (lc, lc_coarse)
        if key not in self.mesh_checksums:
            path = self.mesh_path(lc, lc_coarse)
            if not path.is_file():
                raise FileNotFoundError(path)
            self.mesh_checksums[key] = sha256_file(path)
        return self.mesh_checksums[key]

    def is_campaign_source_mesh(self, lc, lc_coarse):
        """Whether the campaign source MAP is this mesh's own invert output.

        Since the 2026-09-17 promotion the source is the 2500/25000
        re-inversion: an invert on that mesh (at any --maxiter) would
        republish the source cache, and a prepare there is a same-mesh DG0
        transfer simulation.py refuses. Its cache is the one the invert job
        published from the MAP itself. Recognised by basename, so an
        overridden TIMING_INVERSION that is not a per-mesh MAP has no
        source mesh.
        """
        return mesh_inversion_source_mesh(self.inversion.name) == (
            int(lc), int(lc_coarse)
        )

    def boundary_path(self, lc, lc_coarse):
        return self.mesh_dir / (
            f"boundary_ids_antarctica_{lc_coarse}_{lc}_buffered{BUFFER_M}.json"
        )

    def cache_status_path(self, lc, lc_coarse):
        return self.timing_dir / (
            f"status_cache_{CACHE_TAG}_{lc}_{lc_coarse}.txt"
        )

    def cache_audit_paths(self, lc, lc_coarse):
        stem = f"cache_audit_{CACHE_TAG}_{lc}_{lc_coarse}"
        return (
            self.timing_dir / f"{stem}.json",
            self.timing_dir / f"status_{stem}.txt",
        )

    def cache_validation(
        self,
        lc,
        lc_coarse,
        *,
        require_mesh_inversion=False,
        require_source=False,
    ):
        """``require_mesh_inversion``: only the per-mesh invert MAP may be the
        cache's source; ``require_source``: only the imported source MAP may
        (a control lane must not quietly run from an invert-published cache)."""
        if self.dry_run and self.assume_valid_caches:
            return True, "assumed valid for dry-run command inspection"
        cache, manifest_path = cache_paths(
            self.cache_dir, lc, lc_coarse
        )
        try:
            with open(manifest_path) as stream:
                manifest = json.load(stream)
        except (OSError, json.JSONDecodeError) as exc:
            return False, f"cache manifest unavailable: {exc}"
        try:
            mesh_sha256 = self.mesh_checksum(lc, lc_coarse)
        except FileNotFoundError as exc:
            return False, f"cache source unavailable: {exc}"

        mesh_map = mesh_inversion_map_path(
            self.root, lc, lc_coarse, maxiter=self.maxiter
        )
        candidates = []
        if mesh_map.is_file() and not require_source:
            candidates.append(
                (
                    sha256_file(mesh_map),
                    mesh_inversion_basename(
                        lc, lc_coarse, maxiter=self.maxiter
                    ),
                )
            )
        if not require_mesh_inversion:
            try:
                candidates.append(
                    (self.source_checksum(), SOURCE_INVERSION_BASENAME)
                )
            except FileNotFoundError as exc:
                if not candidates:
                    return False, f"cache source unavailable: {exc}"

        details = []
        for source_sha256, _basename in candidates:
            valid, detail = validate_cache_manifest(
                manifest,
                lc=lc,
                lc_coarse=lc_coarse,
                cache_path=cache,
                source_sha256=source_sha256,
                mesh_sha256=mesh_sha256,
                solver_fingerprint=self.solver_fingerprint,
            )
            if valid:
                return True, detail
            details.append(detail)
        if require_mesh_inversion and not mesh_map.is_file():
            return False, f"mesh inversion MAP missing: {mesh_map}"
        if require_source and details:
            return False, (
                f"{details[0]} (control lanes need the cache published by "
                "make timing-prepare, not by an invert)"
            )
        return False, details[0] if details else "cache provenance mismatch"

    def warm_start_validation(self, lc, lc_coarse):
        """Validate the invert's warm start.

        Prefer the pristine prepare copy: its manifest still names the
        imported source MAP after an earlier invert (say a maxiter=2 smoke)
        has republished the cache path, so a production invert never waits
        on a re-prepare merely because a smoke invert published first.
        """
        pristine, pristine_manifest = pristine_cache_paths(
            self.cache_dir, lc, lc_coarse
        )
        if not pristine.is_file() or (
            self.dry_run and self.assume_valid_caches
        ):
            return self.cache_validation(lc, lc_coarse)
        if not pristine_manifest.is_file():
            valid, detail = self.cache_validation(lc, lc_coarse)
            if valid:
                return True, detail
            return False, (
                f"{detail}; pristine copy {pristine.name} has no manifest "
                "(prepared before manifests were copied) -- re-run "
                "make timing-prepare FORCE_TIMING=1"
            )
        try:
            with open(pristine_manifest) as stream:
                manifest = json.load(stream)
        except (OSError, json.JSONDecodeError) as exc:
            return False, f"pristine manifest unavailable: {exc}"
        try:
            mesh_sha256 = self.mesh_checksum(lc, lc_coarse)
            source_sha256 = self.source_checksum()
        except FileNotFoundError as exc:
            return False, f"cache source unavailable: {exc}"
        valid, detail = validate_cache_manifest(
            manifest,
            lc=lc,
            lc_coarse=lc_coarse,
            cache_path=None,
            source_sha256=source_sha256,
            mesh_sha256=mesh_sha256,
            solver_fingerprint=self.solver_fingerprint,
        )
        return valid, f"pristine prepare copy: {detail}"

    def lane_cache(self, lc, lc_coarse):
        """``(cache, manifest, detail)`` for a lane from the transferred state.

        Prefer the pristine prepare copy (its manifest keeps naming the
        source MAP after an invert has republished the cache path); fall
        back to the published cache, which must then descend from the
        source MAP itself. ``cache`` is None when neither validates.
        """
        pristine, pristine_manifest = pristine_cache_paths(
            self.cache_dir, lc, lc_coarse
        )
        if (
            pristine.is_file()
            and pristine_manifest.is_file()
            and not (self.dry_run and self.assume_valid_caches)
        ):
            valid, detail = self.warm_start_validation(lc, lc_coarse)
            if valid:
                return pristine, pristine_manifest, detail
        valid, detail = self.cache_validation(
            lc, lc_coarse, require_source=True
        )
        cache, manifest = cache_paths(self.cache_dir, lc, lc_coarse)
        if not valid:
            return None, None, detail
        return cache, manifest, detail

    def lane_exports(self, ncores_unused=None):
        """Interval, contract and tripwire environment shared by every lane."""
        exports = {
            "ISMIP7_MATRIX_STEPS": self.matrix_steps,
            "ISMIP7_MATRIX_DT_2500": f"{self.matrix_dt_2500:g}",
            "ISMIP7_TIMING_CONTRACT": self.contract,
            # Print the budget line every step, not only at steps 1 and 5.
            "ISMIP7_OUTPUT_INTERVAL": "1",
        }
        exports.update(contract_exports(self.contract))
        exports.update(TRIPWIRE_DEFAULTS)
        return exports

    def inversion_map_path(self, lc, lc_coarse):
        return mesh_inversion_map_path(
            self.root, lc, lc_coarse, maxiter=self.maxiter
        )

    def inversion_paths(self, lc, lc_coarse):
        ncores = inversion_cores(lc)
        return (
            self.inversion_map_path(lc, lc_coarse),
            mesh_inversion_timing_json_path(
                self.root, lc, lc_coarse, ncores, maxiter=self.maxiter
            ),
            mesh_inversion_status_path(
                self.root, lc, lc_coarse, maxiter=self.maxiter
            ),
            ncores,
        )

    def inversion_result(self, lc, lc_coarse):
        map_path, timing_json, status_path, ncores = self.inversion_paths(
            lc, lc_coarse
        )
        if self.dry_run and self.assume_valid_caches:
            return (
                "passed",
                "assumed valid for dry-run command inspection",
                status_path,
            )
        if not inversion_required(lc):
            return (
                "passed",
                "skipped by policy: lanes use the transferred prepare cache",
                status_path,
            )
        status, active = self.reconcile_status(status_path)
        if active:
            return "active", status["state"], status_path
        if map_path.is_file() and timing_json.is_file():
            try:
                with open(timing_json) as stream:
                    record = json.load(stream)
            except (OSError, json.JSONDecodeError) as exc:
                return "failed", f"inversion timing JSON unreadable: {exc}", status_path
            if (
                int(record.get("lc", -1)) == int(lc)
                and int(record.get("lc_coarse", -1)) == int(lc_coarse)
                and int(record.get("maxiter", -1)) == int(self.maxiter)
                and int(record.get("ncores", -1)) == int(ncores)
                and os.path.realpath(str(record.get("map_path", "")))
                == os.path.realpath(map_path)
            ):
                # The invert's own verdict on its publishing solve comes
                # first: a controls-only MAP fails the cache check below too,
                # but "cache not republished" hides the cause.
                if "final_solve" not in record:
                    return (
                        "failed",
                        "inversion record predates the publish gate "
                        "(2026-09-14); its MAP may carry a mixed state from "
                        "other controls -- re-run make timing-inversion "
                        "FORCE_TIMING=1",
                        status_path,
                    )
                phase = record.get("phase")
                if phase != "finished":
                    return (
                        "failed",
                        f"inversion record stopped at phase={phase!r} "
                        "(job ended before the mixed state was published)",
                        status_path,
                    )
                final_solve = record.get("final_solve") or {}
                if final_solve and not final_solve.get("published", False):
                    return (
                        "failed",
                        "inversion final solve did not publish a mixed state: "
                        f"{final_solve.get('reason')} after "
                        f"{final_solve.get('iterations')} Newton iterations",
                        status_path,
                    )
                cache_ok, cache_detail = self.cache_validation(
                    lc, lc_coarse, require_mesh_inversion=True
                )
                if cache_ok:
                    return "passed", "inversion MAP and cache ready", status_path
                return (
                    "failed",
                    f"inversion MAP present but cache not republished: "
                    f"{cache_detail}",
                    status_path,
                )
            return "failed", "inversion timing JSON does not match mesh", status_path
        if status and status.get("state") in {
            "failed", "finished", "submission_failed", "not_runnable"
        }:
            if status.get("state") == "finished" and not map_path.is_file():
                return "failed", "status finished but MAP missing", status_path
            return (
                "failed",
                status.get("category", status["state"]),
                status_path,
            )
        return "missing", "no inversion record or active status", status_path

    def invert(self):
        for lc, lc_coarse in self.selected_rows():
            if not inversion_required(lc):
                print(
                    f"INVERSION SKIPPED {lc}/{lc_coarse}: policy -- lanes use "
                    "the transferred prepare cache"
                )
                continue
            if self.is_campaign_source_mesh(lc, lc_coarse):
                print(
                    f"INVERSION IS SOURCE {lc}/{lc_coarse}: its invert output "
                    f"{self.inversion_map_path(lc, lc_coarse).name} is the "
                    "campaign source MAP and is never overwritten (promote "
                    "another MAP by changing TIMING_INVERSION)"
                )
                continue
            state, detail, status_path = self.inversion_result(lc, lc_coarse)
            # A live allocation owns the MAP, cache and status paths; --force
            # retries a failed invert or replaces a passed one, it must not
            # race a running one.
            if state == "active":
                print(
                    f"INVERSION ACTIVE {lc}/{lc_coarse}: {detail}"
                    + ("; scancel it before FORCE_TIMING=1" if self.force else "")
                )
                continue
            # --assume-valid-caches makes inversion_result report "passed" so
            # scout/scale dry-runs can print lanes; invert must still emit its
            # own DRY RUN submit lines in that mode.
            if (
                state == "passed"
                and not self.force
                and not (self.dry_run and self.assume_valid_caches)
            ):
                print(f"INVERSION PASSED {lc}/{lc_coarse}: {detail}")
                continue
            if state == "failed" and not self.force:
                print(
                    f"INVERSION FAILED {lc}/{lc_coarse}: {detail}; retry with "
                    "FORCE_TIMING=1 after inspection"
                )
                continue
            # Prepare must have produced a usable imported-MAP cache first,
            # unless --follow-prepare queues the invert behind an active
            # prepare allocation via Slurm afterok.
            prep_ok, prep_detail = self.warm_start_validation(lc, lc_coarse)
            dependency = None
            if not prep_ok:
                if not self.follow_prepare:
                    print(
                        f"INVERSION WAITING CACHE {lc}/{lc_coarse}: "
                        f"{prep_detail}"
                    )
                    continue
                prep_status, prep_active = self.reconcile_status(
                    self.cache_status_path(lc, lc_coarse)
                )
                prep_job = (
                    prep_status.get("job_id") if prep_status else None
                )
                if prep_active and prep_job:
                    dependency = f"afterok:{prep_job}"
                    print(
                        f"INVERSION FOLLOW PREPARE {lc}/{lc_coarse}: "
                        f"dependency={dependency}"
                    )
                elif (
                    prep_status
                    and prep_status.get("state") == "finished"
                    and prep_job
                ):
                    # Stamp says finished but cache validation failed - do not
                    # afterok a dead success that left a bad cache.
                    print(
                        f"INVERSION WAITING CACHE {lc}/{lc_coarse}: "
                        f"{prep_detail} (prepare stamped finished)"
                    )
                    continue
                else:
                    print(
                        f"INVERSION WAITING CACHE {lc}/{lc_coarse}: "
                        f"{prep_detail}; no active prepare job to follow"
                    )
                    continue
            mesh = self.mesh_path(lc, lc_coarse)
            boundary = self.boundary_path(lc, lc_coarse)
            missing = [path for path in (mesh, boundary) if not path.is_file()]
            if missing and not (
                self.dry_run and self.assume_valid_caches
            ):
                message = ",".join(os.fspath(path) for path in missing)
                print(
                    f"INVERSION NOT RUNNABLE {lc}/{lc_coarse}: missing "
                    f"{message}"
                )
                if not self.dry_run:
                    atomic_write_status(
                        status_path,
                        "not_runnable",
                        reason="missing_input",
                    )
                continue
            map_path, timing_json, _, ncores = self.inversion_paths(
                lc, lc_coarse
            )
            cache, manifest = cache_paths(self.cache_dir, lc, lc_coarse)
            map_raw = map_path.with_name(map_path.stem + ".parallel.h5")
            # The prepare cache is the warm start; the imported MAP is not
            # read by the invert (controls/geometry/mixed state come from the
            # cache). After the invert, its full-state MAP is published AS the
            # timing cache (no second cold prepare), so warm-start from the
            # pristine copy the prepare job keeps beside it, never from a
            # cache this stage itself published.
            pristine, _ = pristine_cache_paths(self.cache_dir, lc, lc_coarse)
            if pristine.is_file():
                warm_start = pristine
            else:
                warm_start = cache
                if not (self.dry_run and self.assume_valid_caches):
                    print(
                        f"INVERSION WARM START {lc}/{lc_coarse}: no pristine "
                        f"prepare copy at {pristine}; using {cache} (if that "
                        "cache was published by an earlier invert, re-run "
                        "make timing-prepare FORCE_TIMING=1 first)"
                    )
            exports = {
                "ISMIP7_LC": lc,
                "ISMIP7_LC_COARSE": lc_coarse,
                "ISMIP7_BUFFER_M": BUFFER_M,
                "ISMIP7_MESH": mesh,
                "ISMIP7_BNDIDS": boundary,
                "ISMIP7_WARM_START": warm_start,
                "ISMIP7_SKIP_CONTINUATION": "1",
                "ISMIP7_MAP_OUT": map_path,
                "ISMIP7_MAP_OUT_RAW": map_raw,
                "ISMIP7_INVERSION_TIMING_JSON": timing_json,
                "ISMIP7_MAXITER": self.maxiter,
                "ISMIP7_FRICTION": "budd",
                "ISMIP7_GEOMETRY_SPACE": "dg0",
                "ISMIP7_N_FLOW": "3.0",
                "ISMIP7_A4_FACTOR": "1.0",
                "ISMIP7_MISFIT_NORM": "sigma",
                "ISMIP7_LOG_VEL_WEIGHT": "auto",
                "ISMIP7_LOG_VEL_EPS": "1.0",
                "ISMIP7_DHDT_WEIGHT": "1",
                "ISMIP7_DHDT_NET_SIGMA": "10",
                "ISMIP7_GAMMA_THETA": "1e5",
                "ISMIP7_GAMMA_PHI": "1e5",
                "ISMIP7_DIAGNOSTIC_LINEAR_SOLVER": SOLVER_MODE,
                "ISMIP7_SNES_DIVERGENCE_TOL": SNES_DIVERGENCE_TOL_DEFAULT,
                "ISMIP7_RESCUE_ENABLED": "1",
                "ISMIP7_TIMING_CACHE": cache,
                "ISMIP7_TIMING_CACHE_MANIFEST": manifest,
                "ISMIP7_TIMING_INVERSION_STATUS": status_path,
                "ISMIP7_TIMING_INVERSION_MAXITER": self.maxiter,
            }
            if self.monitor:
                exports.update({
                    "ISMIP7_SNES_MONITOR": "1",
                    "ISMIP7_SNES_LOG": self.logs_dir / (
                        f"inversion_snes_{CACHE_TAG}_{lc}_{lc_coarse}"
                        f"_{ncores}.log"
                    ),
                })
            self._submit(
                f"timing_inv_{lc}_{lc_coarse}",
                ncores,
                inversion_memory(lc),
                exports,
                self.root
                / "scripts/batch_runners/timing_matrix_inversion.script",
                status_path,
                dependency=dependency,
            )

    def _submit(
        self,
        job_name,
        ncores,
        memory,
        exports,
        script,
        status_path,
        dependency=None,
    ):
        # submit.sh composes the sbatch line from this cluster's site file
        # (account, node feature, extra flags, per-node limits), so a lane is
        # the same request here at every site. It submits from the checkout
        # this manager runs in, whatever the site file would default to.
        submit = self.root / "scripts/batch_runners/submit.sh"
        command = [
            "bash",
            os.fspath(submit),
            "script",
            os.fspath(Path(script).relative_to(self.root)),
            "--cd",
            self.root.name,
            "--name",
            job_name,
            "--tasks",
            str(ncores),
            "--mem",
            memory,
            "--time",
            self.walltime,
        ]
        if self.partition:
            command.extend(["--partition", self.partition])
        else:
            command.extend(["--queue", self.queue])
        if self.constraint:
            command.extend(["--constraint", self.constraint])
        if dependency:
            command.extend(["--dependency", dependency])
        if self.dry_run:
            command.append("--dry-run")
        command.extend(f"{key}={value}" for key, value in exports.items())
        env = dict(os.environ, ISMIP7_REPO=os.fspath(self.root.parent))
        if not self.dry_run:
            print("SUBMIT:", shlex.join(command))
            if shutil.which("sbatch") is None:
                raise RuntimeError(
                    "sbatch is unavailable; use --dry-run off-cluster"
                )
            atomic_write_status(
                status_path, "submitting", timestamp=_timestamp()
            )
        result = subprocess.run(
            command, check=False, capture_output=True, text=True, env=env
        )
        if (
            self.dry_run
            and result.returncode == 2
            and "no site definition matches" in result.stderr
        ):
            # Off-cluster there is no site to name. A dry run starts nothing,
            # so show the request as the no-scheduler site would compose it.
            env["ISMIP7_SITE"] = "local"
            result = subprocess.run(
                command, check=False, capture_output=True, text=True, env=env
            )
        composed = result.stderr.strip()
        if result.returncode == _SUBMIT_NOT_RUNNABLE:
            # One node of this site cannot hold the request. The lanes stay
            # single-node and comparable between sites, so this is an outcome
            # to record, not a failure to retry; --force does not change it.
            reason = _not_runnable_reason(composed)
            print("DRY RUN:" if self.dry_run else "NOT RUNNABLE:", composed)
            if not self.dry_run:
                atomic_write_status(
                    status_path,
                    "not_runnable",
                    reason=reason,
                    timestamp=_timestamp(),
                )
            return
        if self.dry_run:
            print("DRY RUN:", composed)
            if result.returncode != 0:
                self.submit_failures += 1
            return
        if composed:
            print(composed)
        if result.returncode != 0:
            atomic_write_status(
                status_path,
                "submission_failed",
                exit_code=result.returncode,
                timestamp=_timestamp(),
            )
            self.submit_failures += 1
            return
        job_id = result.stdout.strip().split(";", 1)[0]
        atomic_write_status(
            status_path, "submitted", job_id=job_id, timestamp=_timestamp()
        )

    def prepare(self):
        try:
            source_sha256 = self.source_checksum()
        except FileNotFoundError as exc:
            raise SystemExit(f"Cannot prepare caches: {exc}") from exc
        print(f"Source inversion sha256: {source_sha256}")
        for lc, lc_coarse in self.selected_rows():
            valid, detail = self.cache_validation(lc, lc_coarse)
            if self.is_campaign_source_mesh(lc, lc_coarse):
                # The source mesh's cache is the MAP itself, published by
                # the invert job; a prepare here would be a same-mesh DG0
                # transfer, which simulation.py refuses. --force cannot
                # rebuild it, only a republish from the MAP can.
                if valid:
                    print(f"CACHE OK {lc}/{lc_coarse}: {detail} (source mesh)")
                else:
                    cache, manifest = cache_paths(
                        self.cache_dir, lc, lc_coarse
                    )
                    print(
                        f"CACHE NOT RUNNABLE {lc}/{lc_coarse}: {detail}; "
                        "this is the campaign source MAP's own mesh, so "
                        "republish its cache from the MAP instead of "
                        "preparing: python scripts/redistribute_checkpoint.py "
                        f"--input {self.inversion} --output {cache} "
                        f"--manifest {manifest} --publish-timing-cache"
                    )
                continue
            if valid and not self.force:
                print(f"CACHE OK {lc}/{lc_coarse}: {detail}")
                continue
            status_path = self.cache_status_path(lc, lc_coarse)
            status, active = self.reconcile_status(status_path)
            if not self.force and active:
                print(f"CACHE ACTIVE {lc}/{lc_coarse}: {status['state']}")
                continue
            if not self.force and status and status.get("state") == "failed":
                print(
                    f"CACHE FAILED {lc}/{lc_coarse}: retry with "
                    "FORCE_TIMING=1 after inspection"
                )
                continue
            mesh = self.mesh_path(lc, lc_coarse)
            boundary = self.boundary_path(lc, lc_coarse)
            missing = [path for path in (mesh, boundary) if not path.is_file()]
            if missing:
                message = ",".join(os.fspath(path) for path in missing)
                print(f"CACHE NOT RUNNABLE {lc}/{lc_coarse}: missing {message}")
                if not self.dry_run:
                    atomic_write_status(
                        status_path, "not_runnable", reason="missing_input"
                    )
                continue
            cache, manifest = cache_paths(self.cache_dir, lc, lc_coarse)
            raw = cache.with_suffix(".parallel.h5")
            ncores = 32 if lc == 500 else 16
            exports = {
                "ISMIP7_LC": lc,
                "ISMIP7_LC_COARSE": lc_coarse,
                "ISMIP7_BUFFER_M": BUFFER_M,
                "ISMIP7_MESH": mesh,
                "ISMIP7_BNDIDS": boundary,
                "ISMIP7_INVERSION": self.inversion,
                "ISMIP7_FRICTION": "budd",
                "ISMIP7_GEOMETRY_SPACE": "dg0",
                "ISMIP7_N_FLOW": "3.0",
                "ISMIP7_A4_FACTOR": "1.0",
                "ISMIP7_DIAGNOSTIC_LINEAR_SOLVER": SOLVER_MODE,
                "ISMIP7_SNES_DIVERGENCE_TOL": SNES_DIVERGENCE_TOL_DEFAULT,
                "ISMIP7_RESCUE_ENABLED": "1",
                "ISMIP7_TIMING_CACHE_RAW": raw,
                "ISMIP7_TIMING_CACHE": cache,
                "ISMIP7_TIMING_CACHE_MANIFEST": manifest,
                "ISMIP7_TIMING_CACHE_PRISTINE": pristine_cache_paths(
                    self.cache_dir, lc, lc_coarse
                )[0],
                "ISMIP7_TIMING_CACHE_STATUS": status_path,
            }
            self._submit(
                f"timing_cache_{lc}_{lc_coarse}",
                ncores,
                MEMORY_BY_LC[lc],
                exports,
                self.root / "scripts/batch_runners/timing_prepare.script",
                status_path,
            )

    def audit(self):
        for lc, lc_coarse in self.selected_rows():
            valid, detail = self.cache_validation(lc, lc_coarse)
            if not valid:
                print(f"AUDIT WAITING CACHE {lc}/{lc_coarse}: {detail}")
                continue
            output_path, status_path = self.cache_audit_paths(lc, lc_coarse)
            if output_path.is_file() and not self.force:
                print(f"AUDIT OK {lc}/{lc_coarse}: {output_path}")
                continue
            status, active = self.reconcile_status(status_path)
            if active and not self.force:
                print(f"AUDIT ACTIVE {lc}/{lc_coarse}: {status['state']}")
                continue
            if not self.force and status and status.get("state") == "failed":
                print(
                    f"AUDIT FAILED {lc}/{lc_coarse}: retry with "
                    "FORCE_TIMING=1 after inspection"
                )
                continue
            cache, manifest = cache_paths(self.cache_dir, lc, lc_coarse)
            ncores = 32 if lc == 500 else 16
            exports = {
                "ISMIP7_TIMING_CACHE": cache,
                "ISMIP7_TIMING_CACHE_MANIFEST": manifest,
                "ISMIP7_TIMING_CACHE_AUDIT_OUTPUT": output_path,
                "ISMIP7_TIMING_CACHE_AUDIT_STATUS": status_path,
            }
            self._submit(
                f"timing_cache_audit_{lc}_{lc_coarse}",
                ncores,
                MEMORY_BY_LC[lc],
                exports,
                self.root
                / "scripts/batch_runners/timing_cache_audit.script",
                status_path,
            )

    def _lane_paths(self, lc, lc_coarse, ncores):
        return (
            self.timing_dir
            / timing_record_basename(self.lane_tag, lc, lc_coarse, ncores),
            self.timing_dir
            / timing_status_basename(self.lane_tag, lc, lc_coarse, ncores),
        )

    def lane_result(self, lc, lc_coarse, ncores):
        record_path, status_path = self._lane_paths(lc, lc_coarse, ncores)
        record, error = read_record(record_path)
        if error:
            return "invalid", error
        if record is not None:
            valid, detail = validate_timing_record(
                record,
                lc=lc,
                lc_coarse=lc_coarse,
                ncores=ncores,
                timing_tag=self.lane_tag,
            )
            return ("passed" if valid else "failed"), detail
        status, active = self.reconcile_status(status_path)
        if active:
            return "active", status["state"]
        if status and status.get("state") in {
            "failed", "finished", "submission_failed", "not_runnable"
        }:
            return "failed", status.get("category", status["state"])
        return "missing", "no record or active status"

    def probe_result(self, lc, lc_coarse, ncores):
        record_path = self.timing_dir / timing_record_basename(
            self.probe_tag, lc, lc_coarse, ncores
        )
        status_path = self.timing_dir / timing_status_basename(
            self.probe_tag, lc, lc_coarse, ncores
        )
        record, error = read_record(record_path)
        if error:
            return "invalid", error, status_path
        if record is not None:
            valid, detail = validate_timing_record(
                record,
                lc=lc,
                lc_coarse=lc_coarse,
                ncores=ncores,
                timing_kind="cache_probe",
                timing_tag=self.probe_tag,
            )
            return ("passed" if valid else "failed"), detail, status_path
        status, active = self.reconcile_status(status_path)
        if active:
            return "active", status["state"], status_path
        if status and status.get("state") in {
            "failed", "finished", "submission_failed", "not_runnable"
        }:
            return (
                "failed",
                status.get("category", status["state"]),
                status_path,
            )
        return "missing", "no record or active status", status_path

    def probe(self):
        if self.only_mesh is None:
            raise SystemExit(
                "The contract probe requires --only-mesh LC/LC_coarse"
            )
        lc, lc_coarse = self.only_mesh
        ncores = 32 if lc == 500 else 16
        state, detail, status_path = self.probe_result(
            lc, lc_coarse, ncores
        )
        if state in {"passed", "active"} and not self.force:
            print(f"PROBE {state.upper()} {lc}/{lc_coarse}: {detail}")
            return
        if state == "failed" and not self.force:
            print(
                f"PROBE FAILED {lc}/{lc_coarse}: {detail}; retry with "
                "FORCE_TIMING=1 after inspection"
            )
            return
        cache, manifest, cache_detail = self.lane_cache(lc, lc_coarse)
        if cache is None:
            print(f"PROBE WAITING CACHE {lc}/{lc_coarse}: {cache_detail}")
            return
        boundary = self.boundary_path(lc, lc_coarse)
        if not boundary.is_file() and not (
            self.dry_run and self.assume_valid_caches
        ):
            print(
                f"PROBE NOT RUNNABLE {lc}/{lc_coarse}: missing {boundary}"
            )
            if not self.dry_run:
                atomic_write_status(
                    status_path,
                    "not_runnable",
                    reason="missing_boundary_ids",
                )
            return
        dt = expected_dt(lc, self.matrix_dt_2500)
        exports = {
            "ISMIP7_LC": lc,
            "ISMIP7_LC_COARSE": lc_coarse,
            "ISMIP7_BUFFER_M": BUFFER_M,
            "ISMIP7_BNDIDS": boundary,
            "ISMIP7_INVERSION": self.inversion,
            "ISMIP7_RESTART": cache,
            "ISMIP7_TIMING_CACHE_MANIFEST": manifest,
            # validate_cache_manifest derives the accepted per-mesh MAP name
            # from this; export it rather than rely on --export=ALL.
            "ISMIP7_TIMING_INVERSION_MAXITER": self.maxiter,
            "ISMIP7_FRICTION": "budd",
            "ISMIP7_GEOMETRY_SPACE": "dg0",
            "ISMIP7_N_FLOW": "3.0",
            "ISMIP7_A4_FACTOR": "1.0",
            "ISMIP7_DIAGNOSTIC_LINEAR_SOLVER": SOLVER_MODE,
            "ISMIP7_SNES_DIVERGENCE_TOL": SNES_DIVERGENCE_TOL_DEFAULT,
            "ISMIP7_RESCUE_ENABLED": "0",
            "ISMIP7_SUBCYCLES": "1",
            "ISMIP7_T_END": (
                f"{expected_t_end(lc, self.matrix_steps, self.matrix_dt_2500):.12g}"
            ),
            "ISMIP7_DT": f"{dt:.12g}",
            "ISMIP7_TIMING_KIND": "cache_probe",
            "ISMIP7_TIMING_TAG": self.probe_tag,
            "ISMIP7_TIMING_EXPERIMENT": (
                f"timing_{self.probe_tag}_lcc{lc_coarse}_n{ncores}"
            ),
            "ISMIP7_TIMING_STATUS": status_path,
        }
        exports.update(self.lane_exports())
        if self.monitor:
            exports.update({
                "ISMIP7_SNES_MONITOR": "1",
                "ISMIP7_SNES_LOG": self.logs_dir / (
                    f"timing_snes_{self.probe_tag}_{lc}_{lc_coarse}"
                    f"_{ncores}.log"
                ),
            })
        self._submit(
            f"timing_probe_{self.contract}_{lc}_{lc_coarse}_{ncores}",
            ncores,
            MEMORY_BY_LC[lc],
            exports,
            self.root / "scripts/batch_runners/timing_transient.script",
            status_path,
        )

    def _submit_lane(self, lane):
        lc, lc_coarse, ncores = lane
        state, detail = self.lane_result(*lane)
        if state in {"passed", "active"} and not self.force:
            print(f"LANE {state.upper()} {lc}/{lc_coarse}/{ncores}: {detail}")
            return
        if state == "failed" and not self.force:
            print(f"LANE FAILED {lc}/{lc_coarse}/{ncores}: {detail}")
            return
        _, status_path = self._lane_paths(lc, lc_coarse, ncores)
        dependency = None
        if self.initial_state == "prepare":
            # Control lane: the transferred prepare state, no invert gate.
            inv_state, inv_detail, inv_status_path = "passed", "control", None
        else:
            inv_state, inv_detail, inv_status_path = self.inversion_result(
                lc, lc_coarse
            )
        if inv_state != "passed":
            if inv_state == "failed":
                print(
                    f"BLOCKED BY INVERSION {lc}/{lc_coarse}/{ncores}: "
                    f"{inv_detail}"
                )
                if not self.dry_run:
                    atomic_write_status(
                        status_path,
                        "blocked_by_inversion",
                        reason=inv_detail,
                    )
                return
            # Scout/scale need the per-mesh invert MAP + republished cache,
            # unless --follow-invert queues behind an active invert allocation.
            if not self.follow_invert:
                print(
                    f"WAITING FOR INVERSION {lc}/{lc_coarse}/{ncores}: "
                    f"{inv_state} ({inv_detail})"
                )
                return
            inv_status, inv_active = self.reconcile_status(inv_status_path)
            inv_job = inv_status.get("job_id") if inv_status else None
            if inv_active and inv_job:
                dependency = f"afterok:{inv_job}"
                print(
                    f"LANE FOLLOW INVERT {lc}/{lc_coarse}/{ncores}: "
                    f"dependency={dependency}"
                )
            elif (
                inv_status
                and inv_status.get("state") == "finished"
                and inv_job
            ):
                print(
                    f"WAITING FOR INVERSION {lc}/{lc_coarse}/{ncores}: "
                    f"{inv_detail} (invert stamped finished)"
                )
                return
            else:
                print(
                    f"WAITING FOR INVERSION {lc}/{lc_coarse}/{ncores}: "
                    f"{inv_state} ({inv_detail}); no active invert job to follow"
                )
                return
        # Campaign lanes (and re-inverted lanes on meshes that are not
        # re-inverted) run from the transferred prepare state; re-inverted
        # lanes on inverted meshes need that mesh's invert-published cache.
        from_prepare = (
            self.initial_state == "prepare" or not inversion_required(lc)
        )
        if from_prepare:
            cache, manifest, cache_detail = self.lane_cache(lc, lc_coarse)
            if cache is None:
                print(
                    f"LANE WAITING CACHE {lc}/{lc_coarse}/{ncores}: "
                    f"{cache_detail}"
                )
                return
        else:
            valid, cache_detail = self.cache_validation(
                lc, lc_coarse, require_mesh_inversion=True
            )
            if not valid:
                print(
                    f"LANE WAITING CACHE {lc}/{lc_coarse}/{ncores}: "
                    f"{cache_detail}"
                )
                return
            cache, manifest = cache_paths(self.cache_dir, lc, lc_coarse)
        boundary = self.boundary_path(lc, lc_coarse)
        if not boundary.is_file() and not (
            self.dry_run and self.assume_valid_caches
        ):
            print(
                f"LANE NOT RUNNABLE {lc}/{lc_coarse}/{ncores}: "
                f"missing {boundary}"
            )
            if not self.dry_run:
                atomic_write_status(
                    status_path,
                    "not_runnable",
                    reason="missing_boundary_ids",
                )
            return
        dt = expected_dt(lc, self.matrix_dt_2500)
        exports = {
            "ISMIP7_LC": lc,
            "ISMIP7_LC_COARSE": lc_coarse,
            "ISMIP7_BUFFER_M": BUFFER_M,
            "ISMIP7_BNDIDS": boundary,
            "ISMIP7_INVERSION": (
                self.inversion if from_prepare
                else self.inversion_map_path(lc, lc_coarse)
            ),
            "ISMIP7_RESTART": cache,
            "ISMIP7_TIMING_CACHE_MANIFEST": manifest,
            "ISMIP7_TIMING_INVERSION_MAXITER": self.maxiter,
            "ISMIP7_FRICTION": "budd",
            "ISMIP7_GEOMETRY_SPACE": "dg0",
            "ISMIP7_N_FLOW": "3.0",
            "ISMIP7_A4_FACTOR": "1.0",
            "ISMIP7_DIAGNOSTIC_LINEAR_SOLVER": SOLVER_MODE,
            "ISMIP7_SNES_DIVERGENCE_TOL": SNES_DIVERGENCE_TOL_DEFAULT,
            "ISMIP7_RESCUE_ENABLED": "0",
            "ISMIP7_SUBCYCLES": "1",
            "ISMIP7_T_END": (
                f"{expected_t_end(lc, self.matrix_steps, self.matrix_dt_2500):.12g}"
            ),
            "ISMIP7_DT": f"{dt:.12g}",
            "ISMIP7_TIMING_KIND": "matrix",
            "ISMIP7_TIMING_TAG": self.lane_tag,
            "ISMIP7_TIMING_EXPERIMENT": (
                f"timing_{self.lane_tag}_lcc{lc_coarse}_n{ncores}"
            ),
            "ISMIP7_TIMING_STATUS": status_path,
        }
        exports.update(self.lane_exports())
        if self.monitor:
            exports.update({
                "ISMIP7_SNES_MONITOR": "1",
                "ISMIP7_SNES_LOG": self.logs_dir / (
                    f"timing_snes_{self.lane_tag}_{lc}_{lc_coarse}"
                    f"_{ncores}.log"
                ),
            })
        self._submit(
            f"timing_{lc}_{lc_coarse}_{ncores}",
            ncores,
            MEMORY_BY_LC[lc],
            exports,
            self.root / "scripts/batch_runners/timing_transient.script",
            status_path,
            dependency=dependency,
        )

    def scout(self):
        for lane in self.selected_lanes(scout_lanes()):
            self._submit_lane(lane)

    def scale(self):
        scouts = {(lc, lc_coarse): ncores
                  for lc, lc_coarse, ncores in scout_lanes()}
        for lane in self.selected_lanes(scaling_lanes()):
            lc, lc_coarse, ncores = lane
            scout_ncores = scouts[(lc, lc_coarse)]
            scout_state, detail = self.lane_result(
                lc, lc_coarse, scout_ncores
            )
            _, status_path = self._lane_paths(*lane)
            decision = scaling_decision(scout_state)
            if decision == "block":
                print(
                    f"BLOCKED BY SCOUT {lc}/{lc_coarse}/{ncores}: {detail}"
                )
                if not self.dry_run:
                    atomic_write_status(
                        status_path,
                        "blocked_by_scout",
                        scout_cores=scout_ncores,
                        reason=detail,
                    )
                continue
            if decision == "wait":
                print(
                    f"WAITING FOR SCOUT {lc}/{lc_coarse}/{ncores}: "
                    f"{scout_state} ({detail})"
                )
                continue
            self._submit_lane(lane)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "stage",
        choices=("prepare", "audit", "probe", "invert", "scout", "scale"),
    )
    parser.add_argument("--root", default=_ROOT)
    parser.add_argument("--cache-dir", default=_ROOT / "results/timing/cache")
    parser.add_argument("--timing-dir", default=_ROOT / "results/timing")
    parser.add_argument(
        "--inversion",
        default=(
            _ROOT / "results/timing/inversion" / SOURCE_INVERSION_BASENAME
        ),
        help="campaign source MAP (the Makefile's TIMING_INVERSION)",
    )
    parser.add_argument(
        "--queue",
        choices=("short", "long", "debug"),
        default="short",
        help="the site's partition of this class (batch_runners/sites/)",
    )
    parser.add_argument(
        "--partition",
        default=None,
        help="a partition by name, instead of the site's --queue class",
    )
    parser.add_argument(
        "--constraint",
        default=None,
        help="a node feature by name, instead of the site's "
        "ISMIP7_CONSTRAINT_TIMING; also skips the site's per-node limits",
    )
    parser.add_argument("--walltime", default="12:00:00")
    parser.add_argument(
        "--maxiter",
        type=int,
        default=None,
        help=(
            "L-BFGS iterations for the timing-inversion stage "
            f"(default {inversion_maxiter()})"
        ),
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--assume-valid-caches",
        action="store_true",
        help="dry-run only: print downstream commands before caches exist",
    )
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--follow-prepare",
        action="store_true",
        help=(
            "invert only: submit behind each mesh's active prepare job with "
            "Slurm --dependency=afterok:<prepare_job_id>"
        ),
    )
    parser.add_argument(
        "--follow-invert",
        action="store_true",
        help=(
            "scout/scale: submit behind each mesh's active invert job with "
            "Slurm --dependency=afterok:<invert_job_id>"
        ),
    )
    parser.add_argument(
        "--only-mesh",
        type=_parse_mesh,
        help="operate only on one LC/LC_coarse pair, for example 2500/25000",
    )
    parser.add_argument(
        "--monitor",
        action="store_true",
        help="enable per-lane SNES/KSP logs under results/logs",
    )
    parser.add_argument(
        "--initial-state",
        choices=("prepare", "invert"),
        default=LANE_INITIAL_STATE_DEFAULT,
        help=(
            "scout/scale: 'prepare' (campaign lanes from the transferred "
            "prepare cache) or 'invert' (lanes from each mesh's per-mesh "
            "invert-published cache, recorded under the _reinverted tag)"
        ),
    )
    parser.add_argument(
        "--matrix-steps",
        type=int,
        default=None,
        help=f"steps per lane (default {matrix_steps()}; encoded in the tag)",
    )
    parser.add_argument(
        "--matrix-dt-2500",
        type=float,
        default=None,
        help=(
            f"timestep at 2.5 km in years, scaled by LC/2500 for other meshes "
            f"up to {MATRIX_DT_MAX:g} yr (default {matrix_dt_2500():g}; encoded "
            f"in the tag)"
        ),
    )
    parser.add_argument(
        "--contract",
        choices=tuple(CONTRACTS),
        default=None,
        help=(
            f"lane physics contract (default {contract_name()}): strict = no "
            "apparent MB and no calving sink; divfront/front/div add the "
            "production closure and/or the fixed 2015 calving front"
        ),
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.assume_valid_caches and not args.dry_run:
        raise SystemExit("--assume-valid-caches is allowed only with --dry-run")
    os.environ["ISMIP7_DIAGNOSTIC_LINEAR_SOLVER"] = SOLVER_MODE
    # A timing campaign is a fixed solver configuration, not an ambient-shell
    # experiment.  Pin PETSc's actual PETSC_UNLIMITED sentinel so an exported
    # legacy -1 cannot silently restore the 1e4 DIVERGED_DTOL cutoff.
    os.environ["ISMIP7_SNES_DIVERGENCE_TOL"] = SNES_DIVERGENCE_TOL_DEFAULT
    os.environ["ISMIP7_FRICTION"] = "budd"
    os.environ["ISMIP7_GEOMETRY_SPACE"] = "dg0"
    os.environ["ISMIP7_N_FLOW"] = "3.0"
    os.environ["ISMIP7_A4_FACTOR"] = "1.0"
    manager = CampaignManager(args)
    getattr(manager, args.stage)()
    if manager.submit_failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
