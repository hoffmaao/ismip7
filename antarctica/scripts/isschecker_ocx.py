#!/usr/bin/env python
r"""isschecker 0.5.1 with core 11 checked (issue #18).

    <venv>/bin/python antarctica/scripts/isschecker_ocx.py \
        --source-path submission/AIS/RICE/icepack2/CORE/C011 --variable-list ismip7

The arguments, the log and the exit status are the checker's own. Run it with
the interpreter of the tools venv (``antarctica/README.md``, "From a finished
run to a submission"); it imports nothing from this repository.

Stock 0.5.1 cannot check core 11. Its experiment table has no ``ocx`` row, so
a C011 set draws one naming error and none of its files is read, and field 5
of a filename takes CMIP models only, where core 11 names the reanalysis,
``ERA5``, as the organisers' conventions document does. This adds, in this
process only:

- the ``ocx`` experiment, the Protocol Overview's C011 row: a start from 1990
  to 2015 and an end in 2025, the start being free inside that window as the
  historical's is (duration -1);
- ``ERA5`` in field 5, for files whose experiment is ``ocx``.

Every other check is the release's own, and a set of any other core is
checked as stock 0.5.1 checks it. The log's version line names the patch. It
reaches into the checker's private functions, so it runs only on the release
it was written against, and stops once a release carries an ``ocx`` row.
"""
import sys

PATCHED_RELEASE = "0.5.1"
OCX_EXPERIMENT = {"experiment": "ocx", "start_year_min": 1990, "start_year_max": 2015,
                  "end_year": 2025, "duration": -1}
# write_ismip7_output.OCX_FORCING, which the unit suite holds this equal to
OCX_FORCING = "ERA5"
PATCH_NOTE = (f"with the {OCX_EXPERIMENT['experiment']} experiment and {OCX_FORCING} "
              f"in field 5 of its files added (icepack/ismip7 issue #18)")


def patch(checker):
    r"""Add the ``ocx`` experiment and its forcing name to ``checker``, the
    ``isschecker.checker`` module, and return the version line for the log."""
    if checker.__version__ != PATCHED_RELEASE:
        raise SystemExit(f"isschecker {checker.__version__}: this patch was written against "
                         f"{PATCHED_RELEASE}; check that a release still needs it (issue #18)")
    load_experiments = checker._load_experiments_csv
    check_naming = checker._check_naming

    def load_with_ocx(*args, **kwargs):
        experiments = load_experiments(*args, **kwargs)
        if any(e["experiment"] == OCX_EXPERIMENT["experiment"] for e in experiments):
            raise SystemExit("isschecker has an ocx experiment of its own; retire "
                             "isschecker_ocx.py (issue #18)")
        return experiments + [dict(OCX_EXPERIMENT)]

    def check_naming_with_ocx_forcing(reporter, file_name, *args, **kwargs):
        parts = file_name.split("_")
        if (len(parts) != checker.ISMIP7_FILENAME_PARTS
                or parts[checker.ISMIP7_FILENAME_EXPERIMENT_IDX] != OCX_EXPERIMENT["experiment"]
                or OCX_FORCING in checker.VALID_ESM_NAMES):
            return check_naming(reporter, file_name, *args, **kwargs)
        checker.VALID_ESM_NAMES.add(OCX_FORCING)
        try:
            return check_naming(reporter, file_name, *args, **kwargs)
        finally:
            checker.VALID_ESM_NAMES.discard(OCX_FORCING)

    checker._load_experiments_csv = load_with_ocx
    checker._check_naming = check_naming_with_ocx_forcing
    return f"{checker._describe_version()} {PATCH_NOTE}"


def main():
    import isschecker.checker as checker

    version = patch(checker)
    args = checker._parse_args()
    print(f"isschecker {version}", flush=True)
    summary = checker.run_checker(source_path=args.source_path, variable_list=args.variable_list,
                                  output_path=args.output_path, version=version)
    return 1 if summary["fatal"] or summary["total_errors"] > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
