r"""The per-basin melt observation table: two published widths, one reader.

The ISMIP7 melt-calibration product ships ``Melt_Paolo_Davison_Adusumilli
_imbie2.csv`` with three columns (basin, melt, uncertainty). The table it
replaces, ``Melt_Paolo_Err_Adusumilli_imbie2_v3.csv``, carries area and two
per-area columns between the melt and its uncertainty, so a reader keyed on
column position silently reads AREA as the uncertainty of every basin: values
near 1e5 km^2 where the uncertainty is tens of Gt/yr, which reweights the
whole calibration without failing. Columns are located by header name instead,
and this pins that.

The integrated targets differ by 23% (865.0 against 1067.4 Gt/yr), so which
table a calibration used is part of its provenance and belongs in the saved
npz.

Serial, no data files, no Firedrake.
"""

import csv
import os
import sys

import numpy as np
import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_SCRIPTS = os.path.join(_ROOT, "antarctica", "scripts")

OLD_HEADER = ("", "BMR (Gt/yr)", "Area (km^2)", "BMR uncert (Gt/yr)",
              "Average BMR (kg/m2/a)", "Average BMR uncert (kg/m2/a)")
NEW_HEADER = ("", "BMR (Gt/yr)", "BMR uncert (Gt/yr)")

# basin, melt, uncertainty; the area column exists only in the old layout and
# is deliberately far from the uncertainty so a positional read stands out.
ROWS = ((0, 39.85, 51.75, 124527.1), (1, 3.30, 3.50, 5455.0))


def _write(path, header, old_layout):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for bid, melt, sigma, area in ROWS:
            if old_layout:
                w.writerow([bid, melt, area, sigma, melt / area * 1e9,
                            sigma / area * 1e9])
            else:
                w.writerow([bid, melt, sigma])


def _load_obs_with(monkeypatch, path):
    r"""Import ``calibrate_melt`` against ``path`` and return its parse.

    The module resolves OBS_CSV at import, so the environment has to carry the
    path before the import and the module has to be dropped afterwards.
    """
    monkeypatch.setenv("ISMIP7_MELT_OBS_CSV", str(path))
    monkeypatch.syspath_prepend(_SCRIPTS)
    sys.modules.pop("calibrate_melt", None)
    try:
        import calibrate_melt
        return calibrate_melt._load_obs()
    finally:
        sys.modules.pop("calibrate_melt", None)


@pytest.mark.parametrize("old_layout", (True, False))
def test_both_published_layouts_read_the_same(tmp_path, monkeypatch,
                                              old_layout):
    path = tmp_path / ("old.csv" if old_layout else "new.csv")
    _write(path, OLD_HEADER if old_layout else NEW_HEADER, old_layout)

    bids, melt, sigma = _load_obs_with(monkeypatch, path)

    assert list(bids) == [row[0] for row in ROWS]
    assert np.allclose(melt, [row[1] for row in ROWS])
    # The uncertainty, never the area column that sits where it used to.
    assert np.allclose(sigma, [row[2] for row in ROWS])
    assert sigma.max() < 1e3


def test_a_table_without_the_named_columns_is_refused(tmp_path, monkeypatch):
    path = tmp_path / "wrong.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["", "melt", "error"])
        w.writerow([0, 39.85, 51.75])

    with pytest.raises(ValueError, match="bmr"):
        _load_obs_with(monkeypatch, path)
