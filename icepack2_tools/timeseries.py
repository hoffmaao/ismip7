r"""The rows a forward run appends to its timeseries CSV.

The year column carries four decimals so a sub-0.1 yr step stays resolvable:
readers such as ``check_ismip6_track.infer_dt`` recover the timestep from it,
and every rate they report is a per-step quantity divided by that timestep.
"""

from icepack2_tools.front import collapse_csv_fields

__all__ = ["timeseries_csv_line"]


def timeseries_csv_line(row, header, collapse_cells):
    r"""One timeseries line for ``row``, newline included.

    ``row`` is ``(year, vaf_mm_sle, mass_gt, *budget_columns)``; ``header`` and
    ``collapse_cells`` go to :func:`icepack2_tools.front.collapse_csv_fields`.
    """
    return (
        f"{row[0]:.4f},{row[1]:.6f},{row[2]:.2f},"
        + ",".join(f"{v:.4f}" for v in row[3:])
        + collapse_csv_fields(header, collapse_cells) + "\n"
    )
