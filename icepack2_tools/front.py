r"""Front bookkeeping the thickness transport does around the level set.

The level set decides where the front IS (see :mod:`icepack2_tools.levelset`).
These are the rules the DG0 transport applies to the cells it touches once the
front has moved. They live here, as pure functions over the cell arrays, so
they can be exercised without standing up a whole run.
"""

import numpy as np

__all__ = ["retreat_slivers", "clear_reference_where_ice_free", "clamp_thickness"]


def clamp_thickness(h, h_clamp, *ice_free):
    r"""Floor the cell thicknesses to ``h_clamp``, except where there is no ice.

    ``h`` is the per-cell thickness array, mutated in place and returned; each
    ``ice_free`` argument is a boolean mask (or ``None``) naming cells that
    hold no ice. The caller measures the mass the floor added by integrating
    before and after.

    The exemption is the whole point. A cell outside the ice domain floored to
    ``h_clamp`` is handed that much fresh ice out of nothing, and whatever rule
    empties it takes it away again on the very next advance, so the run
    fabricates a steady ``h_clamp / dt`` of both clamp mass and calving flux
    there forever. Every rule that empties a cell therefore has to name it
    here:

    * the level set's ice-free extent, not only the cells calved this step:
      under a free law the cells the front just passed are a handful, so
      flooring the rest of the buffer would fabricate ice across every
      never-glaciated cell;
    * the cells beyond a pinned or fixed front;
    * the FLOATING cells the ISMIP7 collapse mask empties, which is the same
      case: the shelf is gone there, and a floored cell would be re-calved
      into ``licalvf`` on every advance for the rest of the run.
    """
    floor = np.full_like(h, h_clamp)
    for mask in ice_free:
        if mask is not None:
            floor[mask] = 0.0
    np.maximum(h, floor, out=h)
    return h


def clear_reference_where_ice_free(a_ref, ice_free):
    r"""Zero the frozen apparent-MB reference in the cells that hold no ice.

    ``a_ref`` is the per-cell reference array, mutated in place; ``ice_free``
    is the level set's ice-free mask for the CURRENT extent.

    The reference is the t=0 flux divergence, so at a t=0 front cell it is the
    terminus outflow: large and positive. Applied against the live extent it
    closes both directions of the same rule - no balancing reference where
    there is no ice:

    * advance - a frozen SINK outside the extent re-empties the cells a free
      front advances into;
    * retreat - a frozen SOURCE inside the t=0 extent regrows the cells a free
      front has just calved, so the front cannot retreat and its calving tally
      counts the regrown ice again every step.

    Zeroing in place is deliberate: a cell that later re-enters the ice keeps
    ``a_ref = 0``, because the frozen reference was only ever defined on the
    t=0 ice.
    """
    a_ref[ice_free] = 0.0
    return a_ref


def retreat_slivers(h_new, h_old, front_hmin):
    r"""Mask of cells holding a retreat remainder the front should take.

    ``h_new`` and ``h_old`` are cell thicknesses after and before one transport
    advance, over the same owned cells. A cell qualifies when it HELD ice
    (``h_old > front_hmin``) and ended the step below the extent threshold
    (``0 < h_new <= front_hmin``). That remainder is what the sub-cell calving
    shed and the melt left behind, so it is the front's own loss and belongs
    to the calving tally.

    A cell that was ice-free is never touched, however little it holds. That
    sub-threshold inflow is how a free front ADVANCES: it accumulates over
    successive steps until it crosses ``front_hmin`` and joins the extent.
    Removing it every step would pin the front wherever the one-step influx
    is under the threshold.
    """
    return (h_new > 0.0) & (h_new <= front_hmin) & (h_old > front_hmin)
