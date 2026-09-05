r"""ISMIP7's calving front: the shared level set, anchored at the ice extent.

The implementation lives in :mod:`icepack_tools.levelset` and is shared with
the CalvingMIP project, so a law tuned in one runs unchanged in the other and
the front is described by one object everywhere:

* ``phi`` is a DG0 signed distance on the transport's own cells, negative in
  ice;
* it is the solution of the eikonal problem ``|grad phi| = 1`` with
  ``phi = 0`` on the ice/water facets of the current extent -- the boundary
  condition sits INSIDE the mesh, at the ice sheet's own edge, so a buffered
  mesh needs no condition on its outer boundary;
* the calving rate retreats it by normal flow (Hahn, Mikula & Frolkovic 2025,
  arXiv:2504.05845), and :meth:`~icepack_tools.levelset.LevelSet.calving_masks`
  reports both the cells the front has passed and the sub-cell mass the front
  cells shed;
* advance is the thickness transport's job.

Everything ISMIP7 adds is here: the ``extent`` anchor as the default (the
sheet's thickness is the authority on where ice is) and the knob names, which
come from :mod:`icepack2_tools.runconfig` (``ISMIP7_CALVING``,
``ISMIP7_CALVING_SIGMA_MAX_GROUNDED`` / ``_FLOATING``).  Their tests are
``icepack_tools/test/levelset_test.py``, which cover both anchors serially
and on three ranks.
"""

from icepack_tools.levelset import (
    ANCHORS, LAWS, LevelSet as _LevelSet,
)

__all__ = ["ANCHORS", "LAWS", "LevelSet"]


class LevelSet(_LevelSet):
    r"""The shared level set with ISMIP7's anchoring default.

    Identical to :class:`icepack_tools.levelset.LevelSet` except that
    ``anchor`` defaults to ``"extent"``: on an ice sheet with melt and a
    buffered ocean the thickness decides where ice is, and re-solving the
    eikonal problem there each step keeps the front and the mass conservation
    describing the same ice.
    """

    def __init__(self, *args, anchor="extent", **kwargs):
        super().__init__(*args, anchor=anchor, **kwargs)
