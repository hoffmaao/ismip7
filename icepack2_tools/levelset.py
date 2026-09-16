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
sheet's thickness is the authority on where ice is), :func:`initial_distance`
for the t=0 anchor the ``fixed`` law freezes on, and the knob names, which
come from :mod:`icepack2_tools.runconfig` (``ISMIP7_CALVING``,
``ISMIP7_CALVING_SIGMA_MAX_GROUNDED`` / ``_FLOATING``).  The shared class's
tests are ``icepack_tools/test/levelset_test.py``, which cover both anchors
serially and on three ranks; what this module adds is tested in
``tests/test_initial_distance.py``.
"""


from icepack_tools.levelset import LevelSet as _LevelSet

__all__ = ["LevelSet", "initial_distance"]


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


class _DistanceOnly(LevelSet):
    r"""Construction stops at the eikonal solve: no advection or extension
    solvers, because a caller after the t=0 distance never advances this
    object."""

    def _build_solvers(self):
        pass


def initial_distance(mesh, h_dg, h_min=1.0):
    r"""The signed distance to the extent of ``h_dg``, as a DG0 Function.

    Same field as ``LevelSet(mesh, h_dg, law="none", h_min=h_min).phi``, but
    without building the advection and extension solvers, which a caller that
    only wants the t=0 distance never uses.

    The shared class also prints its front banner unconditionally in
    ``__init__``, so a `fixed` run logs one naming ``law=none`` above the one
    naming the law actually in force.  A distance-only entry point belongs in
    the toolbox itself (plan Phase 0); until then the extra line stands.
    """
    return _DistanceOnly(
        mesh, h_dg, law="none", h_min=h_min, drag_mask=None,
    ).phi
