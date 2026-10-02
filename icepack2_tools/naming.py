r"""MAP and mesh filename construction shared by the run tools.

A MAP is only valid for the (friction law, flow exponent, geometry space) it
was inverted under - the inversion absorbs whatever those get wrong into
theta/phi - so each of them is in the filename. The name is built in one place
because every copy of the rule is a chance for a gate to look for a file the
run will never write (or to bless a MAP the run will never load).

The standard mesh basename and its parser share one pattern so checkpoint
provenance can be checked against the mesh builder's convention. Pure Python:
importable without Firedrake so the preflight stays fast.
"""

import os
import re

from icepack2_tools.runconfig import geometry_space, n_flow

_FRICTION_TAGS = {"regularized_coulomb": "_rc", "budd": "_budd"}
_MESH_BASENAME_RE = re.compile(
    r"antarctica_(\d+)_(\d+)(?:_buffered(\d+))?$"
)


def mesh_basename(lc_coarse, lc, buffer_m):
    r"""Standard Antarctica mesh basename without its extension."""
    return f"antarctica_{lc_coarse}_{lc}_buffered{int(float(buffer_m))}"


def parse_mesh_basename(name):
    r"""Return ``(lc_coarse, lc, buffer_m)`` from a standard mesh name.

    Directory and ``.msh`` extension components are accepted. Legacy names
    without a buffer tag return ``None`` for ``buffer_m``.
    """
    basename = os.path.basename(str(name))
    stem = basename[:-len(".msh")] if basename.endswith(".msh") else basename
    match = _MESH_BASENAME_RE.fullmatch(stem)
    if match is None:
        raise ValueError(
            f"Mesh filename {name!r} does not match "
            "antarctica_<COARSE>_<FINE>[_buffered<BUFFER_M>]"
        )
    lc_coarse, lc = int(match.group(1)), int(match.group(2))
    buffer_m = None if match.group(3) is None else int(match.group(3))
    return lc_coarse, lc, buffer_m


def map_n_tag():
    r"""Filename tag distinguishing MAPs inverted at different flow exponents
    so n=3 and n=4 MAPs coexist on disk. n=4 keeps the legacy untagged name
    (backward compatible with the `antarctica` MAPs); any other n gets
    `_n<N>` (e.g. `_n3`)."""
    n = n_flow()
    return "" if abs(n - 4.0) < 1e-9 else f"_n{int(round(n))}"


def map_geom_tag():
    r"""Filename tag for the geometry discretization a MAP was inverted under.

    CG1 keeps the legacy untagged name (every MAP on disk before Aug 2026);
    DG0 gets `_dg0`. They are NOT interchangeable: the inversion absorbs the
    front treatment into theta/phi, so a CG1 MAP driven by a DG0 forward has
    friction tuned to a calving front that the lumped lift made ~2x too thick.
    Separate names stop that from happening by accident.
    """
    return "" if geometry_space() == "cg1" else "_dg0"


def map_tag(friction, geometry=True):
    r"""Full MAP tag: friction law + flow exponent + geometry space."""
    return (_FRICTION_TAGS.get(friction, "") + map_n_tag()
            + (map_geom_tag() if geometry else ""))


def map_basename(friction, lc, geometry=True):
    r"""MAP filename (no directory) for this configuration."""
    return f"inversion_icepack2{map_tag(friction, geometry)}_{int(lc)}.h5"
