r"""The t=0 anchor for the `fixed` calving law.

``icepack2_tools.levelset.initial_distance`` exists so that a run with
``ISMIP7_CALVING=fixed`` gets the level set's t=0 distance field without the
side effects of a full construction: no advection/extension solvers, and no
second front banner claiming ``law=none`` above the real one.

Serial, one 8x8 unit mesh, no data files.
"""

import numpy as np
import pytest

from firedrake import Function, FunctionSpace, SpatialCoordinate, UnitSquareMesh

from icepack2_tools.levelset import LevelSet, initial_distance

HMIN = 1.0


@pytest.fixture(scope="module")
def thickness():
    r"""A left-half ice sheet on a unit mesh: h = 100 m for x < 0.5, ice-free
    beyond, so the extent has an interior front to anchor on."""
    mesh = UnitSquareMesh(8, 8)
    Q0 = FunctionSpace(mesh, "DG", 0)
    x, _ = SpatialCoordinate(mesh)
    h = Function(Q0).interpolate(100.0 * (x < 0.5))
    return mesh, h


def test_matches_the_full_construction(thickness):
    r"""The anchor is the same field the discarded full LevelSet produced."""
    mesh, h = thickness
    reference = LevelSet(mesh, h, law="none", h_min=HMIN, drag_mask=None).phi
    phi = initial_distance(mesh, h, h_min=HMIN)
    assert np.allclose(phi.dat.data_ro, reference.dat.data_ro)


def test_signs_follow_the_extent(thickness):
    r"""Negative in ice, positive beyond it: this is what `fixed` freezes."""
    mesh, h = thickness
    phi = initial_distance(mesh, h, h_min=HMIN)
    ice = h.dat.data_ro > HMIN
    assert np.all(phi.dat.data_ro[ice] < 0.0)
    assert np.all(phi.dat.data_ro[~ice] > 0.0)


def test_prints_no_front_banner(thickness, capfd):
    r"""No banner, so a `fixed` run logs exactly one level-set line: the one
    naming the law actually in force."""
    mesh, h = thickness
    capfd.readouterr()
    initial_distance(mesh, h, h_min=HMIN)
    assert "Level-set front" not in capfd.readouterr().out

    LevelSet(mesh, h, law="none", h_min=HMIN, drag_mask=None)
    assert "Level-set front" in capfd.readouterr().out


def test_banner_suppression_is_scoped(thickness):
    r"""The shared module's PETSc is restored, so later printing still works."""
    import icepack_tools.levelset as shared

    mesh, h = thickness
    before = shared.PETSc
    initial_distance(mesh, h, h_min=HMIN)
    assert shared.PETSc is before
