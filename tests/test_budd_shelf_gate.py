r"""The Budd shelf gate: floating ice carries no Weertman friction.

``N = max(p_I - p_W, 0)`` is built from the MODEL surface, so on a shelf the
two pressures cancel to a roundoff residue whose sign is meaningless.  The
old gate was ``conditional(gt(N, 0), ...)``, a sign test on that residue: a
positive residue passed, and with ``N_ref`` equally tiny the delta floor
lifted it straight to ``nhat_cap`` - triple Weertman friction on floating
ice.  Gating on ``He`` alone is not enough either; ice floating by a few
metres sits inside the grounding-zone band where ``He`` is O(0.1).  The
production gate is height above flotation.

Serial, one 6x2 rectangle mesh, no data files.
"""

# dual_friction first: it pulls icepack2 -> irksome, which must be imported
# before any UFL form is assembled.
from icepack2_tools.dual_friction import (                      # noqa: E402
    budd_nhat,
    budd_nhat_ungated,
    effective_pressure,
    grounded_mask,
)

import numpy as np                                              # noqa: E402
import pytest                                                   # noqa: E402

from firedrake import (                                         # noqa: E402
    Function,
    FunctionSpace,
    RectangleMesh,
    SpatialCoordinate,
    conditional,
)
from icepack2.constants import (                                # noqa: E402
    ice_density as RHO_I,
    water_density as RHO_W,
)

NHAT_FLOOR = 0.02
NHAT_CAP = 3.0

# The far shelf: 500 m of ice over a 2 km trough, floating by ~1.7 km.
H_FAR, B_FAR = 500.0, -2000.0
# The grounding zone: floating by 3 m, well inside the 10 m He band.
B_NEAR, HAF_NEAR = -100.0, -3.0
# Grounded ice on the same bed, 50 m above flotation.
B_GND, HAF_GND = -100.0, 50.0
# A surface offset far above the cancellation roundoff (~1e-15 MPa) yet far
# below the 1e-6 MPa N_ref floor: it makes the shelf residue reliably
# POSITIVE, which is exactly what the old sign test waved through.
S_OFFSET = 1e-9


def _haf_to_thickness(b, haf):
    return haf + (RHO_W / RHO_I) * max(-b, 0.0)


@pytest.fixture(scope="module")
def three_regions():
    r"""Three columns of cells: far shelf (x < 1), near-GL shelf (1 < x < 2),
    grounded (x > 2).  Both shelf columns carry the hydrostatic surface, so
    their ``N`` is pure roundoff; the grounded column has ``s = b + H``."""
    mesh = RectangleMesh(6, 2, 3.0, 1.0)
    Q0 = FunctionSpace(mesh, "DG", 0)
    x, _ = SpatialCoordinate(mesh)

    h_near = _haf_to_thickness(B_NEAR, HAF_NEAR)
    h_gnd = _haf_to_thickness(B_GND, HAF_GND)

    H = Function(Q0, name="thickness").interpolate(
        conditional(x < 1.0, H_FAR, conditional(x < 2.0, h_near, h_gnd))
    )
    b = Function(Q0, name="bed").interpolate(
        conditional(x < 1.0, B_FAR, conditional(x < 2.0, B_NEAR, B_GND))
    )
    s_float = (1.0 - RHO_I / RHO_W) * H + S_OFFSET
    s = Function(Q0, name="surface").interpolate(
        conditional(x < 2.0, s_float, b + H)
    )
    centroid_x = Function(Q0).interpolate(x).dat.data_ro
    regions = {
        "far": centroid_x < 1.0,
        "near": (centroid_x > 1.0) & (centroid_x < 2.0),
        "grounded": centroid_x > 2.0,
    }
    return Q0, H, b, s, regions


def _dg0(Q0, expr):
    return Function(Q0).interpolate(expr).dat.data_ro


def test_shelf_residue_is_a_positive_roundoff_pressure(three_regions):
    r"""The premise of the whole gate: on the flotation branch ``N`` is not
    zero, it is a positive number ~11 orders below any real N."""
    Q0, H, b, s, regions = three_regions
    N = _dg0(Q0, effective_pressure(H, s))
    shelf = regions["far"] | regions["near"]
    assert np.all(N[shelf] > 0.0)
    assert np.all(N[shelf] < 1e-6)
    assert np.all(N[regions["grounded"]] > 0.1)


def test_ungated_nhat_reaches_the_cap_on_the_shelf(three_regions):
    r"""What the old sign test let through: the delta floor divides the local
    overburden by a roundoff ``N_ref``, so the shelf lands on ``nhat_cap``."""
    Q0, H, b, s, regions = three_regions
    N = effective_pressure(H, s)
    nhat_ungated = _dg0(
        Q0, budd_nhat_ungated(N, None, H, nhat_floor=NHAT_FLOOR, nhat_cap=NHAT_CAP)
    )
    assert np.allclose(nhat_ungated[regions["far"]], NHAT_CAP)


def test_production_gate_zeroes_both_shelf_regions(three_regions):
    r"""The fix: exactly zero on floating ice, grounded friction untouched."""
    Q0, H, b, s, regions = three_regions
    N = effective_pressure(H, s)
    nhat = _dg0(
        Q0, budd_nhat(N, None, H, b, nhat_floor=NHAT_FLOOR, nhat_cap=NHAT_CAP)
    )
    assert np.all(nhat[regions["far"]] == 0.0)
    assert np.all(nhat[regions["near"]] == 0.0)
    assert np.all(nhat[regions["grounded"]] > 0.0)


def test_the_gate_does_not_scale_grounded_friction(three_regions):
    r"""The gate is exact, not a smooth ramp. Grounded ice gets the ungated
    N_hat unchanged: an interim form multiplied through by ``He``, which cut
    friction everywhere inside the He band, and the Weertman branch already
    carries its own ``exp(theta * He)`` gate on theta."""
    Q0, H, b, s, regions = three_regions
    N = effective_pressure(H, s)
    ungated = _dg0(
        Q0, budd_nhat_ungated(N, None, H, nhat_floor=NHAT_FLOOR, nhat_cap=NHAT_CAP)
    )
    nhat = _dg0(
        Q0, budd_nhat(N, None, H, b, nhat_floor=NHAT_FLOOR, nhat_cap=NHAT_CAP)
    )
    g = regions["grounded"]
    assert np.allclose(nhat[g], ungated[g], rtol=0.0, atol=0.0)
    # He is well below 1 there (50 m HAF over a 10 m band is not the issue;
    # the band is what the interim form scaled), so this is a real distinction.
    assert np.all(_dg0(Q0, grounded_mask(H, b))[g] < 1.0)


def test_he_alone_would_not_have_closed_the_near_gl_band(three_regions):
    r"""Why ``He`` is not the gate either: 3 m of flotation leaves ``He``
    around 0.35, so an He-only form still puts a third of the capped Weertman
    friction on floating ice."""
    Q0, H, b, s, regions = three_regions
    N = effective_pressure(H, s)
    He = grounded_mask(H, b)
    he_only = _dg0(
        Q0,
        He * budd_nhat_ungated(N, None, H, nhat_floor=NHAT_FLOOR, nhat_cap=NHAT_CAP),
    )
    assert np.all(he_only[regions["near"]] > 0.5)
