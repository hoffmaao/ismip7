r"""The horizontal-force-balance calving criterion, as UFL on the front cells.

A shelf front holds when the resistive stress it must carry is below what a
crevasse field can bear, and calves when it is not. Buck (2023) derived that
balance for an unbuttressed shelf; Coffey et al. (2024) and Coffey and Lai
(2025) wrote it as a buttressing number, where ``B = 0`` is the front in
equilibrium and no unbuttressed shelf should exist; Slater and Wagner (2025)
added the ice's tensile strength and basal friction. This is that criterion,
with the zero-stress (Nye) threshold beside it for comparison, which on a
shelf is twice as large.

Two quantities, both depth-averaged and in MPa:

``R_xx``, the resistive stress the front carries, is the front-normal
component of the membrane stress, ``n . M n``. The papers write it as
``R_xx = 2 tau_xx + tau_yy`` with ``x`` normal to the front, and icepack2's
membrane stress is ``M = tau_h + tr(tau_h) I``, so in front-aligned
coordinates ``M_nn`` is exactly their ``R_xx``. Taking the largest principal
value of ``M`` instead is never smaller and moves a front that should hold
into retreat, so the normal component is the only one offered here.

``R_crit``, the stress at which the crevasse field penetrates, is Slater and
Wagner's equation 21,

    R_crit / (rho_i g H) = (1 - (rho_i/rho_c)(1 - a)^2) / 2
                           + rho_c sigma~^2 / (2 (rho_c - rho_i))

with ``a`` the height-above-flotation fraction, ``rho_c`` the density of the
water filling a basal crevasse and ``sigma~`` the tensile strength over the
overburden. At zero strength and seawater in the crevasse this is Coffey and
Lai's ``B = 0`` and, afloat, Buck's result.

The rate is ``c = |u| min(R_xx / R_crit, ratio_max)^p``: the front holds where
the stress is at the threshold, retreats where crevasses would penetrate, and
advances where they would not. ``R_crit`` scales with the thickness, so on ice
thinning to nothing the ratio is unbounded and one step could carry the front
across many cells; ``ratio_max`` bounds what a step may remove and never binds
on thick ice.

Densities are icepack2's own, 917 and 1024, which are what the forward's
flotation test and driving stress use. The CalvingMIP protocol and the
`calving` project that tunes these laws use 1028, a 0.4 percent difference
that moves height above flotation by about 0.3 m on a 300 m shelf, far below
the scatter a tuned threshold carries.

Units: MPa, m, yr.
"""

from firedrake import Constant, conditional, dot, gt, max_value, min_value, sqrt
from icepack2.constants import gravity as g, ice_density as rho_I, water_density as rho_W

__all__ = ["resistive_stress", "critical_stress", "height_above_flotation",
           "buttressing_number", "hfb_calving_rate", "thickness_calving_rate",
           "density_in_model_units", "HFB_MODES"]

#: Seconds in icepack's year, which is how its densities carry their units.
_YEAR = 365.25 * 24 * 60 * 60


def density_in_model_units(rho_si):
    r"""A density in kg/m3, which is how the papers quote it, in icepack2's
    own MPa, m, yr. ``density_in_model_units(1024)`` is its water density."""
    return float(rho_si) / _YEAR ** 2 * 1e-6

#: ``hfb`` is the horizontal-force-balance threshold of the papers above;
#: ``zero_stress`` is the Nye criterion, which on a shelf is twice as large,
#: kept because Bassis et al. (2026) contest the force-balance depths for
#: closely spaced crevasses and the two can then be put to the same run.
HFB_MODES = ("hfb", "zero_stress")

#: A thickness below which the overburden is not taken literally, so a cell
#: thinning to nothing cannot divide by zero.
_H_FLOOR = 1.0


def height_above_flotation(h, b):
    r"""``h - (rho_w / rho_i) max(-b, 0)``, cell-wise, as the forward's own
    grounding test computes it."""
    return h - (rho_W / rho_I) * max_value(-b, Constant(0.0))


def resistive_stress(M, n):
    r"""The front-normal resistive stress ``n . M n`` [MPa], the papers'
    ``R_xx``. Negative under compression, so a caller wanting a stress the
    ice must carry floors it at zero."""
    return dot(n, dot(M, n))


def critical_stress(h, b, sigma_max=0.0, rho_c=None, mode="hfb"):
    r"""``R_crit`` [MPa]: the resistive stress at which the crevasse field
    penetrates the full thickness.

    ``sigma_max`` is the ice's tensile strength [MPa]; zero is the Buck,
    Coffey and Lai limit, in which an unbuttressed shelf sits exactly at the
    threshold. ``rho_c`` is the density of the water in a basal crevasse
    [kg/m3 in icepack2 units]; seawater by default, fresh water for a
    meltwater-filled crevasse. ``mode`` selects the force balance or the
    zero-stress threshold.
    """
    if mode not in HFB_MODES:
        raise ValueError(f"mode must be one of {HFB_MODES}, got {mode!r}")
    rho_c = rho_W if rho_c is None else rho_c
    thickness = max_value(h, Constant(_H_FLOOR))
    overburden = rho_I * g * thickness
    fraction = max_value(height_above_flotation(h, b), Constant(0.0)) / thickness
    ratio = Constant(float(rho_I / rho_c))
    if mode == "zero_stress":
        return overburden * (Constant(1.0) - ratio * (Constant(1.0) - fraction))
    strength = (Constant(float(rho_c / (rho_c - rho_I)))
                * (Constant(float(sigma_max)) / overburden) ** 2 / 2)
    return overburden * (
        (Constant(1.0) - ratio * (Constant(1.0) - fraction) ** 2) / 2 + strength)


def buttressing_number(M, n, h, b, sigma_max=0.0, rho_c=None, mode="hfb"):
    r"""Coffey and Lai's ``B = 1 - R_xx / R_crit``, floored at zero: the
    margin a front has left. Zero is a front at the threshold, and the
    theory says an unbuttressed shelf should sit there."""
    r_xx = max_value(resistive_stress(M, n), Constant(0.0))
    return max_value(Constant(1.0) - r_xx / critical_stress(
        h, b, sigma_max, rho_c, mode), Constant(0.0))


def hfb_calving_rate(u, M, h, b, n, sigma_max=0.0, rho_c=None, mode="hfb",
                     exponent=1.0, ratio_max=5.0, grounded_gate=False):
    r"""``c = |u| min(R_xx / R_crit, ratio_max)^p`` [m/yr along ``n``].

    ``grounded_gate`` zeroes the rate on grounded ice, for a run whose
    protocol says grounded ice does not calve.
    """
    r_xx = max_value(resistive_stress(M, n), Constant(0.0))
    r_crit = critical_stress(h, b, sigma_max, rho_c, mode)
    ratio = min_value(r_xx / r_crit, Constant(float(ratio_max)))
    speed = sqrt(dot(u, u) + Constant(1e-30))
    rate = speed * ratio ** Constant(float(exponent))
    if grounded_gate:
        haf = height_above_flotation(h, b)
        rate = rate * conditional(gt(haf, Constant(0.0)),
                                  Constant(0.0), Constant(1.0))
    return rate


def thickness_calving_rate(u, h, h_critical, grounded_gate=False, b=None):
    r"""``c = max(0, 1 + (Hc - H) / Hc) |u|`` [m/yr along the front normal].

    The minimum-thickness rule, in the rate form CalvingMIP's experiment 5
    prescribes (there with ``Hc = 375 m``). ``Hc`` is
    the thickness the front settles at, not a cutoff: at ``H = Hc`` the rate is
    exactly the speed the ice arrives with, so the front is stationary; below
    it the front retreats and above it the front advances. That makes the rule
    self-limiting, which a stress threshold is not.

    It is also the one rule here that reads no inferred field. Over most of an
    Antarctic front the observed speed is a few metres a year against a 3 m/yr
    error floor, so the inversion has no leverage there and the rheology a
    stress law would read is the regularizer's extrapolation rather than
    anything the data constrained; cells the front later advances into were
    never in the inversion's domain at all. A geometric rule is unaffected by
    both.

    On ``Hc``, note what does and does not transfer. Wilner et al. (2023)
    calibrate a minimum-thickness threshold against ten Antarctic shelves and
    get 55 to 440 m, but theirs is a position law (calve where ``h <= hmin``)
    and they find the values "largely dependent on the original thickness of
    the ice shelf" - unlike the von Mises strength in the same study, they are
    shelf-specific. The transferable quantity for a rate form is the observed
    front thickness, because that is where this rate holds a front still:
    BedMachine v4.1 gives a 144.7 m mean over 77762 Antarctic front cells and
    a 152.4 m floating-front median. ``b`` is only needed for
    ``grounded_gate``.
    """
    from firedrake import Constant, conditional, dot, gt, max_value, sqrt

    thickness = max_value(h, Constant(0.0))
    hc = Constant(float(h_critical))
    factor = max_value(Constant(1.0) + (hc - thickness) / hc, Constant(0.0))
    speed = sqrt(dot(u, u) + Constant(1e-30))
    rate = factor * speed
    if grounded_gate:
        if b is None:
            raise ValueError("grounded_gate needs the bed")
        rate = rate * conditional(
            gt(height_above_flotation(h, b), Constant(0.0)),
            Constant(0.0), Constant(1.0))
    return rate
