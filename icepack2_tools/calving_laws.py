r"""The horizontal-force-balance calving criterion, as UFL on the front cells.

A shelf front holds when the resistive stress it must carry is below what a
crevasse field can bear, and calves when it is not. Buck (2023) derived that
balance for an unbuttressed shelf; Coffey et al. (2024) and Coffey and Lai
(2025) wrote it as a buttressing number, where ``B = 0`` is the front in
equilibrium and no unbuttressed shelf should exist; Slater and Wagner (2025)
added the ice's tensile strength and basal friction. This is that criterion.

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
flotation test and driving stress use.

Units: MPa, m, yr.
"""

from firedrake import Constant, dot, max_value, min_value, sqrt
from icepack2.constants import gravity as g, ice_density as rho_I, water_density as rho_W

__all__ = ["resistive_stress", "critical_stress", "height_above_flotation",
           "buttressing_number", "hfb_calving_rate", "density_in_model_units"]

#: Seconds in icepack's year, which is how its densities carry their units.
_YEAR = 365.25 * 24 * 60 * 60


def density_in_model_units(rho_si):
    r"""A density in kg/m3, which is how the papers quote it, in icepack2's
    own MPa, m, yr. ``density_in_model_units(1024)`` is its water density."""
    return float(rho_si) / _YEAR ** 2 * 1e-6


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


def critical_stress(h, b, sigma_max=0.0, rho_c=None):
    r"""``R_crit`` [MPa]: the resistive stress at which the crevasse field
    penetrates the full thickness.

    ``sigma_max`` is the ice's tensile strength [MPa]; zero is the Buck,
    Coffey and Lai limit, in which an unbuttressed shelf sits exactly at the
    threshold. ``rho_c`` is the density of the water in a basal crevasse
    [kg/m3 in icepack2 units]; seawater by default, fresh water for a
    meltwater-filled crevasse.
    """
    rho_c = rho_W if rho_c is None else rho_c
    thickness = max_value(h, Constant(_H_FLOOR))
    overburden = rho_I * g * thickness
    fraction = max_value(height_above_flotation(h, b), Constant(0.0)) / thickness
    ratio = Constant(float(rho_I / rho_c))
    strength = (Constant(float(rho_c / (rho_c - rho_I)))
                * (Constant(float(sigma_max)) / overburden) ** 2 / 2)
    return overburden * (
        (Constant(1.0) - ratio * (Constant(1.0) - fraction) ** 2) / 2 + strength)


def buttressing_number(M, n, h, b, sigma_max=0.0, rho_c=None):
    r"""Coffey and Lai's ``B = 1 - R_xx / R_crit``, floored at zero: the
    margin a front has left. Zero is a front at the threshold, and the
    theory says an unbuttressed shelf should sit there."""
    r_xx = max_value(resistive_stress(M, n), Constant(0.0))
    return max_value(Constant(1.0) - r_xx / critical_stress(
        h, b, sigma_max, rho_c), Constant(0.0))


def hfb_calving_rate(u, M, h, b, n, sigma_max=0.0, rho_c=None,
                     exponent=1.0, ratio_max=5.0):
    r"""``c = |u| min(R_xx / R_crit, ratio_max)^p`` [m/yr along ``n``]."""
    r_xx = max_value(resistive_stress(M, n), Constant(0.0))
    r_crit = critical_stress(h, b, sigma_max, rho_c)
    ratio = min_value(r_xx / r_crit, Constant(float(ratio_max)))
    speed = sqrt(dot(u, u) + Constant(1e-30))
    return speed * ratio ** Constant(float(exponent))
