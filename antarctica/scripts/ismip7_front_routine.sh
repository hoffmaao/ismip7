#!/bin/bash
# The ISMIP7 submission's front, sourced by a runner before the core drivers
# (Rice, 26 Sep 2026: "for the ismip7 submission we need to focus on just one
# simulation" - the resistive-stress law).
#
#   * calving: the horizontal-force-balance (resistive-stress) law from
#     icepack_tools.calving, at the Slater & Wagner tensile strength of
#     150 kPa (sigma_max=0.15 MPa; 0 is the Buck / Coffey & Lai limit);
#   * advance: none. The law retreats the front; ice the transport carries
#     past the t=0 extent is removed each step and booked as calving, the
#     retreat-only front of most ISMIP6 models (Seroussi et al. 2020);
#   * ice-shelf collapse: the protocol's path C mask on the front cells
#     (ISMIP7_FRACTURE=mask_front) wherever a scenario ships one; the
#     historical trees carry none, so a historical runs with `none`.
# The friction is the MAP's (sub-element grounding when the MAP records it);
# the fixed-front flag is irrelevant under a law and left unset here.
export ISMIP7_CALVING=hfb
export ISMIP7_CALVING_PARAMS="${ISMIP7_CALVING_PARAMS:-sigma_max=0.15}"
export ISMIP7_FRONT_ADVANCE=none
export ISMIP7_FRACTURE="${ISMIP7_FRACTURE:-mask_front}"
unset ISMIP7_FIXED_FRONT
