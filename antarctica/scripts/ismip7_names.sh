# Shell mirror of the MAP filename rule in icepack2_tools/naming.py, for the
# runners and launchers that have to know the name before Python starts.
# Sourced, never executed; defines functions and exports nothing.
#
# Keep these in step with naming.py's _FRICTION_TAGS, map_n_tag and
# map_geom_tag. That module's docstring gives the reason the rule lives in one
# place: every copy of it is a chance for a gate to look for a file the run
# will never write.

# <friction law> -> the filename tag, "" for a law with no tag.
ismip7_friction_tag() {
    case "${1:-}" in
        regularized_coulomb) printf '_rc' ;;
        budd)                printf '_budd' ;;
        *)                   printf '' ;;
    esac
}

# The flow-exponent tag, from ISMIP7_N_FLOW. n=4 keeps the legacy untagged
# name; anything else gets _n<N>, so n=3 and n=4 MAPs coexist on disk.
ismip7_n_tag() {
    awk -v n="${ISMIP7_N_FLOW:-3.0}" \
        'BEGIN { if ((n - 4.0) ^ 2 < 1e-18) printf ""; else printf "_n%d", int(n + 0.5) }'
}

# The geometry-space tag, from ISMIP7_GEOMETRY_SPACE. CG1 keeps the legacy
# untagged name; DG0 gets _dg0.
ismip7_geom_tag() {
    case "$(printf '%s' "${ISMIP7_GEOMETRY_SPACE:-dg0}" | tr '[:upper:]' '[:lower:]')" in
        cg1) printf '' ;;
        *)   printf '_dg0' ;;
    esac
}

# <friction law> <lc> -> the MAP basename the runners write.
#
# The `logvelnet` element records the objective the runners configure (ISSM's
# logarithmic velocity term plus the integrated net mass-balance term) and is
# NOT part of naming.py's rule, so a forward cannot find this MAP by the
# derived name: the runners point ISMIP7_INVERSION at it explicitly.
ismip7_map_basename() {
    printf 'inversion_icepack2%s%s%s_logvelnet_%s.h5' \
        "$(ismip7_friction_tag "${1:-}")" "$(ismip7_n_tag)" "$(ismip7_geom_tag)" "${2:-}"
}
