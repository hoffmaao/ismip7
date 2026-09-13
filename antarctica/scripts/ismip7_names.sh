# Shell mirror of the MAP filename rule in icepack2_tools/naming.py, for the
# runners and launchers that have to know the name before Python starts.
# Sourced, never executed; defines functions and exports nothing.
#
# Keep the friction tags in step with naming.py's _FRICTION_TAGS. That module's
# docstring gives the reason the rule lives in one place: every copy of it is a
# chance for a gate to look for a file the run will never write.

# <friction law> -> the filename tag, "" for a law with no tag.
ismip7_friction_tag() {
    case "${1:-}" in
        regularized_coulomb) printf '_rc' ;;
        budd)                printf '_budd' ;;
        *)                   printf '' ;;
    esac
}

# <friction law> <lc> -> the MAP basename the runners write.
#
# The `logvelnet` element records the objective the runners configure (ISSM's
# logarithmic velocity term plus the integrated net mass-balance term) and is
# NOT part of naming.py's rule, so a forward cannot find this MAP by the
# derived name: the runners point ISMIP7_INVERSION at it explicitly.
ismip7_map_basename() {
    printf 'inversion_icepack2%s_n3_dg0_logvelnet_%s.h5' \
        "$(ismip7_friction_tag "${1:-}")" "${2:-}"
}
