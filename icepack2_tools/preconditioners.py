r"""Repository-specific Firedrake preconditioners.

Keep Firedrake imports out of :mod:`icepack2_tools.solverconfig`; PETSc loads
this module lazily only when the corresponding Python PC is selected.
"""

from firedrake.petsc import PETSc
from firedrake.slate.slate import AssembledVector
from firedrake.slate.static_condensation.la_utils import (
    LAContext,
    SchurComplementBuilder,
)
from firedrake.slate.static_condensation.scpc import SCPC


def frozen_linearization(F, z):
    r"""``(J, pre_jacobian_callback)``: a Jacobian of ``F`` that stays at the
    Newton iterate it was formed at, for a matrix-free operator.

    A matrix-free Jacobian is the form's action at whatever the state Function
    holds when it is applied, and Firedrake's residual callback writes every
    point it evaluates into that Function. The NLEQ-ERR line search evaluates
    the residual at a trial point and then solves with the *same* KSP for its
    simplified Newton step, ``J(x_k)^{-1} F(x_trial)``: by then the operator
    has silently become ``J(x_trial)``, while everything SCPC assembled (the
    condensed matrix and its factorization or hierarchy) is still ``x_k``'s.
    An exact condensed MUMPS solve then needs 7-36 outer iterations where it
    should need one (quartz, 2500/25000 x 16: 1134 of a lane's 1246 outer
    iterations), each a mixed-Jacobian action and three Slate sweeps, and the
    step the line search judges is not the one NLEQ-ERR defines. An assembled
    Jacobian (``full_mumps``) has always been frozen by construction.

    The Jacobian is built on a copy of the state that only
    ``pre_jacobian_callback`` refreshes, which Firedrake calls with the
    iterate each time SNES re-forms the Jacobian."""
    from firedrake import Function, derivative
    from ufl import replace

    z_lin = Function(z.function_space(), name="linearization_state")
    # replace() expands the derivative first, so this is J(z_lin), not the
    # derivative of F(z_lin) with respect to a z that is no longer in it.
    J = replace(derivative(F, z), {z: z_lin})

    def pre_jacobian_callback(X):
        with z_lin.dat.vec_wo as v:
            X.copy(v)

    return J, pre_jacobian_callback


def rigid_body_modes(V):
    r"""Orthonormal translations and in-plane rotation of a 2-D vector space:
    the modes a membrane-stress operator without basal drag does not see."""
    from firedrake import (
        Constant, Function, SpatialCoordinate, VectorSpaceBasis, as_vector,
    )

    x, y = SpatialCoordinate(V.mesh())
    modes = (Constant((1.0, 0.0)), Constant((0.0, 1.0)), as_vector((-y, x)))
    basis = VectorSpaceBasis([Function(V).interpolate(mode) for mode in modes])
    basis.orthonormalize()
    return basis


class ISMIP7SCPC(SCPC):
    r"""SCPC for ``(u, M, tau)`` with the retained field ordered first.

    Firedrake's stock three-field helper correctly slices arbitrary eliminated
    fields, but passes their original indices into the newly sliced two-field
    tensor.  Eliminating fields ``1,2`` therefore asks that tensor for blocks
    ``(1,2)`` instead of its renumbered blocks ``(0,1)``.  The stock/default
    order (eliminate ``0,1`` and retain ``2``) does not expose the bug.

    Reordering this model's mixed state would touch checkpoint layout and every
    use of ``z.subfunctions``.  This narrow override keeps the established
    ``(u, M, tau)`` order and changes only the two local indices passed to the
    installed Slate ``SchurComplementBuilder``.  Reconstruction remains the
    upstream implementation and uses the original field numbers.
    """

    def initialize(self, pc):
        super().initialize(pc)
        # Work done on the condensed system since this PC was built: one
        # solve per outer Krylov iteration, and the iterations of those
        # solves (one each under ``preonly``). SNES's own linear-iteration
        # count misses every solve the NLEQ-ERR line search makes for its
        # simplified Newton step -- 62 % of the Krylov work of the first
        # scpc_gamg lane -- so the transient reads its solver work here.
        self.condensed_solves = 0
        self.condensed_iterations = 0

        prefix = (pc.getOptionsPrefix() or "") + "condensed_field_"
        kind = PETSc.Options().getString(prefix + "near_nullspace", "none")
        if kind == "rigid_body":
            # Upstream's ``condensed_field_nullspace`` hook sets a true
            # nullspace, which the KSP would project out of the velocity.
            # GAMG reads the near-nullspace once, at its first setup, which
            # has not happened yet: the condensed KSP is only configured.
            self.near_nullspace = rigid_body_modes(self.weight.function_space())
            _, P = self.condensed_ksp.getOperators()
            P.setNearNullSpace(self.near_nullspace.nullspace())
        elif kind != "none":
            raise ValueError(
                f"{prefix}near_nullspace must be rigid_body or none, not {kind!r}"
            )

    def sc_solve(self, pc):
        super().sc_solve(pc)
        self.condensed_solves += 1
        self.condensed_iterations += self.condensed_ksp.getIterationNumber()

    def condensed_system(self, A, rhs, elim_fields, prefix, pc):
        elim_fields = sorted(map(int, elim_fields))
        if elim_fields != [1, 2]:
            raise ValueError(
                "ISMIP7SCPC requires pc_sc_eliminate_fields=1,2; "
                f"got {elim_fields}"
            )

        blocks = A.blocks
        rhs_blocks = AssembledVector(rhs).blocks
        # Original field 0 (velocity) is retained. Original fields 1:3 become
        # local fields 0:2 after slicing into Aee/Afe/Aef.
        Aff = blocks[0:1, 0:1]
        Aef = blocks[1:3, 0:1]
        Afe = blocks[0:1, 1:3]
        Aee = blocks[1:3, 1:3]
        bf = rhs_blocks[0:1]
        be = rhs_blocks[1:3]

        builder = SchurComplementBuilder(
            prefix,
            Aee,
            Afe,
            Aef,
            pc,
            0,
            1,
            non_zero_saddle_mat=Aff,
        )
        reduced_rhs, reduced_operator = builder.build_schur(
            be, non_zero_saddle_rhs=bf
        )
        return (
            LAContext(
                lhs=reduced_operator,
                rhs=reduced_rhs,
                field_idx=[0],
            ),
            builder,
        )
