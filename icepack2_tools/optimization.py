r"""Stopping rules for the inversions, shared by the TAO and scipy paths."""


class FunctionalDecreaseStop:
    r"""Relative-decrease stopping rule, for a TAO monitor to apply.

    PETSc's TAO tests only gradient norms. This stops on the relative decrease
    of the functional between successive iterates,

        (J_old - J_new) / max(|J_old|, |J_new|, 1) <= ftol,

    and never before ``min_iter`` iterations. scipy's L-BFGS-B ``ftol`` is the
    same expression, so both optimizer paths read one knob. 1e-10 converges a
    production inversion; the L-surface sweeps pass 1e-4.
    """

    def __init__(self, ftol, min_iter=3):
        self.ftol = float(ftol)
        self.min_iter = int(min_iter)
        self.criterion = None
        self._last = None

    def update(self, it, f):
        r"""Record iterate ``it``'s functional ``f``; True when the rule says stop."""
        f = float(f)
        if self._last is None:
            self._last = f
            return False
        self.criterion = (self._last - f) / max(abs(self._last), abs(f), 1.0)
        self._last = f
        return self.ftol > 0 and it >= self.min_iter and self.criterion <= self.ftol
