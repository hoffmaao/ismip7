r"""Small transformations of UFL forms."""
import ufl


def with_quadrature_degree(form, degree):
    r"""``form`` with ``quadrature_degree`` set in every integral's metadata
    (an integral that already names a degree keeps it)."""
    return ufl.Form([
        itg.reconstruct(metadata={"quadrature_degree": degree, **itg.metadata()})
        for itg in form.integrals()])
