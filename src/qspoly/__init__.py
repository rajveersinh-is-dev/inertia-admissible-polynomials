"""(q,s)-admissible Boolean polynomials and the graphs they produce."""
from .boolean import (  # noqa: F401
    AdmissiblePolynomial,
    coefficients_to_g,
    from_g,
    full,
    is_admissible,
    mobius,
    popcount,
    subsets_of,
    zeta,
)
from .graphs import (  # noqa: F401
    basis,
    leading_sign,
    neps_adjacency,
    neps_eigenvalues,
    neps_inertia,
    sign_profile,
)

__all__ = [
    "AdmissiblePolynomial", "from_g", "mobius", "zeta", "coefficients_to_g",
    "full", "popcount", "subsets_of", "negative_levels", "is_admissible",
    "basis", "neps_eigenvalues", "leading_sign", "neps_inertia", "sign_profile",
    "neps_adjacency",
]
__version__ = "1.0.0"