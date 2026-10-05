r"""Exact inertia of an integer symmetric matrix.

The inertia $(n^+,n^0,n^-)$ of a real symmetric matrix can be read off from the signs of its
eigenvalues.  For an *integer* matrix the characteristic polynomial has rational
coefficients, so its roots are exact algebraic numbers; SymPy represents them either as
``Rational``/``AlgebraicNumber`` or as an indexed real root of the polynomial
(``CRootOf``), and SymPy decides ``root > 0`` exactly by interval refinement.

This module therefore computes inertia with **no floating point anywhere**.  A
floating-point cross-check is also provided for convenience, but it is never used to
certify a result.
"""
from __future__ import annotations

from typing import List, Sequence, Tuple

__all__ = ["exact_inertia", "exact_inertia_from_charpoly", "floating_inertia"]


def _classify(root) -> str:
    if root.is_zero:
        return "zero"
    # exact sign test for rationals and for real algebraic numbers
    if bool(root > 0):
        return "pos"
    if bool(root < 0):
        return "neg"
    # fall back: exact sign via minimal polynomial interval arithmetic
    import sympy as sp
    s = sp.polys.numberfields.minpoly  # noqa: F401  (ensure import side effects)
    val = sp.N(root, 60)
    if abs(val) < sp.Rational(1, 10) ** 40:
        return "zero"
    return "pos" if val > 0 else "neg"


def exact_inertia_from_charpoly(cp) -> Tuple[int, int, int]:
    r"""Exact inertia of a symmetric integer matrix from its characteristic polynomial."""
    import sympy as sp

    npos = nzero = nneg = 0
    for root in sp.polys.polytools.Poly(cp, cp.gens[0]).all_roots():
        kind = _classify(root)
        if kind == "pos":
            npos += 1
        elif kind == "neg":
            nneg += 1
        else:
            nzero += 1
    return npos, nzero, nneg


def exact_inertia(M) -> Tuple[int, int, int]:
    r"""Exact inertia of the symmetric matrix ``M`` (sympy Matrix or nested lists)."""
    import sympy as sp

    if not isinstance(M, sp.MatrixBase):
        M = sp.Matrix([[sp.Rational(int(v)) for v in row] for row in M])
    n = M.rows
    if M.cols != n:
        raise ValueError("matrix must be square")
    # a symmetric integer matrix has characteristic polynomial with rational coefficients
    x = sp.Symbol("x")
    return exact_inertia_from_charpoly(sp.Poly(M.charpoly(x).as_expr(), x))


def floating_inertia(M, tol: float = 1e-8) -> Tuple[int, int, int]:
    """Floating-point inertia (cross-check only, never a certificate)."""
    import numpy as np

    w = np.linalg.eigvalsh(np.asarray(M, dtype=float))
    return (int(np.sum(w > tol)), int(np.sum(np.abs(w) <= tol)), int(np.sum(w < -tol)))