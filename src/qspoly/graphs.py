r"""The graph construction: NEPS graphs from a $(q,s)$-admissible polynomial.

For a Boolean multilinear polynomial
$$P(x)=\sum_{T\subseteq[q]}c_T\prod_{i\in T}x_i \in\mathbb Z[x_1,\dots,x_q],$$
define the *support* of $P$ on the cube by
$$\mathcal X_P=\{x\in\{0,1\}^q: P(x)=1\}.$$
Since $P(\mathbf 1)=0$, $(0,\dots,0)\notin\mathcal X_P$.  Set
$$\mathcal F_P=\{S\subseteq[q]: P(1_{[q]\setminus S})=1\},\qquad
 \mathcal B_P=\{1_S: S\in\mathcal F_P\}.$$
Then $G_m=\mathrm{NEPS}(K_m,\dots,K_m;\mathcal B_P)$ (with $q$ factors) has order $m^q$, and
$$\lambda_T(m)=\sum_{T\subseteq R\subseteq[q]}c_R\,m^{q-|R|}
 \quad\text{with multiplicity}\quad (m-1)^{|T|}. \tag{$\ast$}$$
Every $\lambda_T(m)$ is an integer, so the inertia is computed exactly.

The sign of $\lambda_T(m)$ for all sufficiently large $m$ is the sign of the leading
non-zero coefficient $c_R$ over the $R\supseteq T$ with $c_R\neq0$ of smallest $|R|$;
in particular the graph is non-singular for all large $m$, and
$n^+(G_m)=\Theta(m^q)$, $n^-(G_m)=\Theta(m^s)$.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

from .boolean import AdmissiblePolynomial, full, popcount, subsets_of

__all__ = [
    "basis",
    "neps_eigenvalues",
    "leading_sign",
    "neps_inertia",
    "sign_profile",
    "neps_adjacency",
]


def basis(P: AdmissiblePolynomial) -> List[Tuple[int, ...]]:
    r"""$\mathcal F_P$ as a list of 0-based index tuples."""
    q = P.q
    fq = full(q)
    return [tuple(i for i in range(q) if (S >> i) & 1)
            for S in range(1 << q) if P.g[fq ^ S] == 1]


def neps_eigenvalues(P: AdmissiblePolynomial, m: int) -> Dict[int, int]:
    r"""``lambda_T(m)`` for every ``T subset [q]``, from $(\ast)$."""
    q, c, fq = P.q, P.c, full(P.q)
    out: Dict[int, int] = {}
    for T in range(1 << q):
        rest = fq ^ T
        val = 0
        for U in subsets_of(rest):
            R = T | U
            val += c.get(R, 0) * (m ** (q - popcount(R)))
        out[T] = val
    return out


def leading_sign(P: AdmissiblePolynomial, T: int) -> int:
    r"""Sign of $\lambda_T(m)$ for all sufficiently large $m$."""
    c, fq = P.c, full(P.q)
    cands = [R for R in subsets_of(fq) if (R & T) == T and c.get(R, 0) != 0]
    if not cands:
        return 0
    R = min(cands, key=popcount)
    return (c[R] > 0) - (c[R] < 0)


def neps_inertia(P: AdmissiblePolynomial, m: int) -> Tuple[int, int, int]:
    r"""Exact inertia ``(n^+, n^0, n^-)`` of ``G_m``, from $(\ast)$."""
    p = z = n = 0
    for T, v in neps_eigenvalues(P, m).items():
        mult = (m - 1) ** popcount(T)
        if v > 0:
            p += mult
        elif v < 0:
            n += mult
        else:
            z += mult
    return p, z, n


def sign_profile(P: AdmissiblePolynomial) -> Dict[int, int]:
    r"""``{(sign, |T|): number of T}`` for all sufficiently large $m$."""
    prof: Dict[int, int] = {}
    for T in range(1 << P.q):
        prof[leading_sign(P, T)] = prof.get(leading_sign(P, T), 0) + 1
    return prof


def _smallest_m_non_singular(P: AdmissiblePolynomial, m_max: int = 40) -> int:
    """Smallest ``m`` in ``[2, m_max]`` with ``lambda_T(m) != 0`` for all ``T``."""
    for m in range(2, m_max + 1):
        if all(v != 0 for v in neps_eigenvalues(P, m).values()):
            return m
    raise ValueError(f"no non-singular m found up to {m_max}")


def neps_adjacency(P: AdmissiblePolynomial, m: int):
    r"""Dense integer adjacency matrix of ``G_m`` via
    $A(G_m)=\sum_{S\in\mathcal F_P}\bigotimes_{i=1}^q M_i(S)$,
    $M_i(S)=J_m-I_m$ if $i\in S$, $I_m$ otherwise."""
    import numpy as np

    F = basis(P)
    n = m ** P.q
    Jm = np.ones((m, m), dtype=np.int64) - np.eye(m, dtype=np.int64)
    Im = np.eye(m, dtype=np.int64)
    A = np.zeros((n, n), dtype=np.int64)
    for S in F:
        term = np.ones((1, 1), dtype=np.int64)
        for i in range(P.q):
            term = np.kron(term, Jm if i in S else Im)
        A += term
    A = np.minimum(A, 1)
    np.fill_diagonal(A, 0)
    return A