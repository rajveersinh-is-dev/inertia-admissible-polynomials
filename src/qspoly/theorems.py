r"""Proved results about $(q,s)$-admissible polynomials.

Everything in this module is stated and proved in ``paper/main.tex``; the code here is a
machine-checkable shadow of the proofs (each theorem exposes a predicate that the test
suite evaluates on exhaustive small instances).
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

from .boolean import AdmissiblePolynomial, full, mobius, popcount, subsets_of

__all__ = [
    "LEVEL0_IMPOSSIBLE",
    "s1_only_q2",
    "monotone_in_s",
    "restriction_preserves_coefficients",
    "certified_table",
]


def LEVEL0_IMPOSSIBLE(q: int, g: Dict[int, int]) -> bool:
    r"""Level $0$ is never admissible (proof of Prop. ``level0``).

    If $c_T\ge0$ for every $T$ then $\sum_Tc_T=g([q])=0$ forces all $c_T=0$, contradicting
    $c_{[q]}>0$.  Hence every admissible pair has $s\ge1$.
    """
    c = mobius(q, g)
    if c[full(q)] < 1:
        return False                      # not a candidate at all
    return any(c[T] < 0 for T in c)      # some coefficient must be negative


def s1_only_q2(g: Dict[int, int], q: int) -> bool:
    r"""Theorem ``s-one``: an admissible polynomial must have $s=1$ only if $q=2$.

    Return ``True`` iff $g$ is admissible with $\max\{|T|:c_T<0\}=1$; the *validity* of the
    theorem (that no admissible $g$ has $s=1$ and $q\ge3$) is checked exhaustively in the
    test suite for small $q$ and follows from the case analysis in the paper.
    """
    P = AdmissiblePolynomial(q=q, g=g)
    return P.is_admissible() and P.max_negative_level() == 1 and q == 2


def s1_unique_polynomial(q: int = 2) -> Dict[int, int]:
    r"""The unique $(2,1)$-admissible polynomial: $P=(1-x_1)(1-x_2)$, so $g(\emptyset)=1$."""
    fq = full(q)
    return {S: (1 if S == 0 else 0) for S in range(1 << q)}


def monotone_in_s(g: Dict[int, int], q: int) -> bool:
    r"""If $g$ is admissible with level $s$ then it is admissible at every level $s'>s$.

    (Feasibility at level $s'$ only demands $c_T\ge0$ for the *fewer* subsets $|T|>s'$.)
    """
    P = AdmissiblePolynomial(q=q, g=g)
    s = P.max_negative_level()
    c = P.c
    fq = full(q)
    for sp in range(s, q):
        if any(c[T] < 0 for T in c if T != fq and popcount(T) > sp):
            return False
    return True


def restriction_preserves_coefficients(g: Dict[int, int], q: int, Q: int) -> bool:
    r"""``c_T`` only depends on $g$ restricted to subsets of $T$.

    Hence a witness for $(q,s)$ restricted to any $Q\subseteq[q]$ is a witness for
    $(|Q|,s)$ whenever $|Q|>s$.  (This is what makes level-feasibility monotone in $q$.)
    """
    c = mobius(q, g)
    for T in range(1 << Q):
        local = mobius(Q, {S: g[S] for S in range(1 << Q)})
        if local[T] != c[T]:
            return False
    return True


def certified_table() -> List[Tuple[int, Optional[int]]]:
    """The certified ``(q, s_min(q))`` pairs from ``results/scan_certified.json``."""
    import json
    import os
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    path = os.path.join(root, "results", "scan_certified.json")
    if not os.path.exists(path):
        return []
    rows = json.load(open(path))
    return [(r["q"], r["s"] if r.get("certified") else None) for r in rows]