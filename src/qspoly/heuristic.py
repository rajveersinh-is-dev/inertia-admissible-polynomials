r"""Heuristic search for $(q,s)$-admissible polynomials.

The MILP formulation of :mod:`qspoly.milp` has $2^q-q-1$ binary variables, so exact
feasibility is limited to moderate $q$.  To probe large $q$ we do local search directly on
$g\in\{0,1\}^{2^q}$.

Objective (minimised)
---------------------
$$J(g)\;=\;\sum_{\substack{T\subseteq[q],\,T\neq[q]\\ |T|>s}}\max\bigl(0,\,1-c_T\bigr)
\;+\;W\cdot\max\bigl(0,\,1-c_{[q]}\bigr),\qquad W=4 .$$
In any admissible solution every $c_T$ ($|T|>s$) lies in $\{0,1\}$: indeed
$c_T\ge0$ by admissibility and
$$c_T\;\le\;\sum_{T'\subseteq T}c_{T'}=g(T)\;\le\;1,$$
the first inequality because every other term of $g(T)$ with $|T'|\le s$ is $\le0$ is not
needed -- rather, all $|T'|>s$ terms other than $T$ are $\ge 0$, so this inequality does not
hold in general.  We therefore only *count* violations, and accept ``$J=0$`` as a
certificate only after :func:`verify` recomputes everything from scratch.
"""
from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Dict, List, Optional

from .boolean import full, mobius, popcount, subsets_of

__all__ = ["verify", "heuristic_search", "greedy_warm_start", "SearchResult"]


def verify(g: Dict[int, int], q: int, s: int) -> bool:
    """Recompute $c$ from $g$ and check admissibility conditions at level $s$.

    This is the *only* way a candidate witness is accepted: incremental updates inside the
    search may drift, so the certificate is always recomputed from the stored $g$.
    """
    if set(g.keys()) != set(range(1 << q)):
        return False
    if set(g.values()) - {0, 1}:
        return False
    fq = full(q)
    if g[fq] != 0:
        return False
    c = mobius(q, g)
    if c[fq] < 1:
        return False
    for T in range(1 << q):
        if T != fq and popcount(T) > s and c[T] < 0:
            return False
    return True


def _score(c: Dict[int, int], q: int, s: int, weight: int = 4) -> int:
    """Score.
    
    Args:
        c:
        q:
        s:
        weight (int):
    
    Returns:
        The computed result
    
    """
    J = 0
    fq = full(q)
    for T in range(1 << q):
        if T != fq and popcount(T) > s and c[T] < 0:
            J += 1 - c[T]
    if c[fq] < 1:
        J += weight * (1 - c[fq])
    return J


@dataclass
class SearchResult:
    q: int
    s: int
    feasible: bool
    g: Optional[Dict[int, int]]
    best_J: int
    restarts: int


def heuristic_search(q: int, s: int, restarts: int = 20, iters: int = 20000,
                     seed: int = 0, init: Optional[Dict[int, int]] = None,
                     verbose: bool = False, time_budget: Optional[float] = None
                     ) -> SearchResult:
    r"""Local search for a $g$ with $c_T\ge0$ for all $|T|>s$ and $c_{[q]}\ge1$.

    A move flips one $g(S)$; this changes $c_T$ by $(-1)^{|T|-|S|}\delta$ for every
    $T\supseteq S$, where $\delta=1-2g(S)$.  Supersets of $S$ are enumerated as
    $S\cup U$ for $U\subseteq[q]\setminus S$, costing $2^{q-|S|}$ updates.

    ``time_budget`` (seconds, wall clock) bounds the whole call; the search then reports
    the best value found so far.  A returned ``feasible=True`` has *always* been confirmed
    by :func:`verify`.
    """
    rng = random.Random(seed)
    fq = full(q)
    best_g: Optional[Dict[int, int]] = None
    best_J: Optional[int] = None
    import time as _time
    t_end = None if time_budget is None else _time.time() + time_budget

    for r in range(restarts):
        g = dict(init) if init is not None else {S: rng.getrandbits(1) for S in range(1 << q)}
        g[fq] = 0
        c = mobius(q, g)
        J = _score(c, q, s)
        for _ in range(iters):
            if J == 0:
                break
            if t_end is not None and _time.time() > t_end:
                break
            S = rng.randrange(1, fq)          # never flip g([q])
            delta = 1 - 2 * g[S]
            parity = popcount(S) & 1
            # affected supersets
            sup = [S | U for U in subsets_of(fq ^ S)]
            old = [c[T] for T in sup] + [c[fq]]
            for T in sup:
                d = (popcount(T) & 1) ^ parity
                c[T] += -delta if (d & 1) else delta
            c[fq] += -delta if ((q & 1) ^ parity) else delta
            newJ = _score(c, q, s)
            if newJ <= J:
                g[S] = 1 - g[S]
                J = newJ
            else:
                for T, v in zip(sup, old[:
                    -1]):
                    c[T] = v
                c[fq] = old[-1]
        if best_J is None or J < best_J:
            best_J, best_g = J, dict(g)
            if verbose:
                print(f"   restart {r}: J={J}", flush=True)
        if best_J == 0 and verify(best_g, q, s):
            return SearchResult(q, s, True, best_g, 0, r + 1)
    ok = bool(best_g is not None and verify(best_g, q, s))
    return SearchResult(q, s, ok, best_g if ok else None,
                        int(best_J if best_J is not None else -1), restarts)


def greedy_warm_start(q: int, s: int) -> Optional[Dict[int, int]]:
    """Structured warm starts: ``g(S)=1`` iff ``|S| <= t`` (or its complement-free variants)."""
    fq = full(q)
    best, bestJ = None, None
    for t in range(0, q + 1):
        for flip in (False, True):
            g = {S: (1 if (popcount(S) <= t) else 0) for S in range(1 << q)}
            if flip:
                g = {S: 1 - g[S] for S in range(1 << q)}
            g[fq] = 0
            c = mobius(q, g)
            J = _score(c, q, s)
            if bestJ is None or J < bestJ:
                bestJ, best = J, g
            if J == 0 and verify(g, q, s):
                return g
    return best