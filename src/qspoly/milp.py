r"""Exact search for $(q,s)$-admissible polynomials as a 0--1 integer program.

Key reformulation (Proposition ``ref`` of the paper)
---------------------------------------------------
A multilinear $P=\sum_Tc_T\prod_{i\in T}x_i$ is $(q,s)$-admissible iff the
$2^q$-tuple $g(S):=P(1_S)$ satisfies

* $g(S)\in\{0,1\}$ for all $S$ (then $P$ is Boolean automatically), and $g([q])=0$;
* $c_{[q]}=\sum_{S\subseteq[q]}(-1)^{q-|S|}g(S)\ \ge 1$;
* $\sum_{S\subseteq T}(-1)^{|T|-|S|}g(S)\ \ge 0$ for every $T$ with $|T|>s$;

and for the $g$ obtained this way the actual value of $s=\max\{|T|:c_T<0\}$ equals the
admissible level: feasibility for a level $s$ is monotone in $s$, so the *first*
feasible level is exactly $\max\{|T|:c_T<0\}$ of every witness.

All constraints are linear in the binary variables $g(S)$, so feasibility is decided by
a 0--1 MILP (solved here with HiGHS via :func:`scipy.optimize.milp`).  Because feasibility
is monotone in $s$, the minimal level is found by binary search.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

from .boolean import full, mobius, popcount, subsets_of

__all__ = ["FeasibilityResult", "feasible_at_level", "minimal_level", "search"]


@dataclass
class FeasibilityResult:
    q: int
    s: int
    feasible: bool
    g: Optional[Dict[int, int]]
    status: str
    time_limit_hit: bool = False


def _build_milp(q: int, s: int):
    """Return (csc_matrix, lower, upper, variables) for the MILP at level ``s``."""
    import numpy as np
    from scipy.sparse import lil_matrix
    from scipy.optimize import Bounds, LinearConstraint, milp

    fq = full(q)
    idx: Dict[int, int] = {}
    for S in range(1 << q):
        if S != fq:  # g([q]) = 0 is a structural requirement
            idx[S] = len(idx)
    nv = len(idx)

    rows: list[Tuple[Dict[int, int], float, float]] = []

    # c_[q] >= 1
    rows.append(({S: (-1) ** (q - popcount(S)) for S in idx}, 1.0, float("inf")))
    # c_T >= 0 for |T| > s
    for T in range(1 << q):
        if T == fq:
            continue
        t = popcount(T)
        if t > s:
            rows.append(({S: (-1) ** (t - popcount(S)) for S in subsets_of(T) if S in idx},
                         0.0, float("inf")))

    A = lil_matrix((len(rows), nv))
    for i, (r, lo, hi) in enumerate(rows):
        for S, co in r.items():
            A[i, idx[S]] = co
    A = A.tocsr()
    cons = LinearConstraint(A, np.array([r[1] for r in rows]), np.array([r[2] for r in rows]))
    return milp, cons, nv, idx


def feasible_at_level(q: int, s: int, time_limit: float = 300.0) -> FeasibilityResult:
    """Decide whether a $(q,s)$-admissible polynomial exists (exactly)."""
    milp, cons, nv, idx = _build_milp(q, s)
    import numpy as np
    from scipy.optimize import Bounds

    res = milp(c=np.zeros(nv), constraints=[cons], integrality=np.ones(nv),
               bounds=Bounds(0, 1), options={"time_limit": time_limit, "presolve": True})
    status = {0: "optimal", 1: "time_limit", 2: "infeasible", 3: "unbounded"}.get(
        res.status, f"status_{res.status}")
    ok = res.status == 0 and res.x is not None
    g = None
    if ok:
        g = {S: int(round(res.x[idx[S]])) for S in idx}
        g[full(q)] = 0
    return FeasibilityResult(q=q, s=s, feasible=bool(ok), g=g, status=status,
                             time_limit_hit=(res.status == 1))


def minimal_level(q: int, s_hi: Optional[int] = None, time_limit: float = 300.0
                  ) -> Tuple[Optional[int], Dict[int, FeasibilityResult], bool]:
    """Smallest ``s`` for which level ``s`` is feasible, by binary search.

    Feasibility is monotone: if level $s$ is feasible then level $s+1$ is too (it has
    strictly fewer ``c_T >= 0`` constraints), hence binary search is valid.  Level $0$
    is always infeasible: the constraint sum_T c_T = g([q]) = 0 with all c_T >= 0
    forces all c_T = 0.

    Returns ``(s, log, certified)``.  ``certified`` is ``False`` if any MILP probe hit the
    time limit: a time-limited probe is *unknown*, never "infeasible", so the minimal level
    may then have been over-estimated.
    """
    if s_hi is None:
        s_hi = q - 1
    if s_hi < 0:
        return None, {}, True
    hi = s_hi
    rhi = feasible_at_level(q, hi, time_limit)
    certified = rhi.status in ("optimal", "infeasible")
    if not rhi.feasible:
        return None, {}, certified
    lo = 1  # level 0 is impossible
    log: Dict[int, FeasibilityResult] = {}
    while lo < hi:
        mid = (lo + hi) // 2
        r = feasible_at_level(q, mid, time_limit)
        log[mid] = r
        certified = certified and r.status in ("optimal", "infeasible")
        if r.feasible:
            hi = mid
        else:
            lo = mid + 1
    r = feasible_at_level(q, lo, time_limit)
    log[lo] = r
    certified = certified and r.status in ("optimal", "infeasible")
    return lo, log, certified


def search(q_values, time_limit: float = 300.0) -> list[dict]:
    """Scan ``q`` values; return rows with ``s``, the ratio ``q/s`` and a witness."""
    rows = []
    for q in q_values:
        s, _ = minimal_level(q, time_limit=time_limit)
        if s is None:
            rows.append({"q": q, "s": None, "ratio": None, "g": None})
        else:
            g = feasible_at_level(q, s, time_limit).g
            rows.append({"q": q, "s": s, "ratio": q / s, "g": g})
        print(f"q={q:3d}  s={rows[-1]['s']}  ratio={rows[-1]['ratio']}", flush=True)
    return rows