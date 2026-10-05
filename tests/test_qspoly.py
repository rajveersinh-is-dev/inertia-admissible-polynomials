"""Test suite for :mod:`qspoly`.

Run with ``pytest tests/`` or directly with ``python tests/test_qspoly.py``.
Every test that certifies a mathematical statement recomputes the relevant quantity from
scratch; nothing trusts an incremental computation.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import numpy as np

from qspoly import (AdmissiblePolynomial, coefficients_to_g, from_g, mobius, zeta,
                    full, popcount, subsets_of)
from qspoly.graphs import (basis, leading_sign, neps_adjacency, neps_eigenvalues,
                           neps_inertia)
from qspoly.heuristic import greedy_warm_start, heuristic_search, verify
from qspoly.milp import feasible_at_level
from qspoly.exact import exact_inertia


def _published_7_3() -> AdmissiblePolynomial:
    def m(ts):
        return sum(1 << (i - 1) for i in ts)
    c = {}
    for t in [(1, 3), (2, 3), (1, 4), (2, 4), (5, 6), (5, 7)]:
        c[m(t)] = 1
    for t in [(1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4),
              (1, 5, 6), (2, 5, 6), (3, 5, 7), (4, 5, 7), (5, 6, 7)]:
        c[m(t)] = -1
    for t in [(1, 2, 5, 6), (3, 4, 5, 7)]:
        c[m(t)] = 1
    c[m(tuple(range(1, 8)))] = 1
    return from_g(7, coefficients_to_g(7, c))


# ------------------------------------------------------------------ definitions
def test_subsets_of_enumeration():
    for T in [0, 1, 5, 7, 12345]:
        subs = list(subsets_of(T))
        assert len(subs) == 2 ** popcount(T)
        assert set(subs) == set(subsets_of(T))
        for S in subs:
            assert (S & ~T) == 0


def test_zeta_mobius_inverse():
    import random
    rng = random.Random(7)
    for q in range(1, 7):
        g = {S: rng.randint(0, 1) for S in range(1 << q)}
        c = mobius(q, g)
        assert zeta(q, c) == g


def test_g_defines_boolean_polynomial():
    P = _published_7_3()
    assert P.is_boolean()
    for S in range(1 << 7):
        x = [(S >> i) & 1 for i in range(7)]
        assert P.evaluate(x) == P.g[S]


# ------------------------------------------------------------------ admissibility
def test_published_7_3_is_admissible():
    P = _published_7_3()
    assert P.g[full(7)] == 0                    # P(1,...,1) = 0
    assert P.c[full(7)] > 0                      # c_[7] = 1 > 0
    assert P.max_negative_level() == 3
    assert P.is_admissible()
    assert abs(P.qs_ratio() - 7 / 3) < 1e-12


def test_simple_2_1_is_admissible():
    P = from_g(2, {0: 1, 1: 0, 2: 0, 3: 0})
    assert P.is_admissible()
    assert P.max_negative_level() == 1
    assert P.c[full(2)] == 1


def test_level_zero_is_impossible():
    """sum_T c_T = g([q]) = 0 with all c_T >= 0 forces c = 0, so c_[q] > 0 is impossible."""
    import random
    rng = random.Random(11)
    for q in range(1, 5):
        for _ in range(50):
            g = {S: rng.randint(0, 1) for S in range(1 << q)}
            g[full(q)] = 0
            c = mobius(q, g)
            if c[full(q)] >= 1:
                assert any(c[T] < 0 for T in c)


def test_s1_forces_q2_exhaustive_small_q():
    """Exhaustive check for q <= 3 that s = 1 occurs only at q = 2, and is then unique.

    The unique witness is $g(\\emptyset)=1$, $g(S)=0$ otherwise, i.e. $P=(1-x_1)(1-x_2)$;
    with ``free = [0,1,2]`` that is ``bits = 1``.
    """
    found = []
    for q in range(1, 4):
        fq = full(q)
        free = [S for S in range(1 << q) if S != fq]
        for bits in range(1 << len(free)):
            g = {fq: 0}
            for i, S in enumerate(free):
                g[S] = (bits >> i) & 1
            P = from_g(q, g)
            if P.is_admissible() and P.max_negative_level() == 1:
                found.append((q, bits))
    assert found == [(2, 1)], found
    g = {0: 1, 1: 0, 2: 0, 3: 0}
    assert from_g(2, g).polynomial_expression().replace(" ", "") in ("1-1*x1-1*x2+1*x1*x2",
                                                                    "1-x1-x2+x1*x2")


def test_monotonicity_of_level_feasibility():
    """A witness at level s is automatically a witness at every level s' > s."""
    P = _published_7_3()
    c, fq = P.c, full(7)
    s = P.max_negative_level()
    for sp in range(s, 7):
        assert all(c[T] >= 0 for T in c if T != fq and popcount(T) > sp)


def test_restriction_preserves_coefficients():
    """c_T depends only on g restricted to 2^T (this is why level feasibility is monotone in q)."""
    P = _published_7_3()
    c = P.c
    for Q in [0, 1, 3, 5, 7]:
        local = mobius(Q, {S: P.g[S] for S in range(1 << Q)})
        for T in range(1 << Q):
            assert local[T] == c[T], (Q, T)


# ------------------------------------------------------------------ spectra
def test_eigenvalue_formula_totals_order():
    """sum_T (m-1)^{|T|} = m^q, i.e. the multiplicities account for all vertices."""
    P = _published_7_3()
    for m in [2, 3, 4]:
        tot = sum((m - 1) ** popcount(T) for T in range(1 << 7))
        assert tot == m ** 7


def test_neps_inertia_is_consistent_and_non_singular():
    P = _published_7_3()
    for m in [2, 3, 4, 5]:
        p, z, n = neps_inertia(P, m)
        assert p + z + n == m ** 7
        assert z == 0, "admissible polynomials give non-singular NEPS graphs"
        assert n > 0 and p > 0


def test_neps_inertia_matches_numerical_and_exact():
    P = _published_7_3()
    A = neps_adjacency(P, 2)
    assert A.shape == (128, 128)
    assert set(np.unique(A)) <= {0, 1}
    assert np.all(np.diag(A) == 0)
    w = np.linalg.eigvalsh(A.astype(float))
    assert int(np.sum(w > 1e-7)) == neps_inertia(P, 2)[0]
    assert int(np.sum(w < -1e-7)) == neps_inertia(P, 2)[2]
    # exact inertia of a small case
    Q = from_g(2, {0: 1, 1: 0, 2: 0, 3: 0})
    B = neps_adjacency(Q, 3)
    assert exact_inertia([[int(x) for x in row] for row in B]) == (5, 0, 4)


def test_neps_adjacency_matches_definition():
    """Brute-force NEPS adjacency for a small case, independent of the Kronecker formula."""
    Q = from_g(2, {0: 1, 1: 0, 2: 0, 3: 0})
    m, F = 4, basis(Q)
    verts = [(i, j) for i in range(m) for j in range(m)]
    A = neps_adjacency(Q, m)
    for ui, u in enumerate(verts):
        for vi, v in enumerate(verts):
            if ui == vi:
                continue
            adj = False
            for S in F:
                if all((u[i] == v[i]) if i not in S else (u[i] != v[i]) for i in range(2)):
                    adj = True
                    break
            assert bool(A[ui, vi]) == adj


def test_trace_moments():
    """Exact integer identity tr(A^k) = sum_T lambda_T^k (m-1)^{|T|}."""
    from collections import Counter
    Q = from_g(2, {0: 1, 1: 0, 2: 0, 3: 0})
    m, K = 3, 5
    A = neps_adjacency(Q, m).astype(np.int64)
    lam = neps_eigenvalues(Q, m)
    mult = Counter()
    for T, v in lam.items():
        mult[v] += (m - 1) ** popcount(T)
    Ak = np.eye(A.shape[0], dtype=object)
    for k in range(1, K + 1):
        Ak = Ak @ A
        assert int(np.trace(Ak)) == sum(v ** k * c for v, c in mult.items())


def test_leading_sign_matches_finite_m_for_large_m():
    P = _published_7_3()
    m = 5000            # far beyond the asymptotic threshold
    lam = neps_eigenvalues(P, m)
    for T in range(1 << 7):
        s = leading_sign(P, T)
        v = lam[T]
        assert (v > 0) - (v < 0) == s, (T, v, s)


# ------------------------------------------------------------------ search
def test_milp_witness_passes_independent_check():
    for q, s in [(4, 2), (5, 3), (6, 3), (7, 3)]:
        r = feasible_at_level(q, s, time_limit=120)
        assert r.feasible, (q, s, r.status)
        assert verify(r.g, q, s)
        P = from_g(q, r.g)
        assert P.max_negative_level() == s
        assert P.is_admissible()


def test_milp_certifies_known_infeasibility():
    """s = 2 at q = 6 is exactly infeasible (proved in the paper, 3s = q boundary)."""
    r = feasible_at_level(6, 2, time_limit=180)
    assert r.status == "infeasible", r.status


def test_heuristic_witness_is_reverified():
    q, s = 7, 3
    ws = greedy_warm_start(q, s)
    h = heuristic_search(q, s, restarts=3, iters=8000, seed=3, init=ws)
    if h.feasible:                      # must always survive from-scratch verification
        assert verify(h.g, q, s)


def test_verify_rejects_malformed_or_inadmissible_candidates():
    """``verify`` must reject each failure mode of admissibility."""
    good = {0: 1, 1: 0, 2: 0, 3: 0}
    assert verify(good, 2, 1)
    # (a) g([q]) != 0
    assert not verify({0: 1, 1: 0, 2: 0, 3: 1}, 2, 1)
    # (b) not {0,1}-valued
    assert not verify({0: 2, 1: 0, 2: 0, 3: 0}, 2, 1)
    # (c) missing entries
    assert not verify({0: 1, 1: 0, 2: 0}, 2, 1)
    # (d) c_[q] < 1 : take g = 1 - (indicator of {[2]}) ... use g(emptyset)=1, g({1})=1
    g2 = {0: 1, 1: 1, 2: 0, 3: 0}
    c2 = mobius(2, g2)
    assert c2[3] <= 0
    assert not verify(g2, 2, 1)
    # (e) a negative coefficient above the declared level
    P = _published_7_3()
    s = P.max_negative_level()
    assert not verify(P.g, 7, s - 1)      # lowering the level re-exposes a negative c_T


if __name__ == "__main__":
    fns = [(k, v) for k, v in sorted(globals().items())
           if k.startswith("test_") and callable(v)]
    npass = nfail = 0
    for name, fn in fns:
        try:
            fn()
            print(f"PASS  {name}")
            npass += 1
        except Exception as exc:      # noqa: BLE001
            print(f"FAIL  {name}: {type(exc).__name__}: {exc}")
            nfail += 1
    print(f"\n{npass} passed, {nfail} failed")
    sys.exit(1 if nfail else 0)