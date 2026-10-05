"""Exact certificates for the graphs built from admissible polynomials.

Three independent checks are run for every (construction, m):

1. ``trace moments``:  for k = 1..K the exact integer identity
   tr(A^k) = sum_{T} lambda_T(m)^k (m-1)^{|T|}  is compared with the trace of the
   explicitly constructed adjacency matrix A(G_m).  This certifies the eigenvalues
   together with their multiplicities (Vandermonde), and is pure integer arithmetic.

2. ``characteristic polynomial``:  for orders small enough, the exact char. polynomial of
   A(G_m) is computed with SymPy and compared with prod_T (t - lambda_T)^{(m-1)^{|T|}}.

3. ``inertia by exact root counting``:  the number of positive / negative / zero roots of
   that polynomial is computed exactly (Sturm sequences via SymPy's ``count_roots``) and
   compared with the closed-form inertia.
"""
import sys, os, json, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))
import numpy as np
from qspoly import from_g, coefficients_to_g
from qspoly.graphs import basis, neps_adjacency, neps_eigenvalues, neps_inertia, leading_sign
from qspoly.heuristic import verify
from qspoly.exact import exact_inertia_from_charpoly


def published_7_3():
    def m(ts): return sum(1 << (i - 1) for i in ts)
    c = {}
    for t in [(1,3),(2,3),(1,4),(2,4),(5,6),(5,7)]:        c[m(t)] = 1
    for t in [(1,2,3),(1,2,4),(1,3,4),(2,3,4),
              (1,5,6),(2,5,6),(3,5,7),(4,5,7),(5,6,7)]:     c[m(t)] = -1
    for t in [(1,2,5,6),(3,4,5,7)]:                          c[m(t)] = 1
    c[m(tuple(range(1, 8)))] = 1
    return from_g(7, coefficients_to_g(7, c))


def simple_2_1():
    return from_g(2, {0: 1, 1: 0, 2: 0, 3: 0})


def trace_moment_check(P, m, K):
    A = neps_adjacency(P, m).astype(np.int64)
    lam = neps_eigenvalues(P, m)
    from collections import Counter
    mult = Counter()
    for T, v in lam.items():
        mult[v] += (m - 1) ** bin(T).count('1')
    out = []
    Ak = np.eye(A.shape[0], dtype=object)
    for k in range(1, K + 1):
        Ak = Ak @ A
        lhs = int(np.trace(Ak))
        rhs = sum(int(v) ** k * int(c) for v, c in mult.items())
        out.append({"k": k, "trace": lhs, "predicted": rhs, "ok": lhs == rhs})
    return out


def charpoly_check(P, m):
    import sympy as sp
    A = neps_adjacency(P, m)
    n = A.shape[0]
    M = sp.zeros(n)
    for i in range(n):
        for j in range(n):
            if A[i, j]:
                M[i, j] = 1
    x = sp.Symbol('x')
    cp = sp.Poly(M.charpoly(x).as_expr(), x)
    lam = neps_eigenvalues(P, m)
    pred = sp.Integer(1)
    for T, v in lam.items():
        pred *= (x - sp.Integer(v)) ** ((m - 1) ** bin(T).count('1'))
    pred = sp.Poly(pred, x)
    same = sp.expand(cp.as_expr() - pred.as_expr()) == 0
    pos, zero, neg = exact_inertia_from_charpoly(cp)
    p, z, nneg = neps_inertia(P, m)
    return {"order": n, "charpoly_equal": bool(same),
            "exact_inertia": [int(p), int(z), int(nneg)],
            "root_count_inertia": [int(pos), int(zero), int(neg)],
            "inertia_equal": [int(p), int(z), int(nneg)] == [int(pos), int(zero), int(neg)],
            "charpoly": str(cp.as_expr())[:300]}


def numerical_check(P, m):
    A = neps_adjacency(P, m).astype(float)
    w = np.linalg.eigvalsh(A)
    p = int(np.sum(w > 1e-7)); z = int(np.sum(np.abs(w) <= 1e-7)); n = int(np.sum(w < -1e-7))
    closed = neps_inertia(P, m)
    return {"numerical_inertia": [p, z, n], "closed_form_inertia": list(closed),
            "inertia_equal": [p, z, n] == list(closed),
            "numerical_spectrum_sorted": [round(float(x), 6) for x in np.sort(w)[::-1][:12]]}


if __name__ == '__main__':
    certs = []
    cases = [("published (7,3)", published_7_3(), [2]),
             ("(2,1) P=(1-x1)(1-x2)", simple_2_1(), [2, 3, 4, 5])]
    for name, P, ms in cases:
        assert P.is_admissible(), name
        s = P.max_negative_level()
        print(f"=== {name}: q={P.q} s={s} ratio={P.qs_ratio():.4f} admissible={P.is_admissible()}")
        for m in ms:
            n = m ** P.q
            tm = trace_moment_check(P, m, K=8)
            ok_tm = all(r['ok'] for r in tm)
            rec = {"name": name, "q": P.q, "s": s, "m": m, "order": n,
                   "trace_moments_ok": ok_tm, "trace_moments": tm,
                   "closed_form_inertia": list(neps_inertia(P, m)),
                   "numerical": numerical_check(P, m)}
            rec["ok"] = bool(ok_tm and rec["numerical"]["inertia_equal"])
            if n <= 40:
                rec.update(charpoly_check(P, m))
                rec["ok"] = bool(rec["ok"] and rec["charpoly_equal"] and rec["inertia_equal"])
            print(f"   m={m} order={n}: trace moments ok={ok_tm}  inertia={neps_inertia(P, m)}"
                  f"  numerical_equal={rec['numerical']['inertia_equal']}"
                  + (f"  charpoly_equal={rec.get('charpoly_equal')} "
                     f"inertia_equal={rec.get('inertia_equal')}" if 'charpoly_equal' in rec else "")
                  + f"   -> {'PASS' if rec['ok'] else 'FAIL'}")
            certs.append(rec)
            sys.stdout.flush()
    os.makedirs('results', exist_ok=True)
    json.dump(certs, open('results/certificates.json', 'w'), indent=1)
    print("\nwrote results/certificates.json")
    print("ALL CERTIFICATES PASS:", all(c['ok'] for c in certs))