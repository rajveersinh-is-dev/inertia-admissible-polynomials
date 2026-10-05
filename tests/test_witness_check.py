"""Independent re-verification of witnesses (fast version: heuristic budget trimmed).

The exhaustive/expensive checks live in test_qspoly.py and exp_certify_graphs.py; this
script exists to show that MILP and heuristic witnesses agree with a *third*, brute-force
recomputation of the Moebius transform that shares no code with qspoly.
"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))
from qspoly import from_g
from qspoly.heuristic import heuristic_search, greedy_warm_start, verify
from qspoly.milp import feasible_at_level


def independent_check(g, q, s):
    """Re-derive c_T from g by brute force; shares no code with qspoly.boolean."""
    def subsets(T):
        r = 0
        while True:
            yield r
            if r == T:
                return
            r = (r - T) & T
    c = {}
    for T in range(1 << q):
        c[T] = sum((-1) ** (bin(T).count('1') - bin(U).count('1')) * g[U]
                   for U in subsets(T))
    full = (1 << q) - 1
    ok = (g[full] == 0 and c[full] >= 1 and
          all(c[T] >= 0 for T in c if T != full and bin(T).count('1') > s))
    return ok, c


if __name__ == '__main__':
    for (q, s) in [(4, 2), (5, 3), (6, 3), (7, 3), (7, 2)]:
        r = feasible_at_level(q, s, time_limit=45)
        if not r.feasible:
            print(f"q={q} s={s}: MILP status {r.status} (skipped, not a claim)")
            continue
        ok, c = independent_check(r.g, q, s)
        neg = sorted({bin(T).count('1') for T in c if c[T] < 0})
        print(f"q={q:3d} s={s:3d} MILP witness: independent={ok}  levels={neg}  "
              f"c_[q]={c[(1 << q) - 1]}")

    print("\nheuristic on the q >= 3s boundary (small budget):")
    for (q, s) in [(6, 2), (9, 3), (12, 4)]:
        t0 = time.time()
        ws = greedy_warm_start(q, s)
        h = heuristic_search(q, s, restarts=3, iters=6000, seed=1, init=ws)
        if h.feasible:
            ok, c = independent_check(h.g, q, s)
            print(f"  q={q:3d} s={s:3d} q/s={q/s:.3f}: FEASIBLE, independent={ok}, "
                  f"levels={sorted({bin(T).count('1') for T in c if c[T] < 0})} "
                  f"[{time.time()-t0:.1f}s]")
        else:
            print(f"  q={q:3d} s={s:3d} q/s={q/s:.3f}: not found (best J={h.best_J}, "
                  f"floor 4 = only c_[q]>=1 unsatisfied) [{time.time()-t0:.1f}s]")