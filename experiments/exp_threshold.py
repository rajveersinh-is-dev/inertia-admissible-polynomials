"""Boundary search: is there an admissible (q,s) with q >= 3s?

For each s we test the critical level q = 3s (and a few above) with
  * exact MILP when the model fits in memory (2^q <= 2^19), else
  * the heuristic local search,
and record whether a witness was found.  A "not found" verdict is evidence only.
Budgets are wall-clock bounded so this script always terminates quickly.
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))
from qspoly import from_g
from qspoly.milp import feasible_at_level
from qspoly.heuristic import heuristic_search, greedy_warm_start, verify

TL = float(sys.argv[1]) if len(sys.argv) > 1 else 30.0     # MILP seconds per probe
SMIN = int(sys.argv[2]) if len(sys.argv) > 2 else 1
SMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 6
BUDGET = float(sys.argv[4]) if len(sys.argv) > 4 else 12.0  # heuristic seconds per (q,s)
MILP_MAX_Q = 19

rows = []
t0 = time.time()
print("searching for admissible (q,s) with q >= 3s (ratio >= 3)")
for s in range(SMIN, SMAX + 1):
    for q in range(max(2, 3 * s), 3 * s + 5):
        t1 = time.time()
        witness = None
        method = None
        if q <= MILP_MAX_Q:
            r = feasible_at_level(q, s, time_limit=TL)
            if r.feasible:
                witness, method = r.g, "milp-exact"
            elif r.status == "infeasible":
                rows.append({"q": q, "s": s, "ratio": q / s, "result": "infeasible (exact)"})
                print(f"  q={q:3d} s={s:3d} q/s={q/s:.4f}: infeasible (EXACT) "
                      f"[{time.time()-t1:.1f}s]", flush=True)
                break
        if witness is None:
            method = "heuristic"
            ws = greedy_warm_start(q, s)
            h = heuristic_search(q, s, restarts=2, iters=4000, seed=q, init=ws,
                                 time_budget=BUDGET)
            if h.feasible and verify(h.g, q, s):
                witness = h.g
        if witness is not None:
            P = from_g(q, witness)
            rows.append({"q": q, "s": s, "ratio": q / s, "result": "FEASIBLE",
                         "method": method, "levels": P.negative_levels(),
                         "g": {str(k): v for k, v in sorted(witness.items())}})
            print(f"  q={q:3d} s={s:3d} q/s={q/s:.4f}: *** FEASIBLE via {method} "
                  f"levels={P.negative_levels()} [{time.time()-t1:.1f}s]", flush=True)
            break
        rows.append({"q": q, "s": s, "ratio": q / s,
                     "result": "not found (heuristic only; not a proof)"})
        print(f"  q={q:3d} s={s:3d} q/s={q/s:.4f}: not found by heuristic "
              f"[{time.time()-t1:.1f}s]", flush=True)
    json.dump(rows, open('results/threshold_scan.json', 'w'), indent=1)
feasible_rows = [r for r in rows if r['result'] == 'FEASIBLE']
print(f"\ntotal {time.time()-t0:.0f}s; FEASIBLE with q>=3s: "
      f"{[(r['q'], r['s']) for r in feasible_rows] or 'none'}")