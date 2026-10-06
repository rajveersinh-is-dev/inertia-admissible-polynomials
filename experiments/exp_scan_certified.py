"""Definitive certified scan of the minimal admissible level s for q = 2..14."""
import sys, os, json, time

from qspoly import from_g
from qspoly.heuristic import verify
from qspoly.milp import minimal_level, feasible_at_level
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))

QMIN, QMAX = int(sys.argv[1]), int(sys.argv[2])
TL = float(sys.argv[3]) if len(sys.argv) > 3 else 900.0
rows = []
t0 = time.time()
for q in range(QMIN, QMAX + 1):
    s, log, cert = minimal_level(q, time_limit=TL)
    rec = {"q": q, "s": s, "ratio": (q / s) if s else None, "certified": cert}
    if s is not None:
        g = feasible_at_level(q, s, time_limit=TL).g
        rec["verified"] = bool(g is not None and verify(g, q, s))
        rec["negative_levels"] = from_g(q, g).negative_levels()
        rec["g"] = {str(k): v for k, v in sorted(g.items())}
    rows.append(rec)
    print(f"q={q:3d}  s_min={s}  ratio={rec['ratio']}  certified={cert}  "
          f"verified={rec.get('verified')}  levels={rec.get('negative_levels')}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    json.dump(rows, open('results/scan_certified.json', 'w'), indent=1)
print("best certified ratio:", max((r['ratio'] for r in rows if r['ratio'] and r['certified']),
                                   default=None))