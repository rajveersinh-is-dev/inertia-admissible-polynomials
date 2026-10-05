"""Reproduce every numerical result reported in the paper and README.

    python experiments/run_all.py            # fast tier (~2 min)
    python experiments/run_all.py --full     # + exact MILP certification scan (~30 min)

Deterministic: all random seeds are fixed and every accepted witness is re-verified by
recomputing the Moebius transform from scratch (``qspoly.heuristic.verify``).
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXPS = os.path.join(ROOT, "experiments")
TESTS = os.path.join(ROOT, "tests")


def run(script: str, *args: str) -> None:
    path = script if os.path.isabs(script) else None
    if path is None:
        for base in (EXPS, TESTS):
            cand = os.path.join(base, script)
            if os.path.exists(cand):
                path = cand
                break
    assert path is not None, f"script not found: {script}"
    print(f"\n=== {script} {' '.join(args)}", flush=True)
    t0 = time.time()
    r = subprocess.run([sys.executable, path, *args], cwd=ROOT)
    if r.returncode != 0:
        raise SystemExit(f"{script} failed with code {r.returncode}")
    print(f"--- {script} finished in {time.time()-t0:.1f}s", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true",
                    help="also run the certified MILP scan (slow)")
    ap.add_argument("--qmax", type=int, default=7,
                    help="largest q for the certified scan in --full mode")
    args = ap.parse_args()

    run("test_qspoly.py")                    # the real test suite
    run("exp_certify_graphs.py")             # exact certificates
    run("test_s1_theorem.py", "4")           # exhaustive: s=1 forces q=2
    run("test_witness_check.py")             # independent re-verification of witnesses
    run("exp_threshold.py", "15", "1", "3", "5")  # boundary q >= 3s for s = 1,2,3
    if args.full:
        run("exp_scan_certified.py", "2", str(args.qmax), "900")
    run("make_figures.py")
    print("\nAll experiments completed.")


if __name__ == "__main__":
    main()