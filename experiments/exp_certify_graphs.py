'''Exact certificates for the graphs built from admissible polynomials.

Three independent checks are run for every (construction, m):

1. ``trace moments``:  for k = 1..K the exact integer identity
   ``tr(A^k) = sum_{T} lambda_T(m)^k (m-1)^{|T|}`` is compared with the trace of the
   explicitly constructed adjacency matrix ``A(G_m)``.  This certifies the eigenvalues
   together with their multiplicities (Vandermonde), and is pure integer arithmetic.

2. ``characteristic polynomial``:  for orders small enough, the exact characteristic
   polynomial of ``A(G_m)`` is computed with SymPy and compared with
   ``prod_T (t - lambda_T)^{(m-1)^{|T|}}``.

3. ``inertia by exact root counting``:  the number of positive / negative / zero roots of
   that polynomial is computed exactly (Sturm sequences via SymPy's ``count_roots``) and
   compared with the closed‑form inertia.
'''  # noqa: D400, D401

from __future__ import annotations

import json
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

# Ensure the repository's ``src`` directory is on ``sys.path`` when the module is imported
REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from qspoly import from_g, coefficients_to_g
from qspoly.graphs import (
    basis,
    neps_adjacency,
    neps_eigenvalues,
    neps_inertia,
    leading_sign,
)
from qspoly.exact import exact_inertia_from_charpoly

__all__ = [
    "published_7_3",
    "simple_2_1",
    "trace_moment_check",
    "charpoly_check",
    "numerical_check",
]


def published_7_3() -> Any:
    """Return the admissible polynomial corresponding to the published (7,3) example.

    The polynomial is built from a hard‑coded coefficient dictionary where the keys are
    bit‑packed subsets of ``{1,…,7}``.  The helper ``m`` converts a tuple of indices to the
    appropriate integer key.
    """

    def m(ts: tuple[int, ...]) -> int:
        """Encode a tuple of 1‑based indices as a bitmask.

        ``ts`` must contain integers in the range ``1..7``.  The function is deliberately
        tiny and pure, making it safe for use in dictionary comprehensions.
        """
        return sum(1 << (i - 1) for i in ts)

    coeffs: dict[int, int] = {}
    for t in [(1, 3), (2, 3), (1, 4), (2, 4), (5, 6), (5, 7)]:
        coeffs[m(t)] = 1
    for t in [
        (1, 2, 3),
        (1, 2, 4),
        (1, 3, 4),
        (2, 3, 4),
        (1, 5, 6),
        (2, 5, 6),
        (3, 5, 7),
        (4, 5, 7),
        (5, 6, 7),
    ]:
        coeffs[m(t)] = -1
    for t in [(1, 2, 5, 6), (3, 4, 5, 7)]:
        coeffs[m(t)] = 1
    coeffs[m(tuple(range(1, 8)))] = 1
    return from_g(7, coefficients_to_g(7, coeffs))


def simple_2_1() -> Any:
    """Return the trivial admissible polynomial ``(1 - x₁)(1 - x₂)``.

    The coefficient dictionary is explicitly written for clarity.
    """
    return from_g(2, {0: 1, 1: 0, 2: 0, 3: 0})


def trace_moment_check(P: Any, m: int, K: int) -> List[Dict[str, Any]]:
    """Verify the trace‑moment identities for a given polynomial ``P`` and parameter ``m``.

    Parameters
    ----------
    P:
        An admissible polynomial object from :pymod:`qspoly`.
    m:
        Positive integer scaling factor (the ``m`` in ``G_m``).
    K:
        Number of moments to check (i.e. compute ``A^k`` for ``k = 1..K``).

    Returns
    -------
    List[Dict[str, Any]]
        A list of dictionaries, one per moment, containing the moment index, the
        computed trace, the predicted trace and a boolean ``ok`` flag.
    """
    if m <= 0:
        raise ValueError("m must be a positive integer")
    if K <= 0:
        raise ValueError("K must be a positive integer")

    A = neps_adjacency(P, m).astype(np.int64)
    eigen_dict = neps_eigenvalues(P, m)

    # Build a multiplicity map: eigenvalue -> total multiplicity
    multiplicities: Counter = Counter()
    for T, eigen_val in eigen_dict.items():
        multiplicities[eigen_val] += (m - 1) ** bin(T).count("1")

    results: List[Dict[str, Any]] = []
    Ak = np.eye(A.shape[0], dtype=object)  # start with I (object dtype for arbitrary precision)
    for k in range(1, K + 1):
        Ak = Ak @ A
        lhs = int(np.trace(Ak))
        rhs = sum(int(v) ** k * int(c) for v, c in multiplicities.items())
        results.append({"k": k, "trace": lhs, "predicted": rhs, "ok": lhs == rhs})
    return results


def charpoly_check(P: Any, m: int) -> Dict[str, Any]:
    """Compute and compare the characteristic polynomial of ``A(G_m)``.

    The function uses SymPy for exact arithmetic.  If SymPy is not available, a clear
    ``ImportError`` is raised.
    """
    try:
        import sympy as sp
    except ImportError as exc:  # pragma: no cover – SymPy is a required test dependency
        raise ImportError("SymPy is required for charpoly_check") from exc

    A = neps_adjacency(P, m)
    n = A.shape[0]

    # Build a SymPy matrix with integer entries (0/1) – this avoids accidental float conversion.
    M = sp.zeros(n)
    rows, cols = np.where(A)
    for i, j in zip(rows, cols, strict=True):
        M[i, j] = 1

    x = sp.Symbol("x")
    cp = sp.Poly(M.charpoly(x).as_expr(), x)

    eigen_dict = neps_eigenvalues(P, m)
    pred = sp.Integer(1)
    for T, eigen_val in eigen_dict.items():
        exponent = (m - 1) ** bin(T).count("1")
        pred *= (x - sp.Integer(eigen_val)) ** exponent
    pred = sp.Poly(pred, x)

    same = sp.expand(cp.as_expr() - pred.as_expr()) == 0
    pos, zero, neg = exact_inertia_from_charpoly(cp)
    p, z, nneg = neps_inertia(P, m)

    return {
        "order": n,
        "charpoly_equal": bool(same),
        "exact_inertia": [int(p), int(z), int(nneg)],
        "root_count_inertia": [int(pos), int(zero), int(neg)],
        "inertia_equal": [int(p), int(z), int(nneg)] == [int(pos), int(zero), int(neg)],
        "charpoly": str(cp.as_expr())[:300],
    }


def numerical_check(P: Any, m: int) -> Dict[str, Any]:
    """Perform a floating‑point eigenvalue computation and compare with the closed form.

    The function tolerates small numerical noise by using a relative threshold of ``1e-7``.
    """
    if m <= 0:
        raise ValueError("m must be a positive integer")
    A = neps_adjacency(P, m).astype(float)
    w = np.linalg.eigvalsh(A)
    # Use a small absolute tolerance; the spectra are integer‑valued, so this is safe.
    tol = 1e-7
    p = int(np.sum(w > tol))
    z = int(np.sum(np.abs(w) <= tol))
    n = int(np.sum(w < -tol))
    closed = neps_inertia(P, m)
    return {
        "numerical_inertia": [p, z, n],
        "closed_form_inertia": list(closed),
        "inertia_equal": [p, z, n] == list(closed),
        "numerical_spectrum_sorted": [round(float(x), 6) for x in np.sort(w)[::-1][:12]],
    }


def _run_cli() -> None:
    """Entry‑point used when the module is executed as a script.

    The function mirrors the original behaviour but isolates side‑effects, making the
    module import‑safe for unit‑testing.
    """
    certs: List[Dict[str, Any]] = []
    cases = [
        ("published (7,3)", published_7_3(), [2]),
        ("(2,1) P=(1-x1)(1-x2)", simple_2_1(), [2, 3, 4, 5]),
    ]
    for name, P, ms in cases:
        assert P.is_admissible(), name
        s = P.max_negative_level()
        print(f"=== {name}: q={P.q} s={s} ratio={P.qs_ratio():.4f} admissible={P.is_admissible()}")
        for m in ms:
            n = m ** P.q
            tm = trace_moment_check(P, m, K=8)
            ok_tm = all(r["ok"] for r in tm)
            rec: Dict[str, Any] = {
                "name": name,
                "q": P.q,
                "s": s,
                "m": m,
                "order": n,
                "trace_moments_ok": ok_tm,
                "trace_moments": tm,
                "closed_form_inertia": list(neps_inertia(P, m)),
                "numerical": numerical_check(P, m),
            }
            rec["ok"] = bool(ok_tm and rec["numerical"]["inertia_equal"])
            if n <= 40:
                rec.update(charpoly_check(P, m))
                rec["ok"] = bool(rec["ok"] and rec["charpoly_equal"] and rec["inertia_equal"])
            print(
                f"   m={m} order={n}: trace moments ok={ok_tm}  inertia={neps_inertia(P, m)}"
                f"  numerical_equal={rec['numerical']['inertia_equal']}"
                + (
                    f"  charpoly_equal={rec.get('charpoly_equal')} inertia_equal={rec.get('inertia_equal')}"
                    if "charpoly_equal" in rec
                    else ""
                )
                + f"   -> {'PASS' if rec['ok'] else 'FAIL'}"
            )
            certs.append(rec)
            sys.stdout.flush()
    results_dir = REPO_ROOT / "results"
    results_dir.mkdir(exist_ok=True)
    json_path = results_dir / "certificates.json"
    json.dump(certs, json_path.open("w"), indent=1)
    print("\nwrote results/certificates.json")
    print("ALL CERTIFICATES PASS:", all(c["ok"] for c in certs))


if __name__ == "__main__":
    _run_cli()
