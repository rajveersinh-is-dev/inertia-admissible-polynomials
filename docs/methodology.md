# Methodology and reproducibility notes

## Determinism

All experiments are deterministic:

* MILP feasibility is solved by HiGHS through `scipy.optimize.milp`.  HiGHS is
  deterministic for a fixed model; we never rely on parallelism.
* Local search uses `random.Random(seed)` with explicit integer seeds passed per
  experiment.
* No experiment reads the network, no wall-clock value enters a result, and no
  result depends on dictionary iteration order (all group-by loops are over
  `range(2**q)`).

## Certification policy (the important part)

Every numerical claim in the paper is one of three kinds, and the repository records which:

1. **Certified** — a mathematical statement.  Either a proof, or an exact computation whose
   truth does not depend on a search succeeding: e.g. `$s_{\min}(6)=2$` is certified because
   HiGHS returned `infeasible` for level 2 and `optimal` for level 1.
2. **Witness** — a feasible `g` together with `qspoly.heuristic.verify(g, q, s)`, which
   *recomputes the Möbius transform from scratch*.  A witness is accepted only if this
   recomputation agrees.
3. **Uncertified** — anything else, including any value obtained when a solver hit its
   time limit.  A time-limited probe is **unknown**, never "infeasible".  Uncertified
   values appear in the JSON output with `"certified": false` and are **not** used in any
   theorem.

This distinction exists because of a concrete bug we hit: an early version of the local
search maintained `c_{[q]}` with an inconsistent sign in its incremental update, and
reported "feasible" witnesses whose recomputed `c` had `max{|T| : c_T < 0} = q`.  The
from-scratch re-check caught it immediately.  Any search-based method for this problem must
be re-verified; that is why `verify()` is not optional in the code path.

## Floating point policy

Floating point appears in exactly two places, and neither is allowed to certify anything:

* `numpy.linalg.eigvalsh`, used to cross-check the closed-form spectrum;
* Heuristic objective values (used only to guide the search).

All certificates use exact arithmetic:

* spectra via the closed-form integer formula `lambda_T(m)`;
* characteristic polynomials via SymPy over $\mathbb{Z}$;
* inertia via exact algebraic roots with exact sign decisions (SymPy refines intervals).

## Complexity and practical limits

The reformulated program of Prop. `ref` has $2^q - q - 1$ binary variables and
$\sum_{j>s}\binom{q}{j}$ rows, with $\sum_{T:|T|>s}2^{|T|}$ nonzeros.  Consequences:

| $q$ | variables | our experience |
|---|---|---|
| $\le 7$ | $\le 119$ | settles in seconds, fully certified |
| $8$–$12$ | $\le 4083$ | often times out mid-search; uncertified |
| $\ge 15$ | $\ge 32768$ | infeasible for this solver; heuristics only |

Binary search costs $O(\log q)$ MILP calls per $q$ because level feasibility is monotone
(Prop. `monotone`).

## How the spectrum formula is checked

For each construction and each $m$:

* **Trace moments.** $\operatorname{tr}A^k$ for $k\le 8$ against
  $\sum_T\lambda_T^k(m-1)^{|T|}$.  Exact integers; by Vandermonde this determines the
  spectrum *with* multiplicities, so it is a complete certificate, not a spot check.
* **Characteristic polynomial.** For order $\le 40$ we compare with `det(tI - A)` computed
  by SymPy.
* **Inertia.** Exact root classification, compared with the closed form.

Reproduce with `python experiments/exp_certify_graphs.py`; results land in
`results/certificates.json`.

## Environment

`pip install -e .` (Python $\ge 3.10$; numpy, scipy, sympy, matplotlib).  Versions used for
the reported numbers: Python 3.14.6, numpy 2.4.6, scipy 1.18.0, sympy 1.14.0,
matplotlib 3.11.1.  The paper compiles with `tectonic -X compile paper/main.tex`.