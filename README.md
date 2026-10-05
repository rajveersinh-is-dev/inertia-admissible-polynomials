# Inertia Admissible Polynomials

**Boolean admissibility, the $q/s$ ceiling for graph positive inertia, and a conjectured barrier at $3$.**

For a graph $G$ let $\operatorname{In}(G) = (n^+(G), n^0(G), n^-(G))$ be the numbers of
positive, zero and negative eigenvalues of $A(G)$. Elphick, Kumar, Pragada and Spier
(arXiv:2609.06319) show that a **$(q,s)$-admissible Boolean polynomial** — a multilinear
$P\in\mathbb{Z}[x_1,\dots,x_q]$ whose Möbius coefficients $c_T$ are non-negative above
level $s$ with $c_{[q]}>0$ — produces graphs with
$$n^+ = \Theta\big((n^-)^{q/s}\big).$$
So the largest attainable exponent $q/s$ decides whether $n^+$ admits *any* polynomial
bound in $n^-$, which is an open question in that paper (their Problem 4.1).

This repository gives an **exact linearisation** of admissibility, **solves the $s=1$ case
completely**, produces **certified exact data for $q\le 7$**, and states the sharp
conjecture that the exponent never reaches $3$.

---

## Research question

> **Main question.** Which pairs $(q,s)$ are $(q,s)$-admissible, and what is the largest
> attainable ratio $q/s$?
>
> **Conjecture (barrier).** Every admissible pair satisfies $q<3s$, so
> $n^+(G)=O((n^-(G))^{3})$ uniformly over graphs.

The conjecture is **not proved** here. It is stated because the certified data and a
targeted search are consistent with it, because it is sharp in spirit ($7/3$ is attained),
and because proving it would settle the currently open question above.

## Main result

**Theorem (`s-one`).** An admissible polynomial has $s=1$ only when $q=2$; in that case it
is the unique polynomial $P=(1-x_1)(1-x_2)$. Consequently no exponent $>2$ can arise from
level $1$. *(Complete proof by case analysis; verified exhaustively for $q\le4$.)*

**Theorem (`certified`).** The minimal admissible level is exactly
$$s_{\min}(2),\dots,s_{\min}(7) = 1,2,2,3,3,3,$$
so the best ratio in that range is $7/3$, and in every case $q<3s_{\min}(q)$. Each entry
rests on a MILP *infeasibility certificate* plus a from-scratch–verified witness.

**Proposition (`ref`).** Admissibility is equivalent to feasibility of a
$2^q-1$-variable $0$–$1$ linear program over the Möbius transform of $g(S)=P(1_S)$:
$$c_{[q]}=\sum_{S\subseteq[q]}(-1)^{q-|S|}g(S)\ge 1,\qquad
c_T=\sum_{S\subseteq T}(-1)^{|T|-|S|}g(S)\ge 0\ \ (|T|>s),\qquad g\in\{0,1\},\ g([q])=0 .$$
Level feasibility is monotone in $s$ and in $q$, which is what licenses binary search and
makes negative answers checkable.

## Why this is interesting

* It is a **clean reformulation**, not a heuristic: the whole search reduces to the signs of
  an alternating sum over sub-subsets of $\{0,1\}^q$.
* It turns an **open question about eigenvalues of graphs** into a **combinatorial
  optimisation problem on $2^q$ binary variables**, which is exactly the right abstraction
  for both branch-and-bound and local search.
* The barrier $q<3s$ would give, for the first time, a *polynomial* bound on the positive
  inertia index in terms of the negative one — the only known bound there is exponential.

## Main theorem, stated plainly

For a $(q,s)$-admissible $P$ with coefficients $c_T$ the constructed graphs $G_m$ (order
$m^q$) have, for all sufficiently large $m$,
$$n^+(G_m)=(m-1)^q+O(m^{q-1}),\qquad n^-(G_m)=\Theta(m^s),\qquad n^0(G_m)=0,$$
with every eigenvalue an *integer*
$$\lambda_T(m)=\sum_{T\subseteq R\subseteq[q]}c_R\,m^{\,q-|R|}\qquad\text{with multiplicity }(m-1)^{|T|}.$$
Exact inertia for $m=2$: $(104,0,24)$ on $128$ vertices.

## Computation and verification

Every certificate is exact; floating point is used only as a non-authoritative cross-check.

| check | what it proves | scope |
|---|---|---|
| trace moments $\operatorname{tr}A^k=\sum_T\lambda_T^k(m-1)^{\|T\|}$ | spectrum **with multiplicities** (Vandermonde) | all constructions, $k\le 8$ |
| characteristic polynomial $\prod_T(t-\lambda_T)^{(m-1)^{\|T\|}}=\det(tI-A)$ | closed-form spectrum | order $\le 40$ |
| exact inertia from exact algebraic roots | signs of all eigenvalues | order $\le 40$ |
| `verify()` from-scratch re-check | every heuristic/MILP witness is genuine | all $q$ |

Reproduce: `python experiments/run_all.py` (fast, ~2 min) or `--full` (adds the certified
MILP scan).

## Repository structure

```
src/qspoly/
  boolean.py     multilinear Boolean polynomials, zeta/Moebius, admissibility
  graphs.py      NEPS construction, closed-form spectrum, exact inertia
  exact.py       exact inertia of an integer symmetric matrix (no floating point)
  milp.py        exact 0-1 feasibility (HiGHS), certified/uncertified verdicts
  heuristic.py   local search for large q, with mandatory from-scratch verification
  theorems.py    machine-checkable shadows of the proved statements
experiments/     run_all.py plus the individual experiments
tests/           19 tests; run with `pytest tests/` or `python tests/test_qspoly.py`
paper/main.tex   full paper with proofs (compiles with tectonic)
figures/         generated figures
results/         machine-generated data and certificates
docs/            methodology and mathematical notes
```

## Reproducing results

```bash
pip install -e .            # numpy, scipy, sympy, matplotlib
python experiments/run_all.py          # certificates + figures (fast)
python experiments/run_all.py --full   # + certified MILP scan (slow)
pytest tests/                          # 19 tests
tectonic -X compile paper/main.tex     # builds paper/main.pdf
```

## Honest status and limitations

* **Proved:** the linearisation (Prop. `ref`), level $0$ impossible, monotonicity, the
  complete $s=1$ classification, the exact spectrum/inertia theorem, the certified table
  for $q\le 7$.
* **Not proved:** the barrier conjecture $q<3s$.
* **Not improved:** we found no admissible pair beating the known exponent $7/3$ (or the
  $8/3$ that Elphick–Kumar–Pragada–Spier report without exhibiting).
* **Certified only to $q\le 7$.** Runs up to $q=12$ completed but are *uncertified*,
  because MILP probes hit time limits; a time-limited probe is never recorded as
  "infeasible". Uncertified values are kept out of the paper's claims.
* A tempting generalisation (Kneser graphs on $r$-subsets instead of $2$-subsets) **fails**,
  and the paper explains the algebraic trap: the identity
  $K=J_N+(r-1)I_N-C^{\mathsf T}C$ holds only for $r=2$.

## What we tried that failed (negative results)

1. **Exponential construction via Kneser graphs.** Predicted $\pp=\binom{k}{r}$, i.e.
   $\approx 2^k$. Fails because the $b$-vertices form an independent set of that size, and
   $\pp(G)\le n-\alpha(G)$ then caps it at $k$. Documented in §7 of the paper.
2. **Local search at the boundary $q=3s$.** Drives the violation count to its floor $4$,
   i.e. finds $g$ with all top-level $c_T\ge0$ but $c_{[q]}=0$ rather than $\ge1$. Making
   the last condition hold is the real obstacle, and it is where we would look next.

## Citation

If you use this, please cite the framework you build on:

```bibtex
@article{Akbari2026,
  author  = {Akbari, Saieed and Elphick, Clive and Kumar, Hitesh and Pragada, Shivaramakrishna and Tang, Quanyu},
  title   = {A new conjecture on the inertia of graphs},
  journal = {Discrete Mathematics}, volume = {349}, year = {2026}, pages = {114953},
  note    = {arXiv:2508.01163}
}
@article{Elphick2026,
  author  = {Elphick, Clive and Kumar, Hitesh and Pragada, Shivaramakrishna and Spier, Tom\'as Jung},
  title   = {Resolution of a problem of Mohar on non-positive inertia},
  year    = {2026}, note = {arXiv:2609.06319}
}
```

Full bibliography in `paper/main.tex`.

## License

MIT — see `LICENSE`.