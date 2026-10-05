## [1.0.0] - 2026-10-05

Initial release.

### Mathematics
- Exact linearisation of $(q,s)$-admissibility as a $2^q-1$-variable $0$--$1$ linear program.
- Proposition: level $0$ is never admissible.
- Proposition: level feasibility is monotone in $s$ and in $q$; binary search is valid.
- Theorem: $s=1$ forces $q=2$, with the unique witness $(1-x_1)(1-x_2)$.
- Theorem: closed-form spectrum of the NEPS graphs built from an admissible polynomial,
  with exact integer eigenvalues and multiplicities, giving $n^+=\Theta(m^q)$,
  $n^-=\Theta(m^s)$, $n^0=0$.
- Theorem: certified values of $s_{\min}(q)$ for $q=2,\dots,7$ (i.e. $1,2,2,3,3,3$).
- Conjecture: every admissible pair satisfies $q<3s$ (a "barrier"), which would give
  $n^+(G)=O((n^-(G))^3)$ and answer the open Problem 4.1 of
  Elphick--Kumar--Pragada--Spier positively.

### Computation
- Exact MILP feasibility backend (HiGHS via `scipy.optimize.milp`) with a
  certified/uncertified distinction; time-limited probes are never recorded as infeasible.
- Local search for large $q$, accepted only after a from-scratch re-verification of the
  Moebius transform.
- Exact inertia of integer symmetric matrices from exact algebraic roots (no floating point
  in any certificate).
- Three independent cross-checks per construction: exact trace moments, exact characteristic
  polynomial, exact inertia.

### Repository
- Paper (`paper/main.tex`) with complete proofs, `README.md`, methodology and mathematical
  notes, 19 tests, generated figures, machine-readable results.