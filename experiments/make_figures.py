"""Figures: admissible region, inertia growth, and the coefficient profile."""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from qspoly import from_g, coefficients_to_g
from qspoly.graphs import neps_inertia

FIG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'figures')
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3,
                     "figure.dpi": 140, "savefig.bbox": "tight"})


def published_7_3():
    def m(ts): return sum(1 << (i - 1) for i in ts)
    c = {}
    for t in [(1,3),(2,3),(1,4),(2,4),(5,6),(5,7)]:      c[m(t)] = 1
    for t in [(1,2,3),(1,2,4),(1,3,4),(2,3,4),
              (1,5,6),(2,5,6),(3,5,7),(4,5,7),(5,6,7)]:   c[m(t)] = -1
    for t in [(1,2,5,6),(3,4,5,7)]:                      c[m(t)] = 1
    c[m(tuple(range(1, 8)))] = 1
    return from_g(7, coefficients_to_g(7, c))


def load_scan():
    p = 'results/scan_certified.json'
    if not os.path.exists(p):
        return []
    return json.load(open(p))


# ---------------------------------------------------------------- Fig 1: s_min(q)
def fig1():
    rows = [r for r in load_scan() if r.get("certified") and r.get("s")]
    if not rows:
        return
    qs = [r["q"] for r in rows]; ss = [r["s"] for r in rows]
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
    ax[0].plot(qs, ss, "o-", color="#2b6cb0", lw=1.6, ms=5, label=r"$s_{\min}(q)$ (exact)")
    ax[0].plot(qs, [q / 3 for q in qs], "--", color="#c53030", lw=1.4,
               label=r"barrier $q=3s$")
    ax[0].set_xlabel("$q$ (number of variables)"); ax[0].set_ylabel("$s$ (level)")
    ax[0].set_title("Minimal admissible level"); ax[0].legend(frameon=False, fontsize=7.5)
    ax[0].set_ylim(0, max(ss) + 1)

    rat = [q / s for q, s in zip(qs, ss)]
    ax[1].bar(range(len(qs)), rat, color="#2b6cb0", alpha=0.85)
    ax[1].axhline(3, ls="--", color="#c53030", lw=1.4, label="conjectured ceiling $3$")
    ax[1].axhline(7/3, ls=":", color="#2f855a", lw=1.6, label=r"best known ratio $7/3$")
    ax[1].axhline(8/3, ls=":", color="#805ad5", lw=1.6,
                  label=r"unpublished $32/12=8/3$")
    ax[1].set_xticks(range(len(qs))); ax[1].set_xticklabels([str(q) for q in qs])
    ax[1].set_xlabel("$q$"); ax[1].set_ylabel(r"$q/s_{\min}(q)$")
    ax[1].set_title("Achievable exponent $n^+ \\sim (n^-)^{q/s}$")
    ax[1].legend(frameon=False, fontsize=7)
    fig.suptitle(r"Certified minimal levels $s_{\min}(q)$ for $(q,s)$-admissible polynomials",
                 y=1.04, fontsize=10)
    fig.savefig(f"{FIG}/fig1_admissible_region.png")
    plt.close(fig)


# ------------------------------------------------- Fig 2: inertia growth (7,3)
def fig2():
    P = published_7_3()
    ms = [2, 3, 4, 5, 6, 7, 8, 9, 10]
    pos = [neps_inertia(P, m)[0] for m in ms]
    neg = [neps_inertia(P, m)[2] for m in ms]
    fig, ax = plt.subplots(figsize=(4.4, 3.4))
    ax.plot(ms, pos, "o-", color="#2b6cb0", label=r"$n^+(G_m)$")
    ax.plot(ms, neg, "s-", color="#c53030", label=r"$n^-(G_m)$")
    ref = np.array(ms, dtype=float) ** 7 * 0.55
    ax.plot(ms, ref, ":", color="#4a5568", lw=1.2, label=r"$\propto m^{7}$")
    ax.plot(ms, np.array(ms, dtype=float) ** 3 * 0.55, ":", color="#4a5568", lw=1.2,
            label=r"$\propto m^{3}$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("$m$"); ax.set_ylabel("eigenvalue count")
    ax.set_title(r"$(7,3)$-admissible construction: $n^+=\Theta(m^7)$, $n^-=\Theta(m^3)$",
                 fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(f"{FIG}/fig2_inertia_growth.png")
    plt.close(fig)


# ------------------------------------------- Fig 3: coefficient level profile
def fig3():
    P = published_7_3()
    c = P.c
    levels = list(range(0, P.q + 1))
    pos = [sum(1 for T in c if bin(T).count('1') == L and c[T] > 0) for L in levels]
    neg = [sum(1 for T in c if bin(T).count('1') == L and c[T] < 0) for L in levels]
    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    w = 0.4
    ax.bar([L - w/2 for L in levels], pos, width=w, color="#2b6cb0", label=r"$c_T>0$")
    ax.bar([L + w/2 for L in levels], neg, width=w, color="#c53030", label=r"$c_T<0$")
    s = P.max_negative_level()
    ax.axvspan(s + 0.5, P.q + 0.7, color="#f6ad55", alpha=0.25)
    ax.text(s + 0.7, ax.get_ylim()[1] * 0.85, r"forced $c_T\geq 0$", fontsize=7.5,
            color="#744210")
    ax.set_xlabel("level $|T|$"); ax.set_ylabel("number of coefficients $c_T$")
    ax.set_title(r"Möbius coefficients of the $(7,3)$ polynomial ($s=3$)", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(f"{FIG}/fig3_coefficient_profile.png")
    plt.close(fig)


if __name__ == '__main__':
    fig1(); fig2(); fig3()
    print("figures written to", FIG)
    for f in sorted(os.listdir(FIG)):
        print("  ", f)