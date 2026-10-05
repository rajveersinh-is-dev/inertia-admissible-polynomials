"""Exhaustive verification of Theorem s-one (s=1 forces q=2) over all 2^(2^q) polynomials.

For each q we enumerate every function g : 2^[q] -> {0,1} with g([q]) = 0 and test
admissibility.  This is feasible up to q = 4 (2^15 = 32768 candidates).
"""
import sys, os, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))
from qspoly.boolean import mobius, popcount, full

def admissible(g, q, c):
    fq = full(q)
    if g[fq] != 0: return None
    if c[fq] < 1: return None
    neg = [popcount(T) for T in c if c[T] < 0]
    if not neg: return None
    s = max(neg)
    for T in c:
        if T != fq and popcount(T) > s and c[T] < 0:
            return None
    return s

if __name__ == '__main__':
    QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    for q in range(1, QMAX + 1):
        fq = full(q)
        free = [S for S in range(1 << q) if S != fq]
        found = {}
        total = 0
        for bits in range(1 << len(free)):
            g = {fq: 0}
            for i, S in enumerate(free):
                g[S] = (bits >> i) & 1
            total += 1
            c = mobius(q, g)
            s = admissible(g, q, c)
            if s is not None:
                found.setdefault(s, []).append(g)
        print(f"q={q}: examined {total} polynomials g with g([q])=0")
        for s in sorted(found):
            print(f"    admissible with s={s}: {len(found[s])} polynomial(s)")
        if 1 in found:
            print(f"    *** s=1 admissible at q={q} with |q|={q}")
        sys.stdout.flush()