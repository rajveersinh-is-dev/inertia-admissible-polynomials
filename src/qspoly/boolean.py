r"""Boolean multilinear polynomials and $(q,s)$-admissibility.

Definitions
-----------
Let $q\ge 1$.  A *multilinear polynomial* in $q$ variables is
$$P(x_1,\dots,x_q)=\sum_{T\subseteq[q]}c_T\prod_{i\in T}x_i,\qquad c_T\in\mathbb Z,$$
where $[q]=\{1,\dots,q\}$.  Write $1_S$ for the indicator vector of $S\subseteq[q]$, so
$\prod_{i\in T}1_S(i)=1$ iff $T\subseteq S$.

$P$ is a **Boolean polynomial** if $P(x)\in\{0,1\}$ for all $x\in\{0,1\}^q$, equivalently
if $g(S):=P(1_S)\in\{0,1\}$ for all $S\subseteq[q]$.

The $g(S)$ and the $c_T$ are related by the (Boolean-lattice) zeta/Mobius pair
$$g(S)=\sum_{T\subseteq S}c_T,\qquad c_T=\sum_{S\subseteq T}(-1)^{|T|-|S|}g(S).$$

$P$ is **$(q,s)$-admissible** (following Elphick--Kumar--Pragada--Spier) if it is Boolean and
$$P(1,\dots,1)=0\;(=g([q])),\qquad c_{[q]}>0,\qquad
s=\max\{\,|T|:\ c_T<0\,\}.$$

Such a $P$ produces, for all sufficiently large $m$, a graph $G_m$ of order $m^q$ with
$n^+(G_m)=\Theta(m^q)$, $n^-(G_m)=\Theta(m^s)$ and $n^0(G_m)=0$, hence
$n^+(G_m)=\Theta\bigl(n^-(G_m)^{q/s}\bigr)$.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Dict, Iterable, List, Sequence, Tuple

__all__ = [
    "full",
    "popcount",
    "subsets_of",
    "MobiusTransform",
    "AdmissiblePolynomial",
    "from_g",
    "coefficients_to_g",
    "negative_levels",
    "is_admissible",
    "boolean_p",
    "qs_ratio",
]

full = lambda q: (1 << q) - 1  # noqa: E731  (bitmask of [q])


@lru_cache(maxsize=None)
def popcount(mask: int) -> int:
    """Number of elements of a subset given by its bitmask."""
    return bin(mask).count("1")


def subsets_of(mask: int) -> Iterable[int]:
    """Yield all subsets (as bitmasks) of the subset given by ``mask``."""
    sub = 0
    while True:
        yield sub
        if sub == mask:
            return
        sub = (sub - mask) & mask


# --------------------------------------------------------------------------------------
# zeta / Mobius transforms on the Boolean lattice
# --------------------------------------------------------------------------------------
def mobius(q: int, g: Dict[int, int]) -> Dict[int, int]:
    """Mobius transform: ``c_T = sum_{S subset T} (-1)^{|T|-|S|} g(S)``."""
    c: Dict[int, int] = {}
    for T in range(1 << q):
        t = popcount(T)
        c[T] = sum((-1) ** (t - popcount(S)) * g[S] for S in subsets_of(T))
    return c


def zeta(q: int, c: Dict[int, int]) -> Dict[int, int]:
    """Inverse of :func:`mobius`: ``g(S) = sum_{T subset S} c_T``."""
    g: Dict[int, int] = {}
    for S in range(1 << q):
        g[S] = sum(c.get(T, 0) for T in subsets_of(S))
    return g


def from_g(q: int, g: Dict[int, int]) -> "AdmissiblePolynomial":
    """The unique multilinear Boolean polynomial with ``P(1_S) = g(S)``."""
    return AdmissiblePolynomial(q=q, g=g)


def coefficients_to_g(q: int, c: Dict[int, int]) -> Dict[int, int]:
    """``g(S) = sum_{T subset S} c_T``  (assumes ``c`` is the coefficient dictionary)."""
    return zeta(q, c)


# --------------------------------------------------------------------------------------
# the polynomial object
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class AdmissiblePolynomial:
    """A multilinear polynomial stored by its values ``g(S) = P(1_S)`` on the cube.

    Attributes
    ----------
    q : number of variables
    g : ``g[S] = P(1_S)`` for ``0 <= S < 2**q``
    """

    q: int
    g: Dict[int, int]

    def __post_init__(self) -> None:
        if len(self.g) != 2 ** self.q:
            raise ValueError("g must have 2**q entries")
        if set(self.g.values()) - {0, 1}:
            raise ValueError("g must be {0,1}-valued")

    # -- basic data ---------------------------------------------------------------
    @property
    def c(self) -> Dict[int, int]:
        """Moeius transform of ``g``: the coefficients ``c_T`` of ``prod_{i in T} x_i``."""
        return mobius(self.q, self.g)

    def evaluate(self, x: Sequence[int]) -> int:
        """``P(x_1, ..., x_q)``."""
        return self.g[sum(1 << i for i, xi in enumerate(x) if xi)]

    def is_boolean(self) -> bool:
        return set(self.g.values()) <= {0, 1}

    def coefficients(self) -> Dict[int, int]:
        return {T: v for T, v in self.c.items() if v}

    # -- admissibility -----------------------------------------------------------
    def negative_levels(self) -> List[int]:
        """The (sorted) set of ``|T|`` with ``c_T < 0``."""
        return sorted({popcount(T) for T, v in self.c.items() if v < 0})

    def max_negative_level(self) -> int:
        """``s = max{|T| : c_T < 0}``;  ``-1`` if all coefficients are nonnegative."""
        lv = self.negative_levels()
        return lv[-1] if lv else -1

    def is_admissible(self, s: int | None = None) -> bool:
        """Check the three defining conditions; ``s`` is recomputed if omitted."""
        if not self.is_boolean():
            return False
        c = self.c
        if self.g[full(self.q)] != 0:  # P(1,...,1) = 0
            return False
        if c[full(self.q)] <= 0:
            return False
        smax = self.max_negative_level()
        if smax < 0:
            return False
        if all(c[T] >= 0 for T in c if popcount(T) > smax):
            return True
        return False

    def qs_ratio(self) -> float:
        s = self.max_negative_level()
        if s <= 0:
            return float("inf")
        return self.q / s

    def support(self) -> List[Tuple[int, ...]]:
        """The ``T`` with ``c_T != 0``, as tuples of indices (1-based)."""
        return [tuple(i + 1 for i in range(self.q) if (T >> i) & 1)
                for T in sorted(self.c) if self.c[T]]

    def polynomial_expression(self) -> str:
        """A readable multilinear expansion with ``x_i`` symbols."""
        terms = []
        for T, v in sorted(self.coefficients().items()):
            mono = "*".join(f"x{i+1}" for i in range(self.q) if (T >> i) & 1) or "1"
            terms.append(f"{v:+d}*{mono}" if terms else (f"{v}*{mono}" if v != 1 else mono))
        return " ".join(terms)

    def count_table(self) -> Dict[Tuple[int, int], int]:
        """``{(level, sign): count}`` for the coefficients; sign in {-1,0,1}."""
        table: Dict[Tuple[int, int], int] = {}
        for T, v in self.c.items():
            sgn = (v > 0) - (v < 0)
            key = (popcount(T), sgn)
            table[key] = table.get(key, 0) + 1
        return table


# --------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------
def is_admissible(g: Dict[int, int], q: int) -> bool:
    return from_g(q, g).is_admissible()


def boolean_p(P: AdmissiblePolynomial, S: int) -> int:
    """``P(1_S)``."""
    return P.g[S]


def qs_ratio(g: Dict[int, int], q: int) -> float:
    return from_g(q, g).qs_ratio()