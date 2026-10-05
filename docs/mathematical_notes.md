# Mathematical notes

Working notes for the claims in `paper/main.tex`.  Kept separate from the paper so that
rejected ideas and dead ends stay visible.

## 1. The objects

A **multilinear** $P\in\mathbb{Z}[x_1,\dots,x_q]$ is $P=\sum_{T\subseteq[q]}c_T\prod_{i\in T}x_i$.
Write $1_S$ for the indicator of $S$.  Then

$$g(S):=P(1_S)=\sum_{T\subseteq S}c_T,$$
so $(g,c)$ is a zeta/Moebius pair on the Boolean lattice.  $P$ is **Boolean** iff
$g\in\{0,1\}^{2^{[q]}}$.

$(q,s)$-**admissible** (Elphick–Kumar–Pragada–Spier): Boolean, $g([q])=0$, $c_{[q]}>0$, and
$s=\max\{|T|:c_T<0\}$.

## 2. The linearisation

Substituting $g$ for $c$ turns everything linear:

$$c_{[q]}=\sum_{S\subseteq[q]}(-1)^{q-|S|}g(S)\ \ge 1,$$
$$c_T=\sum_{S\subseteq T}(-1)^{|T|-|S|}g(S)\ \ge 0\quad (|T|>s),\qquad g\in\{0,1\},\ g([q])=0.$$

Every constraint is linear in the binary variables $g(S)$: a $0$–$1$ LP.

Two corollaries used everywhere:

* **Level $0$ is impossible.** $\sum_Tc(T)=g([q])=0$ with all $c_T\ge0$ forces $c\equiv0$.
* **Monotonicity.** $c_T$ depends only on $g|_{2^T}$, so (i) restricting $g$ to $Q\subseteq[q]$
  with $|Q|>s$ gives a witness for $(|Q|,s)$; (ii) raising $s$ only removes constraints.
  Hence binary search in $s$ is valid, and $s_{\min}(q)$ is well defined.

## 3. Proof of the $s=1$ theorem (expanded)

Set $a=c(\emptyset)=g(\emptyset)\in\{0,1\}$, so $c(\{i\})=g(\{i\})-a\in\{-a,1-a\}$.
$c_{[q]}\ge1$ and $\sum_Tc(T)=0$ give
$$1\le c_{[q]}=-(a+\textstyle\sum_i c_i)-\sum_{|T|\ge2}c_T\le(q-1)a,$$
so $a=1$, $c_i\in\{-1,0\}$, and at least two $c_i=-1$.  If $c_i=c_j=0$ then
$c_{\{i,j\}}=g(\{i,j\})-1<0$, contradicting $c_T\ge0$ for $|T|\ge2$; hence at most one $c_i=0$
and at least $q-1$ of them are $-1$.

* $q\ge4$: four indices with $c=-1$ give, with $J\ge6$ and $K\ge0$,
  $c(X)=g(X)-1+4-J-K\le-2<0$.
* $q=3$: $(c_1,c_2,c_3)=(-1,-1,-1)$ gives $c_{[3]}=2-J\le-1$; $(-1,-1,0)$ gives
  $c_{[3]}=1-J\le0$.  Both contradict $c_{[3]}\ge1$.
* $q=2$: forced, and $P=(1-x_1)(1-x_2)$.

## 4. The construction and its spectrum

$\mathcal F=\{S: P(1_{[q]\setminus S})=1\}$; $G_m=\mathrm{NEPS}(K_m,\dots,K_m;\{\mathbf 1_S: S\in\mathcal F\})$,
order $m^q$.  With $\RR^m=\langle\mathbf 1\rangle\oplus\mathbf 1^\perp$ the tensor product
splits into $\mathcal U_T=\bigotimes_{i\in T}\mathbf 1^\perp$ of dimension $(m-1)^{|T|}$, and
$$A(G_m)\big|_{\mathcal U_T}=\sum_{S\in\mathcal F}(-1)^{|S\cap T|}(m-1)^{|S\setminus T|}.$$
Differentiating the multilinear expansion of $P$ at $x_i=1/m$ identifies this with
$\lambda_T(m)=\sum_{T\subseteq R}c_Rm^{q-|R|}$.

Leading term analysis: the leading coefficient of $\lambda_T$ is $c_T$ if $c_T\ne0$, else
$c_R$ for the smallest $R\supseteq T$ with $c_R\ne0$.  Every $T$ has some such $R$ because
$c_{[q]}>0$, so $n^0(G_m)=0$ for large $m$.  Negative eigenvalues occur only for $|T|\le s$.

## 5. Certified data

$$s_{\min}(2),\dots,s_{\min}(7)=1,2,2,3,3,3,\qquad q/s_{\min}=2,\tfrac32,2,\tfrac53,2,\tfrac73.$$

All satisfy $q<3s$.  Exact inertias: $(7,3)$ at $m=2$ gives $(104,0,24)$ on $128$ vertices;
$(2,1)$ gives $(2,0,2),(5,0,4),(10,0,6),(17,0,8)$ for $m=2,3,4,5$.

Uncertified (time-limited) runs: $s_{\min}(8)=4$, $s_{\min}(9)=5$, $s_{\min}(10)=5$,
$s_{\min}(11)=7$, $s_{\min}(12)=8$.

## 6. Dead ends (kept deliberately)

**(a) Kneser graphs on $r$-subsets.**  Replace the $2$-subsets of $W(k)$ by $r$-subsets:
points $\{a_i\}$ inducing $K_k$, $b_S$ inducing $KG(k,r)$, $a_i\sim b_S$ iff $i\in S$.
The hope was $\pp=\binom{k}{r}\approx 2^k$.  With incidence matrix $C$ ($k\times N$,
$N=\binom kr$) the block matrix is $\begin{pmatrix}J_k-I_k & C\\ C^{\mathsf T}& K\end{pmatrix}$
and **only for $r=2$** is $K=J_N+(r-1)I_N-C^{\mathsf T}C$ the Kneser adjacency: for
$S\ne T$ with $|S\cap T|=t$ that expression has entry $1-t$, not $\mathbf 1_{t=0}$.
Numerically, $k=6,r=4$ gives inertia $(6,6,9)$, not the predicted $(15,\ldots)$.  And the
predicted $\pp=\binom kr$ is impossible anyway: the $b$-vertices are an independent set of
size $\binom kr$, while $\pp(G)\le n-\alpha(G)=k$.

**(b) The tempting general bound $p\le\lfloor n/2\rfloor$.**  It holds for $n\le7$ and is
refuted at $n=8$ by the Wagner graph $M_8=C_8+\{$antipodal edges$\}$, spectrum
$3,\,1,\,1,\,\sqrt2-1\ (\times2),\,-1,\,-(1+\sqrt2)\ (\times2)$: five positive eigenvalues
on eight vertices.  Not part of this paper's contribution, but it is the natural first
guess and it is wrong, which is why the $q/s$ framing (not a ratio of $n$) is the right one.

**(c) An induction giving $p(G)\le\lfloor n/2\rfloor+\alpha(G)$.**  The induction step is
invalid: one needs an *upper* bound for $\alpha(G-v)$ but only has $\alpha(G-v)\ge\alpha(G)-1$.
This is why the seemingly clean constant $3/4$ (from combining $p\le n-\alpha$ with
$p\le n-2+\alpha$) is not available.

## 7. The recurrence worth exploiting next

Splitting $g$ along the last coordinate, $g_0(S)=g(S)$, $g_1(S)=g(S\cup\{q\})$ on
$2^{[q-1]}$, gives with $\widehat{\cdot}$ the Möbius transform on $2^{[q-1]}$:
$$c_{T\cup\{q\}}=\widehat{g_1}(T),\qquad c_T=\widehat{g_0}(T)-\widehat{g_1}(T).$$
Admissibility then demands $\widehat{g_1}([q-1])\ge1$, $\widehat{g_1}(T)\ge0$ for
$|T|\ge s$, and $\widehat{g_0}(T)\ge\widehat{g_1}(T)$ for $|T|>s$, with
$g_1([q-1])=0$.  This looks like the right handle for a divide-and-conquer construction.