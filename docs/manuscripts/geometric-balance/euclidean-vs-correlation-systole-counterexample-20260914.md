# Euclidean systole need not minimize the correlation norm

2026-09-14.  Explicit period-lattice family showing that the warning in #765 is genuine outside the exponentially elongated regime.  The example uses the rigorous dilute directional asymptotics from `dilute-directional-geodesic-entropy-20260914.md`.

## 1. The lattice family

Let

\[
\Lambda_n=\operatorname{span}_{\mathbb Z}
\{u_n,v_n\},
\]

with

\[
u_n=(10n,0),
\qquad
v_n=(5n,9n).                                                    \tag{1.1}
\]

Its area is

\[
N_n=|\det(u_n,v_n)|=90n^2.                                    \tag{1.2}
\]

Every lattice vector is

\[
a u_n+b v_n
=n(10a+5b,9b).                                                 \tag{1.3}
\]

If `b=0`, the shortest nonzero vector has length `10n`.  If `|b|>=1`, minimizing over `a` gives

\[
|10a+5b|\ge
\begin{cases}
0,&b\text{ even},\\
5,&b\text{ odd},
\end{cases}                                                     \tag{1.4}
\]

so

\[
|a u_n+b v_n|^2/n^2
\ge
\begin{cases}
81b^2,&b\text{ even},\\
25+81b^2,&b\text{ odd}.
\end{cases}                                                     \tag{1.5}
\]

For `b=+/-1` this is at least `106>100`; for `|b|>=2` it is at least `324`.  Therefore

\[
\boxed{\ell(\Lambda_n)=10n,\qquad
\text{the unique Euclidean shortest line is }\mathbb R u_n.}  \tag{1.6}
\]

## 2. Matching dilute correlation costs

The fixed-direction dilute theorem gives, for every fixed integer vector `x`,

\[
\tau_{8,p}(x)
=\|x\|_\infty\log(1/p)+O_x(1),
\qquad p\downarrow0,                                          \tag{2.1}
\]

where the `O(1)` term is the geodesic-entropy correction.

For the Euclidean systole,

\[
\tau_{8,p}(u_n)
=10n\log(1/p)+O(n).                                            \tag{2.2}
\]

But

\[
v_n=(5n,9n),
\qquad
v_n-u_n=(-5n,9n),                                              \tag{2.3}
\]

so both have

\[
\|v_n\|_\infty=\|v_n-u_n\|_\infty=9n.                        \tag{2.4}
\]

Hence

\[
\tau_{8,p}(v_n)
=9n\log(1/p)+O(n),                                             \tag{2.5}
\]

\[
\tau_{8,p}(v_n-u_n)
=9n\log(1/p)+O(n).                                             \tag{2.6}
\]

For all sufficiently small fixed `p`, the `log(1/p)` advantage beats the bounded-per-`n` entropy constants, giving

\[
\boxed{
\tau_{8,p}(v_n)<\tau_{8,p}(u_n),
\qquad
\tau_{8,p}(v_n-u_n)<\tau_{8,p}(u_n).}                         \tag{2.7}
\]

The two cheaper lines are nonparallel and form another basis of the same lattice:

\[
|\det(v_n,v_n-u_n)|=90n^2=N_n.                                \tag{2.8}
\]

Thus the matching correlation norm has **two nonparallel cheapest primitive classes** while the Euclidean norm has a unique shortest class.

## 3. NN behaves differently in the same family

For NN site percolation, the dilute directional theorem instead gives

\[
\tau_{4,p}(x)
=\|x\|_1\log(1/p)+O_x(1).                                     \tag{3.1}
\]

Here

\[
\|u_n\|_1=10n,
\qquad
\|v_n\|_1=\|v_n-u_n\|_1=14n.                                \tag{3.2}
\]

Therefore the Euclidean axis systole `u_n` is also the cheapest NN class at sufficiently small `p`.

The same finite period lattice can therefore have different rare winding directions for NN and matching connectivity.

## 4. Why this does not contradict the exponential-elongation separation lemma

For this family,

\[
\frac{N_n/\ell_n}{\ell_n}
=\frac{9n}{10n}=0.9,                                          \tag{4.1}
\]

so the transverse height is only of the same order as the systole.  The deterministic direction-separation result in `structural-consequences-20260914.md` assumes

\[
h/\ell\to\infty,                                              \tag{4.2}
\]

which holds in the fixed-positive exponential-aspect regime but **not** here.

Under (4.2), every nonparallel period has Euclidean length at least `h`, and norm equivalence makes its correlation cost diverge relative to the shortest vector.  That special geometry legitimately restores the Euclidean shortest direction.

So the two statements fit together:

- general period shape: minimize the **correlation norm**, not Euclidean length;
- exponential elongation: determinant geometry forces every nonparallel class to be so long that the Euclidean systole automatically wins.

## 5. Consequence for #765 theorem statements

A correct general directional centre theorem should be phrased using

\[
\boxed{
\min_{v\in\Lambda\setminus\{0\}}\tau_{G,p}(v),}               \tag{5.1}
\]

with the projective minimizer set retained when it is nonunique.

Only after a geometric condition such as `h/ell -> infinity` has been imposed may this be reduced to

\[
\tau_{G,p}(u_{euclidean\ systole}).                            \tag{5.2}
\]

The counterexample shows that this distinction is mathematical, not merely conservative wording.

## 6. Multi-direction implication

At small matching `p`, the present family has two equal-cost nonparallel candidate slopes.  The exact projective-homology classification says the rank-one sector may carry only **one** line at a time; a second nonparallel essential class forces rank two.

Therefore this lattice family is a natural finite-shape testbed for the projective hard-core-gas picture in `projective-homology-gas-20260914.md`:

- first birth: a race between the two cheap projective lines;
- rank-one plateau: selected line remains frozen;
- second nonparallel winding: immediate transition to rank two.

No independent two-species Poisson model can represent all three statements simultaneously.
