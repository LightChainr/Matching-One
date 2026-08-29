# Hermite--Krawtchouk thermal jet and the rank-gap squeeze

Status: exact finite-`N` generating-function identity, followed by conditional
scaling/operator hypotheses frozen before the norm-5 score-mode reveal.

## 1. The score tower is one full-curve generating function

Fix the intrinsic center `p0` and let `H_r` be the repository's orthonormal
`Bin(N,p0)` Krawtchouk basis.  For a microcanonical response `q_N(k)`, write

\[
c_{r,N}=\mathbb E_{p_0}[q_N(K)H_r(K)].
\]

The positive-score convention obeys the exact identity

\[
\mathbb E_p H_r(K)
=\sqrt{\binom Nr}
\left(\frac{p-p_0}{\sqrt{p_0(1-p_0)}}\right)^r.
\]

Completeness of the finite Krawtchouk basis therefore gives

\[
\boxed{
R_N(p)=\sum_{r=0}^N c_{r,N}\sqrt{\binom Nr}
\left(\frac{p-p_0}{\sqrt{p_0(1-p_0)}}\right)^r.
}
\]

Thus the modes are not several unrelated observables.  They are the exact
Taylor-generating coordinates of the same canonical full curve.

For the matching-odd thermal tower define `c_r` by taking `P4[D_r]` at even
`r` and `P4[S_r]` at odd `r`.  If

\[
c_{r,N}\sim N^{-13/8-r/8},
\qquad z=(p-p_0)N^{3/8},
\]

then the scaling-function derivative jet is

\[
\boxed{
d_{r,N}=\frac{\sqrt{r!}\,N^{13/8+r/8}c_{r,N}}
{[p_0(1-p_0)]^{r/2}},
\qquad
F_T(z)=\sum_{r\ge0}\frac{d_r}{r!}z^r.
}
\]

The Krawtchouk-to-Hermite limit is therefore an operator statement: score
order `r` is the `r`-th thermal Taylor component of one matching-odd spin-4
scaling function.  It is not a new CFT primary at every rung.

## 2. Three mechanisms are three operators on the jet

They have sharply different actions.

### Translation

For `F(z)->F(z+delta)`,

\[
d'_r=\sum_{k\ge0}\frac{\delta^k}{k!}d_{r+k},
\qquad (T d)_r=d_{r+1}.
\]

Translation mixes adjacent modes and generally rotates the signed direction
of a truncated jet.

### Width or thermal-metric renormalization

For `F(z)->F(z/w)`,

\[
d'_r=w^{-r}d_r,
\qquad (W d)_r=-r d_r.
\]

Width flow preserves every mode sign and produces a correction proportional
to the score order.  The stable signed high-mode direction together with a
drifting magnitude is therefore more naturally compatible with width flow
than with a large untracked translation.  This is only a mechanism ranking;
the frozen score below decides it.

### Jordan/log mixing

A rank-2 block gives

\[
d_{r,N}=f_r+\log N\,g_r+\cdots.
\]

The correction vector `g` need not be proportional to either `Td` or `Wd`.
Its identifying property is instead the additive multiplier cocycle.  The
alternating signs of the observed jet are properties of `F_T`; by themselves
they are not evidence of a Jordan block.

## 3. The paired rank gap supplies an independent width observable

The joint threshold statistic

\[
G=K_+-K_-
\]

measures the number of occupation ranks in the topologically neutral window.
This statement has an exact canonical bridge, not only a scaling
interpretation.  For each permutation define

\[
m_k=\begin{cases}
-1,&k<K_-,\\
0,&K_-\le k<K_+,\\
+1,&k\ge K_+,
\end{cases}
\qquad U_k=1-m_k^2.
\]

Canonical binomial mixing gives

\[
U_N(p)=\mathbb E[U_{K}],\qquad K\sim\operatorname{Bin}(N,p).
\]

Every degree-`N` Bernstein basis function has the same unit-interval area,

\[
\int_0^1 \binom Nk p^k(1-p)^{N-k}\,dp=\frac1{N+1}.
\]

Because exactly `K_+-K_-` layers have `U_k=1`, linearity over permutations
proves the finite-size identity

\[
\boxed{
\int_0^1 \mathbb E\!\left[U_{K\sim\operatorname{Bin}(N,p)}\right]dp
=\frac{\mathbb E[K_+-K_-]}{N+1}.
}
\]

Thus the expected rank gap is exactly `(N+1)` times the canonical neutral-
window area.  Calling it a width is therefore an observable identity before
any scaling ansatz is imposed.

The thermal law gives

\[
\mathbb E G\sim A N^{5/8}.
\]

The new source fit indicates the sharper form

\[
\mathbb E G=A N^{5/8}+B+\cdots.
\]

In probability units this is

\[
\frac{\mathbb E G}{N}
=N^{-3/8}\left[A+B N^{-5/8}+\cdots\right].
\]

Hence

\[
\boxed{w_N=N^{-5/8}\mathbb E G=A+B N^{-5/8}+\cdots}
\]

is exactly a finite-size thermal-window width with a relative `N^-5/8`
correction.  This turns the rank-gap observation into a prediction for the
whole Hermite jet:

\[
d_{r,N}\propto w_N^{-r}.
\]

After allowing an arbitrary size-dependent common amplitude, the
division-free held-out residual is

\[
\boxed{
R^{\rm width}_{r,Q}
=d_{r,QN}w_{QN}^{r}d_{0,N}
-d_{r,N}w_N^{r}d_{0,QN}=0,
\quad r=2,\ldots,6.
}
\]

This is the strongest immediate prediction because `w_N` comes from paired
joint rank moments, whereas `d_r` comes from the marginal full-curve score
modes.  They share replicas and must be combined inside the same delete-one
replicates, but the width law is not fitted from the mode vector.

## 4. Why `B=-1/4` is not yet an endpoint theorem

The fitted constant is close to `-1/4`, but the frozen rank convention alone
does not produce it.  Both transitions use the same black occupation rank:
the neutral plateau contains exactly

\[
K_+-K_-=G
\]

integer levels.  Moving both boundaries from integer ranks to half-ranks
shifts both endpoints equally and leaves `G` unchanged.  Changing only the
reverse off-by-one changes `G` by an integer.  With two equally pooled
orientations, integer convention changes can produce integer or half-integer
mean shifts, not a forced quarter.

The exact order-statistic bridge also gives

\[
\mathbb E(T_+-T_-)=\frac{\mathbb E G}{N+1},
\]

which does not generate an `O(1)` quarter in `E[G]`.  A true `-1/4` would need
an additional Euler--Maclaurin or topological boundary-weight derivation.
Until then it is a sharp conjectural value, not an exact consequence of the
stored endpoint definitions.

## 5. Width first, then the q2/Jordan cocycle

Define the width-corrected jet

\[
\widetilde d_{r,N}=w_N^r d_{r,N}.
\]

Score width-only collapse first.  If a residual functional direction remains,
use the already protected model order

\[
R_c(r)=\widetilde d_{r,5N}
-c\widetilde d_{r,2N}
+(c-1)\widetilde d_{r,N},
\]

with

\[
c=\frac85\quad\text{before}\quad
c=\frac{\log5}{\log2}.
\]

The first value is the ordinary relative-`q=2` correction; the second is the
rank-2 Jordan cocycle.  A Jordan claim becomes credible only if a coherent
width-corrected vector over `r=2..6` follows the logarithmic multiplier, not
because several raw coefficients alternate in sign.

## 6. Reproduction

```bash
python3 scripts/hermite_krawtchouk_scaling_jet.py
python3 tests/test_hermite_krawtchouk_scaling_jet.py -v
```

The frozen contract is
`predictions/hermite_krawtchouk_jet_20260829.yaml`.
