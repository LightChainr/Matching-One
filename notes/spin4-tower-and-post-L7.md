# Thermal spin-4 tower versus the historical post-L^-7 annihilator

Status: theory/numerics interface note.  The main purpose is to prevent two distinct observations from being conflated.

## 1. What is now strongly supported

The same-N square-site/matching data support a matching-odd spin-4 central term of the form

\[
\Delta M_N(p_c)\simeq A_4\,\Delta\cos(4\theta)\,N^{-13/8}
=A_4\,\Delta\cos(4\theta)\,L^{-13/4}.
\]

The parameter-free Gaussian doubling test independently checks the joint sign and radial law because multiplication by `1+i` sends

\[
N\mapsto2N,\qquad \theta\mapsto\theta+\pi/4,
\]

and therefore predicts

\[
\Delta M_{2N}/\Delta M_N=-2^{-13/8}.
\]

The fresh-seed score on PR #21 is compatible with this fixed ratio on both frozen lineages.

## 2. Thermal-family interpretation

For percolation the thermal primary has

\[
h=\bar h=5/8,\qquad x_t=5/4.
\]

At `c=0,h=5/8` the first null relation is at chiral level 2.  The next Kac singular level is 10.  Up to level 9 the quotient character is therefore

\[
d_n=p(n)-p(n-2).
\]

Subtracting the `L_-1` image gives quasiprimary counts

```text
level n       0 1 2 3 4 5 6 7 8 9
q_n           1 0 0 1 1 1 2 2 3 4
```

(see `scripts/thermal_spin4_tower.py`).

A non-redundant bulk spin-4 field formed from chiral quasiprimaries requires levels `(n,m)` with `n-m=4` and both `q_n,q_m>0`.

The first allowed pair is

\[
(4,0),\qquad n+m=4,
\]

which gives

\[
x=x_t+4=21/4,
\]

and hence

\[
M(p_c)\sim L^{-13/4},\qquad p_L^*-p_c\sim L^{-4}.
\]

The would-be total levels 5, 6, 7, 8, 9 cannot give another non-redundant spin-4 quasiprimary pair because the opposite chiral levels would be 1 or 2, where no quasiprimary survives.

The next allowed pair is

\[
(7,3),\qquad n+m=10,
\]

which gives

\[
x=x_t+10=45/4.
\]

Therefore, **within the same thermal conformal family**, the first ordinary non-logarithmic radial correction to the leading spin-4 amplitude is separated by six powers of `L`:

\[
\Delta M=L^{-13/4}\left(A+B L^{-6}+\cdots\right)
=L^{-13/4}\left(A+B N^{-3}+\cdots\right).
\]

This is a useful structural prediction.  It explains why adding an arbitrary `N^-1` correction to the current orientation law has no compelling CFT motivation from the same family.

Caveat: other conformal families, logarithmic partners and nonlinear scaling fields can generate earlier corrections.  This note only classifies the ordinary quasiprimary tower of the thermal family.

## 3. Why this does NOT explain the historical L^-7 estimator

Mertens and Ziff (2016) considered the accelerated criterion

\[
L^{13/4}M_L(p^*)=(L-1)^{13/4}M_{L-1}(p^*)
\]

and observed an apparent convergence close to `L^-7.06` on small systems.

Write the central amplitude generally as

\[
M_L(p_c)=aL^{-13/4}\left(1+cL^{-q}+\cdots\right),
\]

while

\[
M'_L(p_c)=bL^{3/4}(1+\cdots).
\]

A first-order expansion of the two-size criterion gives

\[
p^*_{\rm ann}(L)-p_c\propto L^{-(4+q)}.
\]

Hence a true `L^-7` accelerated root would mean

\[
q=3.
\]

But the next ordinary spin-4 quasiprimary in the **same thermal family** gives `q=6`, which would instead imply an `L^-10` contribution after leading annihilation.

Therefore the historical exponent near 7, if asymptotically real, cannot simply be called “the next thermal spin-4 descendant.”  Plausible alternatives are:

1. a different matching-odd conformal family or scalar sector with relative exponent 3;
2. a composite/nonlinear correction;
3. logarithmic mixing which mimics an effective power on small sizes;
4. a relation to the percolation correction-to-scaling length exponent `omega=3/2`, possibly at second order (`2 omega=3`) if its linear contribution is absent for this observable;
5. the 2016 `7.06` being preasymptotic.  The original paper itself emphasized strong small-L distortion of apparent exponents.

The `omega=3/2` connection is deliberately labelled conjectural.  The exact cluster-size correction exponent `Omega=72/91` implies a length exponent `omega=3/2`, but it is not automatic that the same correction field appears in this torus matching observable.

## 4. Modern falsification experiment

The repository now has a fast C++ threshold-rank engine, so the 2016 post-L^-7 observation should be remeasured rather than explained from old small-L data.

For ordinary axis square tori reconstruct the full `M_L(p)` and solve

\[
L^{13/4}M_L(p)=(L-1)^{13/4}M_{L-1}(p)
\]

for a sequence of substantially larger L.

Predeclare competing accelerated exponents

\[
w_{\rm ann}\in\{11/2,6,7,8,10\}
\]

corresponding respectively to relative central corrections

\[
q\in\{3/2,2,3,4,6\}.
\]

Also allow a free exponent and explicit log alternatives, but select only on training sizes and score on a frozen tail.

The primary scientific question is no longer “can we get another digit of p_c?” It is:

> Which correction sector survives after the now-identified `x=21/4` leading term is annihilated?

A clean answer would directly constrain the LCFT/operator content beyond the leading spin-4 field.
