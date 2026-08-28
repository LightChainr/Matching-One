# V_<1,4> as a scalar mechanism for the historical post-L^-7 annihilator

Status: exact spectrum/exponent arithmetic plus a conditional matching-parity mechanism. This note does **not** assume that the historical small-size `L^-7.06` observation is already asymptotic, nor that the V14 lattice coupling is nonzero.

## 1. The scalar missed by a spin-4-only search

On the critical Potts branch at percolation,

\[
h_{r,s}=\frac{(2r-3s)^2-1}{24}.
\]

The diagonal degenerate Potts modules `R_<1,s>` are in the trivial internal-symmetry sector. For `s=4`,

\[
h_{1,4}=\frac{33}{8},\qquad
x_{1,4}=\frac{33}{4},\qquad
s_{\rm conf}=0.
\]

Thus `V_<1,4>` is a standard-spectrum singlet **scalar**. A search restricted to the next spin-4 quasiprimary will not see it.

## 2. Exact exponent arithmetic

The leading candidate matching-odd thermal spin-4 term has

\[
x_T=\frac{21}{4},\qquad
M_L(p_c)\sim L^{2-x_T}=L^{-13/4}.
\]

A linear V14 term contributes

\[
L^{2-33/4}=L^{-25/4}
           =L^{-13/4}L^{-3}.
\]

Hence its relative correction is exactly

\[
q=3.
\]

For the Mertens--Ziff leading-annihilated root, whose leading root power is `4`, this produces

\[
\boxed{p^*_{\rm ann}(L)-p_c\sim L^{-(4+3)}=L^{-7}}.
\]

This supplies a concrete standard-spectrum mechanism for the historical `~L^-7.06` observation. It does not prove that this mechanism dominates the finite-size sequence.

The arithmetic is executable in `scripts/v14_scalar_post_l7.py`.

## 3. Why this is different from the next thermal H4 descendant

The next ordinary nonredundant spin-4 quasiprimary in the same thermal family is at total descendant level 10:

\[
x=\frac{45}{4}.
\]

It contributes

\[
L^{-37/4}=L^{-13/4}L^{-6},
\]

so

\[
q=6,\qquad p^*_{\rm ann}-p_c\sim L^{-10}.
\]

The two candidate mechanisms therefore differ in both radial exponent and angular sector:

```text
V_<1,4>:       H0 scalar, q=3, w_ann=7
thermal next:  H4,        q=6, w_ann=10
```

This makes a stable modern `w≈7` result much more informative than an unexplained effective exponent: it would point toward a scalar matching-odd correction sector rather than the next ordinary thermal H4 descendant.

## 4. Matching parity is the conditional step

The finite matching identity alone does not assign a local CFT parity to V14.

A stronger hypothesis is needed: matching/complement extends to an automorphism of the relevant interchiral/OPE structure, sends the thermal generator `V_<1,2>` to minus itself, and does not rotate degenerate channels with identical continuum quantum numbers. Under that hypothesis the degenerate fusion recursion gives

\[
\eta_s=(-1)^{s-1},
\]

and hence

\[
\eta(V_{\langle1,4\rangle})=-1.
\]

So the hierarchy is:

- exact: V14 exists, is a singlet scalar, `x=33/4`;
- exact: if it contributes linearly to the central matching difference, it gives `q=3 -> w=7`;
- conditional: the stronger interchiral matching automorphism makes it matching-odd;
- unresolved: its lattice coupling is nonzero and explains the observed sequence.

Issue #61 remains the correct place to construct the weaker RG-tangent-space matching action without overclaiming an OPE automorphism.

## 5. Orientation projectors: two angles are not enough for an all-order scalar measurement

For a truncated model

\[
D(\theta)=A_0+A_4\cos4\theta,
\]

two orientations can eliminate H4 and recover `A0`. But with H8/H12/H16/... present, that two-angle combination is **not** an exact scalar projector.

The current exact solution is the four-angle N=1105 design now on `main`:

```text
(33,4), (32,9), (31,12), (24,23)
```

which inverts the finite basis `H0,H4,H8,H12` exactly. The number-theory check also proves `N=1105=5*13*17` is the smallest primitive Gaussian torus with four D4-inequivalent orientations.

This still does not annihilate H16 and higher harmonics; the claim is an exact decomposition **through H12**, not an all-orders scalar theorem.

## 6. Gaussian-semigroup discriminator

Once an H0 component is reconstructed at two Gaussian-related sizes, a pure V14 scalar has no angular sign factor. Under a multiplier of norm `Q`,

\[
H0_{QN}/H0_N=Q^{-25/8}.
\]

For the next thermal H4 descendant, a raw H4 orientation contrast instead transforms as

\[
\Delta H4_{QN}/\Delta H4_N
=r_4(h)Q^{-37/8}.
\]

For `h=1+i`, `Q=2` and `r4=-1`, so the two no-fit predictions are

```text
V14 scalar H0:             +2^(-25/8)
next thermal spin4 H4:     -2^(-37/8)
```

They differ in sign, angular character and radial power.

## 7. Recommended modern test

Issue #47 should preserve its leakage-safe accelerated-root challenge, but interpretation should now explicitly include:

1. `w=6`: relative q=2 scalar/nonlinear correction;
2. `w=7`: V14 H0 candidate (conditional matching parity);
3. `w=8`: nonlinear H4/H12 sideband mechanism;
4. `w=10`: next ordinary thermal H4 quasiprimary;
5. logarithmic/free alternatives only as declared competitors.

If `w≈7` wins on a genuinely larger held-out tail, the next task is **not** another free exponent fit. It is to measure an H0 correction independently using an orientation decomposition and test whether its radial amplitude is consistent with `N^-25/8`.

Reference for the Potts spectrum/interchiral structure: Jacobsen--Ribault--Saleur, arXiv:2208.14298. The parity assignment remains conditional as stated above.
