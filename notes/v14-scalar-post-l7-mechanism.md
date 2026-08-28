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
L^{2-33/4}=L^{-25/4}=L^{-13/4}L^{-3}.
\]

Hence its relative correction is exactly `q=3`. For the Mertens--Ziff leading-annihilated root, whose leading root power is 4, this produces

\[
\boxed{p^*_{\rm ann}(L)-p_c\sim L^{-7}}.
\]

This supplies a concrete standard-spectrum mechanism for the historical `~L^-7.06` observation. It does not prove that this mechanism dominates the finite-size sequence.

## 3. Distinction from the next thermal H4 descendant

The next ordinary nonredundant spin-4 quasiprimary in the same thermal family has `x=45/4`. It contributes

\[
L^{-37/4}=L^{-13/4}L^{-6},
\]

so `q=6` and the accelerated root is `L^-10`.

```text
V_<1,4>:       H0 scalar, q=3, w_ann=7
thermal next:  H4,        q=6, w_ann=10
```

A stable modern `w≈7` result would therefore point toward a scalar correction sector rather than simply the next ordinary thermal H4 descendant.

## 4. Matching parity is conditional

The finite matching identity alone does not assign a local CFT parity to V14. A stronger hypothesis is required: matching/complement extends to the relevant interchiral/OPE structure, sends `V_<1,2>` to minus itself, and does not rotate degenerate equal-quantum-number channels. Under that hypothesis

\[
\eta_s=(-1)^{s-1},\qquad \eta(V_{\langle1,4\rangle})=-1.
\]

Claim hierarchy:

- exact: V14 exists, is a singlet scalar, `x=33/4`;
- exact: if it enters the central matching difference linearly, it gives `q=3 -> w=7`;
- conditional: the stronger interchiral matching automorphism makes it matching-odd;
- unresolved: the lattice coupling is nonzero and explains the observed sequence.

Issue #61 remains the appropriate RG-level theory task.

## 5. Two-angle scalar nulling is only a truncation diagnostic

For a truncated model

\[
D(\theta)=A_0+A_4\cos4\theta,
\]

two orientations can eliminate H4 and recover A0. With H8/H12/H16/... present, that combination is **not** an exact all-order scalar projector.

The exact finite-harmonic solution now on `main` is the four-angle N=1105 design

```text
(33,4), (32,9), (31,12), (24,23)
```

which inverts H0/H4/H8/H12 exactly. N=1105=5*13*17 is also the smallest primitive Gaussian torus with four D4-inequivalent orientations. H16 and higher harmonics remain possible.

Therefore the two-angle retrospective result in this PR is only a bounded power/leakage diagnostic. It is not used as an operator claim.

## 6. Why the adjacent-axis annihilator is the efficient primary test

A scalar H0 contribution cancels from a same-N orientation **difference**. Chasing larger orientation-difference simulations is therefore a poor way to detect V14.

The coupled adjacent-axis statistic in this PR instead measures

\[
F_L(p)=L^{13/4}M_L(p)-(L-1)^{13/4}M_{L-1}(p).
\]

At a nearby fixed p,

\[
F_L=C_q[L^{-q}-(L-1)^{-q}] + T[L^4-(L-1)^4]+\cdots,
\]

so the q=3 scalar mechanism can be challenged without supplying the last digits of pc. This is the primary numerical route.

## 7. Gaussian-semigroup discriminator for a future H0 measurement

If an H0 component is reconstructed at Gaussian-related sizes, a pure V14 scalar transforms with no angular sign:

\[
H0_{QN}/H0_N=Q^{-25/8}.
\]

The next thermal H4 raw contrast transforms as

\[
\Delta H4_{QN}/\Delta H4_N=r_4(h)Q^{-37/8}.
\]

For h=1+i these are respectively `+2^-25/8` and `-2^-37/8`.

## 8. Interpretation order

The modern axis challenge should compare the fixed mechanisms already frozen in the protocol:

1. q=3 / w=7: V14 H0 candidate;
2. q=2 / w=6: lower ordinary/composite correction;
3. q=4 / w=8: nonlinear H4/H12 sideband;
4. q=6 / w=10: next ordinary thermal H4 quasiprimary;
5. q=3/2 historical adversary.

If q=3 wins on the held-out tail, the next high-value step is an independent H0 amplitude measurement, not a free exponent retune.

Reference for the Potts spectrum/interchiral structure: Jacobsen--Ribault--Saleur, arXiv:2208.14298. The matching-parity assignment remains conditional as stated above.
