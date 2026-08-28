# A named scalar mechanism for the historical post-L^-7 annihilator

Status: exact spectrum arithmetic plus a conditional E3 matching-parity mechanism.  This note does not claim that the historical small-size exponent is asymptotic or that the relevant lattice coupling is nonzero.

## 1. The scalar omitted by a spin-4-only search

On the critical Potts branch at percolation,

\[
h_{r,s}=\frac{(2r-3s)^2-1}{24}.
\]

The Potts space of states contains every diagonal degenerate module
`R_<1,s>` in the trivial `S_Q` representation.  For `s=4`,

\[
h_{1,4}=\frac{33}{8},\qquad
x_{1,4}=2h_{1,4}=\frac{33}{4},\qquad s_{\rm conf}=0.
\]

Thus `V_<1,4>` is an ordinary diagonal Potts singlet scalar.  Searching only
for a spin-4 field at `x=33/4` misses it.

## 2. Conditional matching parity

If matching/complement extends to an OPE/interchiral automorphism, reverses the
thermal generator `V_<1,2>`, and does not rotate degenerate equal-quantum-number
channels, the degenerate fusion recursion gives

\[
\eta_s=(-1)^{s-1}.
\]

Consequently

\[
\eta(V_{\langle1,4\rangle})=-1.
\]

This parity statement is E3 under the explicit automorphism assumption, not an
E0 consequence of the finite matching identity.

## 3. Exact exponent arithmetic

The leading proposed thermal spin-4 contribution has

\[
x_4=\frac{21}{4},\qquad M_L(p_c)\sim L^{-13/4}.
\]

A linear `V_<1,4>` correction contributes to a dimensionless torus observable
as

\[
L^{2-x_{1,4}}=L^{-25/4}
=L^{-13/4}L^{-3}.
\]

Therefore it supplies exactly the relative correction

\[
q=x_{1,4}-x_4=3.
\]

After annihilating the leading `L^-13/4` amplitude, division by the thermal
slope `M'_L(p_c)\sim L^{3/4}` gives

\[
p^*_{\rm ann}(L)-p_c\sim L^{-7}.
\]

The arithmetic is executable in `scripts/v14_scalar_post_l7.py`.

## 4. Distinguishing it from the thermal level-(7,3) correction

The next ordinary spin-4 quasiprimary pair inside the same thermal family is
at chiral levels `(7,3)`, total level 10.  It gives relative `q=6` and an
`L^-10` accelerated root.  The two mechanisms also differ angularly:

- `V_<1,4>` is scalar (`H0`) and should cancel in a same-modulus orientation
  difference;
- the thermal `(7,3)` correction remains spin 4 (`H4`).

An unprojected axis sequence can contain both.  Same-N orientation projection
separates them without fitting a free exponent.

More explicitly, for two orientations at the same `N`, write schematically

\[
D_i=A_0N^{-25/8}+c_{4,i}
\left(A_4N^{-13/8}+B_4N^{-37/8}\right)+\cdots .
\]

Then

\[
P_0[D]=\frac{c_{4,1}D_2-c_{4,2}D_1}{c_{4,1}-c_{4,2}}
\]

cancels both displayed `H4` terms and retains the `V_<1,4>` scalar term.
This is a direct scalar projector, not an exponent fit.

Gaussian doubling provides a second parameter-free discriminator.  Under
`N -> 2N` and `theta -> theta+pi/4`, the two candidate child/parent ratios are

```text
V_<1,4> scalar H0:             +2^(-25/8)
thermal level-(7,3) spin4 H4:  -2^(-37/8)
```

They differ in both sign and radial power.

## 5. Claim boundary and numerical target

Evidence levels:

- E0 within the stated Potts spectrum: `V_<1,4>` exists, is a singlet scalar,
  and has `x=33/4`;
- E3 conditional: it is matching odd and is therefore allowed in the central
  matching difference;
- E4 unresolved: the lattice/observable coupling is nonzero and explains the
  historical apparent `L^-7` law.

Issue #47 should keep the preregistered `w_ann=7` model but identify this
scalar as its primary standard-spectrum mechanism.  The decisive comparison is
`V_<1,4>` (`H0`, `w=7`) versus the thermal level-(7,3) sector (`H4`, `w=10`),
with logarithmic/composite/preasymptotic alternatives retained.

References: Jacobsen--Ribault--Saleur, arXiv:2208.14298, especially the Potts
space of states in Eq. (3.10) and the `V_<1,2>` interchiral generator; see also
`notes/critical-potts-kac-convention-correction.md` and
`notes/interchiral-v13-matching-parity.md`.
