# C3 homology character versus embedded physical spin four away from the hexagonal modulus

2026-09-14.  Typing clarification for the existing #156 primitive-sector pilot.  No previous finite numbers are erased or rescored here.

The three-line contrast used in #156 is a valid **homology-representation character**.  The new projective-slope analysis shows that it equals the real part of the physically embedded spin-four harmonic only at the exact hexagonal modulus.  Away from that fixed point, the two are distinct observables and should be named separately.

## 1. The existing three-line character

In the positive-`rho` convention of #156, the distinguished primitive lines are

\[
\ell_0=(1,0),\qquad
\ell_1=(0,1),\qquad
\ell_2=(1,-1).                                                  \tag{1.1}
\]

The 60-degree homology action cycles these three lines at the Eisenstein/hexagonal point.  With probabilities `P_0,P_1,P_2`, the real nontrivial C3 character used by #156 is

\[
\boxed{C_{C3}=P_0-\frac{P_1+P_2}{2}.}                          \tag{1.2}
\]

This is the real part of the character with weights

\[
1,e^{2\pi i/3},e^{-2\pi i/3}.                                 \tag{1.3}
\]

It is an exact representation-theoretic coordinate for the declared three-line orbit.

## 2. Physical spin-four weights depend on the actual modulus

For an embedded torus with

\[
\tau=\frac12+i y,
\]

the projective physical spin-four readout is

\[
Z_4(u,v)=\left(\frac{u+v\tau}{|u+v\tau|}\right)^4.            \tag{2.1}
\]

On the three distinguished lines,

\[
Z_4(1,0)=1,                                                     \tag{2.2}
\]

while the other two are complex conjugates.  Reflection symmetry gives `P_1=P_2`, so their imaginary parts cancel and the three-line physical contribution is

\[
\boxed{
A_4^{(3)}(\tau)
=P_0+c_4(y)(P_1+P_2),}                                         \tag{2.3}
\]

where

\[
c_4(y)=\Re\left[
\left(\frac{1/2+i y}{\sqrt{1/4+y^2}}\right)^4\right].          \tag{2.4}
\]

The C3 contrast (1.2) corresponds instead to the fixed coefficient `-1/2`.

Therefore

\[
\boxed{
A_4^{(3)}(\tau)-C_{C3}
=\left(c_4(y)+\frac12\right)(P_1+P_2).}                        \tag{2.5}
\]

The two coordinates coincide **if and only if** the relevant geometric coefficient is `c_4=-1/2` (apart from a trivial zero sector probability).

## 3. Exact coincidence at the hexagonal fixed point

For

\[
y=\sqrt3/2,
\qquad
\tau=e^{i\pi/3},                                               \tag{3.1}
\]

the two non-axis physical directions make angles `+/- pi/3`.  Hence

\[
c_4=\cos(4\pi/3)=-\frac12.                                   \tag{3.2}
\]

Thus

\[
\boxed{A_4^{(3)}(\tau_{hex})=C_{C3}.}                          \tag{3.3}
\]

This is the geometric reason the three-line C3 real character can be described as a spin-four angular contrast **at the exact Eisenstein modulus**.

The full primitive-sector harmonic nevertheless vanishes there after including the entire sixfold orbit structure, as proved by the 60-degree modular automorphism in `projective-slope-modular-covariance-20260914.md`.

## 4. N30 and N56 are not at the exact fixed point

The PR #213 baselines use

\[
\tau_{30}=\frac12+\frac56 i,
\qquad
\tau_{56}=\frac12+\frac78 i.                                  \tag{4.1}
\]

For these moduli,

\[
\boxed{c_4(5/6)=-0.5570934256055362\ldots,}                    \tag{4.2}
\]

\[
\boxed{c_4(7/8)=-0.4844970414201184\ldots.}                    \tag{4.3}
\]

They are close to, but not equal to, `-1/2`.

Using the already-merged PR #213 primitive probabilities:

### N30

\[
P_0=0.1107291776903850\ldots,
\qquad
P_1=P_2=0.1272155037499346\ldots.                              \tag{4.4}
\]

The three-line C3 character is

\[
C_{C3}=-0.01648632605954963\ldots,                             \tag{4.5}
\]

whereas the physically embedded three-line spin-four contribution is

\[
\boxed{A_4^{(3)}=-0.0310126638579850\ldots.}                   \tag{4.6}
\]

The full all-primitive-sector continuum harmonic computed in the 2026-09-14 control is

\[
A_4^{full}=-0.03047986387219586\ldots.                         \tag{4.7}
\]

Thus the `rank1_other` sectors supply a small but nonzero correction of about `+5.33e-4` to the embedded physical harmonic.

### N56

\[
P_0=0.1247005965715636\ldots,
\qquad
P_1=P_2=0.1201521157922650\ldots.                              \tag{4.8}
\]

Here

\[
C_{C3}=+0.004548480779298590\ldots,                            \tag{4.9}
\]

while

\[
\boxed{A_4^{(3)}=+0.00827390732812383\ldots,}                  \tag{4.10}
\]

and the full primitive-sector value is

\[
A_4^{full}=+0.008136416970836356\ldots.                        \tag{4.11}
\]

The other primitive sectors again shift the physical harmonic slightly.

## 5. Interpretation of the old pilot

The continuum-subtracted N30/N56 pilot in #156 remains a valid measurement of the **frozen C3 homology character** `C_{C3}` and its finite-size residual.  This note does not reinterpret or invalidate that measurement.

What changes is the allowed wording:

- `C_{C3}` may be called the physical embedded spin-four three-line contrast at the exact hexagonal modulus;
- away from that modulus, including N30/N56, it is a fixed homology-character coordinate, not literally `Re e^{i4 theta}` in the physical embedding;
- the full physical spin-four harmonic additionally uses modulus-dependent weights and all primitive rank-one lines.

Therefore a sign/radial transport statement about `C_{C3}` is a statement about that representation coordinate.  It cannot be promoted to an embedded-H4 amplitude without the modulus-dependent map (2.3) and the `rank1_other` contribution.

## 6. Why the distinction matters for H4/H8 alias discussions

On the three exact Eisenstein orbit lines, real C3 characters cannot distinguish certain spin-four/spin-eight aliases because both restrict to the same finite character pattern up to the declared convention.  The existing #156 note already recognized this.

The projective embedded harmonic supplies the missing extra structure:

1. away from the hexagonal point, the physical angle weights move continuously with `tau`;
2. other primitive slopes have distinct `e^{i4 theta}` and `e^{i8 theta}` values;
3. the full Pinson--Arguin sum therefore produces genuinely different spin-four and spin-eight modular-covariant functions.

For example at `tau=i`,

\[
H_4^{cont}=0.7791813140\ldots,
\qquad
H_8^{cont}=0.9992417096\ldots.                                \tag{6.1}
\]

while both vanish at the hexagonal fixed point by the appropriate automorphism selection rules.

This gives a representation-faithful route to break the three-line alias **without** inventing an untyped microscopic source.

## 7. Recommended naming discipline

Use three distinct labels in future notes/data:

```text
C_C3              fixed three-line homology character,
A4_three(tau)     physical spin-4 projection on those three lines,
A4_full(tau)      physical spin-4 sum over every primitive rank-one line.
```

At `tau_hex`, `C_C3=A4_three`, but `A4_full` still contains the remaining primitive orbit structure and is constrained by the full torus automorphism.

This naming prevents a correct representation proxy from silently acquiring a stronger physical-spin interpretation than its definition supports.
