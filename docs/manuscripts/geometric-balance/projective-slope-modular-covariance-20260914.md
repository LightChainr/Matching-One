# A map-resolved modular-covariant spin-4 ray from primitive homology sectors

2026-09-14.  This note turns the projective slope harmonic into a genuine modular-covariant continuum control.  It consumes, rather than replaces, the already-merged Pinson--Arguin primitive-sector evaluator from PR #213.

No continuum field is named.  The covariance follows directly from two facts:

1. modular transformations permute primitive homology sectors with their published probabilities;
2. the embedded angle of the same physical cycle rotates by the coordinate rescaling used to renormalize the torus basis.

This gives #585 a concrete map-resolved modular-covariant ray whose lattice semantics are explicit before any fit.

## 1. Continuum slope harmonic

For

\[
T_\tau=\mathbb C/(\mathbb Z+\tau\mathbb Z),\qquad \Im\tau>0,
\]

let the engine primitive line `(u,v)` represent the physical cycle

\[
z_{u,v}=u+v\tau.                                               \tag{1.1}
\]

The critical `Q=1` primitive-sector probability from PR #213 is

\[
\pi_\tau(\{u,-v\}).                                           \tag{1.2}
\]

For even spin `s`, define

\[
Z_s^{(\tau)}(u,v)
=\left(\frac{u+v\tau}{|u+v\tau|}\right)^s,                    \tag{1.3}
\]

and the unnormalized continuum harmonic

\[
\boxed{
A_s(\tau)=
\sum_{(u,v)\in\mathbb Z^2_{prim}/\pm}
\pi_\tau(\{u,-v\})Z_s^{(\tau)}(u,v).}                         \tag{1.4}
\]

The rank-one probability

\[
P_1(\tau)=\sum\pi_\tau(\{u,-v\})                              \tag{1.5}
\]

is a modular scalar.  Hence the conditional harmonic

\[
H_s(\tau)=A_s(\tau)/P_1(\tau)                                 \tag{1.6}
\]

has the same spin covariance as `A_s`.

## 2. Modular covariance from physical coordinate transport

Let

\[
\gamma=\begin{pmatrix}a&b\\c&d\end{pmatrix}\in SL(2,\mathbb Z),
\qquad
\tau'=\gamma\tau=\frac{a\tau+b}{c\tau+d}.                    \tag{2.1}
\]

The standard normalized complex coordinate on the same physical torus changes by

\[
z'=rac{z}{c\tau+d}.                                        \tag{2.2}
\]

Therefore the unit direction of every physical cycle changes by the **same phase**

\[
\frac{z'}{|z'|}
=\frac{|c\tau+d|}{c\tau+d}\frac{z}{|z|}.                     \tag{2.3}
\]

Meanwhile the modular transformation only permutes the primitive sector labels, with probabilities transported exactly by the Pinson--Arguin laws already tested in PR #213.  Summing over all primitive unoriented lines therefore gives

\[
\boxed{
A_s(\gamma\tau)
=\left(\frac{|c\tau+d|}{c\tau+d}\right)^s A_s(\tau).}         \tag{2.4}
\]

Likewise

\[
\boxed{
H_s(\gamma\tau)
=\left(\frac{|c\tau+d|}{c\tau+d}\right)^s H_s(\tau).}         \tag{2.5}
\]

Thus `A_s` and `H_s` are weight-zero **spin-`s` modular-covariant functions**: their magnitude is modular invariant, while their phase follows the physical frame rotation.

### Generator checks

For `T: tau -> tau+1`, `c=0,d=1`, so the phase is one.  The sector relabelling alone leaves the harmonic unchanged.

For `S: tau -> -1/tau`,

\[
A_s(-1/\tau)
=\left(\frac{|\tau|}{\tau}\right)^sA_s(\tau).                 \tag{2.6}
\]

At `tau=i`, the spin-four phase equals one, consistent with the square fixed point.

## 3. Automorphism selection rules

Equation (2.4) immediately gives zeroes at moduli with a nontrivial torus automorphism whose frame phase is not unity in spin `s`.

### Hexagonal fixed point

At

\[
\tau_\hexagon=e^{i\pi/3}=\frac12+i\frac{\sqrt3}{2},           \tag{3.1}
\]

the torus has a 60-degree automorphism.  A spin-four harmonic acquires

\[
e^{i4\pi/3}\ne1.                                             \tag{3.2}
\]

Since the modulus is fixed, covariance forces

\[
\boxed{A_4(\tau_\hexagon)=H_4(\tau_\hexagon)=0.}              \tag{3.3}
\]

This is an exact continuum selection rule, not a numerical cancellation.

Similarly spin eight also fails the 60-degree invariance condition and vanishes at the hexagonal fixed point.

### Square fixed point

At `tau=i`, the nontrivial automorphism is a 90-degree rotation.  Spin four is invariant, so `A_4(i)` need not vanish.  This is exactly what the primitive-sector sum gives.

## 4. Numerical shape controls from the existing primitive-sector formula

Summing the PR #213 primitive probabilities with the phase (1.3) gives the following high-precision controls.  The quoted values were stable under expanding the primitive-pair box from coordinate cutoff 8 to 12; each individual sector probability itself uses the certified PR #213 evaluator.  The global weighted-sector tail was not separately re-certified in this continuation, so the non-fixed-point decimals should be treated as high-precision numerical controls rather than new rigorous enclosures.

| modulus | `P1` | `H4=A4/P1` |
|---|---:|---:|
| `i` | `0.38094744914033735446...` | `+0.77918131402676301260...` |
| `1/2+i` | `0.37162305143377145961...` | `+0.29680403442759164454...` |
| `1/2+i*sqrt(3)/2` | `0.36789317459719662787...` | `0` exactly by automorphism |
| `1/2+5i/6` | `0.36810613542505968165...` | `-0.08280183604383594154...` |
| `1/2+7i/8` | `0.36790953505947126272...` | `+0.02211526529074907955...` |

The corresponding unnormalized square value is

\[
A_4(i)=0.29682713399631153163484531698494418\ldots.            \tag{4.1}
\]

At `tau=i`, the continuum conditional spin-eight value is

\[
H_8(i)=0.99924170963860929291015300312572463\ldots,            \tag{4.2}
\]

whereas at the hexagonal fixed point spin eight vanishes exactly by the same automorphism argument.

The sign change of `H4` across the nearby `Re tau=1/2` shapes is therefore not an arbitrary fitted curvature: it is a prediction of the primitive-sector map with an exact zero anchored at the hexagonal automorphism point.

## 5. Why this is a stronger #585 positive control than a generic modular function

The map-resolved torus programme asks for candidate solution spaces with explicit connectivity semantics.  Here every ingredient is typed:

```text
lattice event       = ambient rank one,
map label           = primitive projective homology line,
finite readout      = e^{i4 theta_line},
continuum weights   = Pinson--Arguin primitive sector probabilities,
modular law         = exact sector permutation + frame rotation,
matching transport  = exact digital-Alexander slope equality.
```

There is no free amplitude after conditioning on rank one and no post-hoc basis expansion.

Therefore an analysis pipeline intended to compare modular-covariant spin-four subspaces should first recover this ray/shape law from a control observable with known semantics.  Failure localizes an orientation, modulus, period-basis, covariance or sector-dictionary problem before a more ambiguous continuum candidate enters.

## 6. Exact finite primal/matching transport remains separate from universality

For every finite honest torus and every occupation probability `p`, the lattice slope harmonic satisfies

\[
A_{s,4}(p;\tau)=A_{s,8}(1-p;\tau)                              \tag{6.1}
\]

by configurationwise slope preservation under complement.  Equation (6.1) is exact and off criticality.

By contrast, (1.4) is the **critical continuum** percolation prediction.  Agreement of a finite lattice sequence with (1.4) is an ordinary scaling-limit/universality question and should retain finite-size errors.

The two statements should not be conflated: exact complement transport is lattice topology; modular covariance is continuum geometry.

## 7. Boundary relative to original-U

This modular-covariant spin-four ray is deliberately a positive control, not a candidate replacement for original-U.  It has no declared mapping from the #275 physical source, normalizer and moving-root counterterm.  It can validate the modular/angular machinery that a future original-U candidate would use, but cannot satisfy that candidate's missing source-to-forward-column contract by itself.
