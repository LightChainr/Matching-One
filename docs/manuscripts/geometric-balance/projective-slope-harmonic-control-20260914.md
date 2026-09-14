# Projective homology-slope harmonics as an exact lattice/continuum control

2026-09-14.  This note does **not** introduce a new homology-flux programme.  Issue #156 already proposed primitive homology-sector characters, and merged PR #213 supplies a general-`tau` Pinson--Arguin primitive-sector continuum evaluator.  The new contribution here is the exact interface supplied by the 2026-09-14 persistent 4/8 slope theorem.

The resulting observable is a useful positive control for angular/modular pipelines because it has all three ingredients simultaneously:

1. an unambiguous finite-lattice definition;
2. an exact primal/matching complement transport;
3. an explicit critical continuum baseline from already-merged primitive-sector probabilities.

It is **not** a replacement for original-U and must not be scored as though it were the #275 observable.

## 1. Rank-one projective slope

On an honest torus, a rank-one occupied configuration has one projective rational homology line

\[
\ell=[u:v]\in\mathbb P(H_1(T^2;\mathbb Q)),
\]

represented by a primitive integer pair `(u,v)`, modulo overall sign.  The exact classification in `projective-homology-gas-20260914.md` says every essential component in that configuration has the same line.

For a declared physical period basis

\[
\omega_1=1,\qquad \omega_2=\tau,\qquad \Im\tau>0,
\]

the embedded physical vector is

\[
z_{u,v}=u+v\tau.                                               \tag{1.1}
\]

Because the line is unoriented, any even spin gives a well-defined projective harmonic.  In particular

\[
\boxed{
Z_4^{(\tau)}(u,v)
=\frac{(u+v\tau)^4}{|u+v\tau|^4}.}                             \tag{1.2}
\]

Changing `(u,v)` to `(-u,-v)` leaves (1.2) unchanged.

This is an **embedded spin-4 geometric readout**, not a claim that the rank-one sector is a particular CFT spin-4 field.  A modular basis change transports both `tau` and `(u,v)` according to the already-frozen period-basis convention; the phase of (1.2) is referred to the declared physical `omega_1` axis.

## 2. Finite-lattice observable

For graph `G` define

\[
A_{4,G}(p;\tau)
=E_p^G\left[1_{\{r=1\}}Z_4^{(\tau)}(\ell)\right].              \tag{2.1}
\]

When `P_G(r=1)>0`, define the conditional harmonic

\[
H_{4,G}(p;\tau)
=E_p^G[Z_4^{(\tau)}(\ell)\mid r=1]
=\frac{A_{4,G}}{P_G(r=1)}.                                    \tag{2.2}
\]

More generally every even projective harmonic

\[
Z_{2k}^{(\tau)}(u,v)
=\frac{(u+v\tau)^{2k}}{|u+v\tau|^{2k}}                        \tag{2.3}
\]

is legitimate.  Spin four is singled out here because it is the first nontrivial harmonic compatible with the square `C4` geometry and it directly separates axis-like from diagonal-like rank-one sectors.

On the square torus `tau=i`,

\[
Z_4(u,v)
=\frac{(u+iv)^4}{(u^2+v^2)^2}.                                \tag{2.4}
\]

Thus

- axis slopes `(1,0),(0,1)` have `Z_4=+1`;
- diagonal slopes `(1,1),(1,-1)` have `Z_4=-1`.

Reflection symmetry makes the expectation real on the square torus; the imaginary part is an exact zero control.

## 3. Exact 4/8 complement transport

The persistent digital-Alexander theorem gives, configuration by configuration in the rank-one sector,

\[
\ell_8(\omega^c)=\ell_4(\omega).                               \tag{3.1}
\]

If black NN sites have density `p`, the complementary white matching sites have density `1-p`.  Therefore for **every** even projective harmonic and every honest torus,

\[
\boxed{
A_{2k,4}(p;\tau)=A_{2k,8}(1-p;\tau).}                          \tag{3.2}
\]

Since rank-one probabilities also match under complement,

\[
\boxed{
H_{2k,4}(p;\tau)=H_{2k,8}(1-p;\tau).}                          \tag{3.3}
\]

This is exact finite-lattice transport, not universality and not an asymptotic approximation.

It is stronger than equality of a coarse directional flag: the complete projective line is preserved.

## 4. Continuum baseline from the already-merged Pinson--Arguin evaluator

Merged PR #213 fixes the convention

\[
\text{engine }(u,v)\longmapsto\text{paper sector }\{u,-v\}    \tag{4.1}
\]

and supplies the critical `Q=1` continuum probability

\[
\pi_\tau(\{a,b\})                                              \tag{4.2}
\]

for every primitive unoriented homology sector on an arbitrary complex torus.

Therefore the continuum topological spin-4 baseline is not a fitted modular function.  It is the absolutely convergent primitive-sector sum

\[
\boxed{
A_4^{\rm cont}(\tau)
=\sum_{(u,v)\in\mathbb Z^2_{\rm prim}/\pm}
\pi_\tau(\{u,-v\})
Z_4^{(\tau)}(u,v).}                                            \tag{4.3}
\]

The total rank-one continuum probability is

\[
P_1^{\rm cont}(\tau)
=\sum_{(u,v)\in\mathbb Z^2_{\rm prim}/\pm}
\pi_\tau(\{u,-v\}),                                           \tag{4.4}
\]

and the corresponding conditional harmonic is

\[
\boxed{
H_4^{\rm cont}(\tau)=A_4^{\rm cont}(\tau)/P_1^{\rm cont}(\tau).}\tag{4.5}
\]

No field identification enters (4.3)--(4.5).  They are a direct function of the published primitive wrapping-sector probabilities.

## 5. Square-torus exact finite controls and continuum target

The new standard-library census `scripts/rank1_slope_harmonic_control.py` gives at `p=1/2`:

### `L=3`

There are 162 rank-one NN configurations out of 512, with primitive slope counts

\[
(1,0):78,\quad(0,1):78,\quad(1,1):3,\quad(1,-1):3.
\]

Hence

\[
A_4^{(L=3)}=\frac{150}{512}=\frac{75}{256},                   \tag{5.1}
\]

and

\[
\boxed{H_4^{(L=3)}=\frac{25}{27}=0.9259259259\ldots}.          \tag{5.2}
\]

### `L=4`

There are 19,932 rank-one NN configurations out of 65,536, with slope counts

\[
(1,0):9406,\quad(0,1):9406,\quad(1,1):560,\quad(1,-1):560.
\]

Thus

\[
A_4^{(L=4)}=\frac{17692}{65536}=\frac{4423}{16384},            \tag{5.3}
\]

and

\[
\boxed{H_4^{(L=4)}=\frac{4423}{4983}
=0.8876179009\ldots}.                                         \tag{5.4}
\]

In both systems the imaginary harmonic is exactly zero.  Every observed rank-one slope is axis or diagonal, so `Z_8=1` identically at these tiny sizes.

Using the merged PR #213 Pinson--Arguin formula at `tau=i` and summing primitive sectors gives

\[
P_1^{\rm cont}(i)
=0.38094744914033735446061273329420245\ldots,                  \tag{5.5}
\]

\[
A_4^{\rm cont}(i)
=0.29682713399631153163484531698494418\ldots,                  \tag{5.6}
\]

and therefore

\[
\boxed{
H_4^{\rm cont}(i)
=0.77918131402676301260435917980210484\ldots.}                 \tag{5.7}
\]

For comparison the continuum conditional spin-eight harmonic is

\[
H_8^{\rm cont}(i)
=0.99924170963860929291015300312572463\ldots.                  \tag{5.8}
\]

The finite values are controls, not a claimed monotone convergence theorem.  The fact that `H_8=1` at `L=3,4` simply reflects the absence of longer primitive slopes in those tiny honest square quotients; the continuum sum includes them.

## 6. Why this is a useful positive control for modular/angular analysis

This observable cleanly separates four questions that are often conflated.

### Lattice semantics

The input is the exact ambient homology line of a rank-one configuration.  No source normalization, derivative, rank-one denominator or continuum operator naming is needed to define it.

### Angular calibration

The phase `Z_4=e^{i4\theta_\ell}` is a literal geometric spin-four readout of the wrapped slope.  If an angular pipeline cannot reproduce its exact lattice symmetries and the known continuum primitive-sector baseline, that pipeline should not be trusted to identify a more complicated spin-four response.

### Modular/shape dependence

At general `tau`, PR #213 already supplies the sector probabilities with the correct period-basis transport.  Equation (4.3) therefore gives a parameter-free shape prediction for the slope harmonic once the physical embedding is declared.

### Primal/matching transport

Equation (3.2) supplies an exact finite complement relation unavailable for generic microscopic angular observables.  This makes the channel particularly useful for debugging orientation conventions and matching-sector maps.

## 7. Interface to #585 and #589

The map-resolved torus programme in #585 requires every candidate continuum sector to have a defensible lattice-to-map dictionary.  The rank-one primitive homology sector has exactly such a dictionary through PR #213.  The harmonic (4.3) is therefore a natural **positive control** for the projective/subspace scoring machinery before applying it to a source whose map semantics remain uncertain.

Likewise, an angular-channel design such as #589 should first demonstrate that its orientation weighting and phase convention recover the known `Z_4` transformation of this topological observable.  This does not solve the original H4/H8 mechanism question; it calibrates the geometry of the measurement.

## 8. Strict boundary relative to original-U

The topological slope harmonic is **not** original-U and should not be inserted as a surrogate forward column in #275.

A candidate mechanism for original-U still owes the declared physical source-to-observable map, normalizer and moving-root counterterm.  Success of the present positive control proves only that the homology/angular/modular plumbing is correct for an observable whose semantics are already known.

That separation is precisely why this control is useful.
