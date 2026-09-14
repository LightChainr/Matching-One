# Digital-Alexander interpretation of the Jacobsen semi-infinite-cylinder eigenvalue identity

2026-09-14.  Literature bridge for the fixed-width charge-free-energy calculation.

## 1. The numerical sequence is known

The charge-coexistence roots obtained independently from the safe site-frontier transfer are

\[
\lambda^0_{4,w}(p)=\lambda^0_{8,w}(1-p).                     \tag{1.1}
\]

For `w=2,...,8` the roots are

```text
0.5651977173836393
0.5888806999178529
0.5914171708531385
0.5922358232050263
0.5925073562056372
0.5926196333998958
0.5926727605746284
```

They agree to the displayed precision with Table 2 of J. L. Jacobsen, arXiv:1507.03027, where the square-site thresholds are computed on `n x infinity` bases by equating the largest eigenvalues of two topologically distinct transfer sectors.

Therefore this sequence is **not a new threshold estimator**.  The useful contribution of the present construction is a probability/topology dictionary for that eigenvalue identity.

## 2. The two sector languages

Jacobsen's Potts/FK formulation compares

- an `open` sector with a propagating FK cluster;
- a `closed` sector with a propagating dual FK cluster.

At `q=1`, their largest eigenvalues define the semi-infinite-cylinder pseudo-critical point.

The site probability formulation instead starts with finite-torus homology ranks and exact digital Alexander duality

\[
 r_4(\omega)+r_8(\omega^c)=2.                                \tag{2.1}
\]

Consequently

\[
 P_p^{4}(r=2)=P_{1-p}^{8}(r=0).                              \tag{2.2}
\]

After cutting one empty separator row, `r=0` has the same exponential rate as survival in the open cylinder without any horizontal essential occupied component.  The safe frontier kernel therefore has Perron root

\[
\lambda^0_{G,w}(p),\qquad
I^0_{G,w}(p)=-\log\lambda^0_{G,w}(p).                         \tag{2.3}
\]

The finite matching charge fugacity

\[
\theta_{w,m}(p)=\log\frac{P_2(p)}{P_0(p)}                    \tag{2.4}
\]

has thermodynamic rate

\[
\frac1m\theta_{w,m}(p)\to
\Theta_w(p)=I^0_{4,w}(p)-I^0_{8,w}(1-p).                     \tag{2.5}
\]

Hence the semi-infinite-cylinder balance criterion is exactly

\[
\boxed{\Theta_w(p)=0.}                                       \tag{2.6}
\]

This is the site/digital-Alexander form of the open/closed eigenvalue identity.

## 3. Why this dictionary matters

The FK transfer statement says two sector eigenvalues cross.  The probability statement adds three pieces of interpretation.

1. **Finite-event meaning.** The two sectors are the exponential tails of the two charged endpoint events `r=0` and `r=2` of the SAME rank observable.
2. **Complement map.** The upper charged sector is not an independent second model; digital Alexander converts it exactly to rank-zero survival of complementary matching sites.
3. **Balance without concentration.** The root can remain sharply selected by the sign of a difference of two exponentially small void-sector free energies even when the ordinary two-birth distribution has a wide rank-one plateau.

This is also why the root should not be described as equality of black and white winding-component intensities.  Those cylinder intensities are paired by alternation and can be equal while their long-gap/void free energies differ.

## 4. Perron derivative gives a new probabilistic slope observable

Normalize the positive left/right Perron eigenvectors of the safe transfer into its Doob/quasi-stationary row law.  Let

\[
\bar K^0_{G,w}(p)
\]

be the mean number of occupied sites added in one row under that conditioned safe phase.  Differentiating the Bernoulli row weights gives

\[
\boxed{
(I^0_{G,w})'(p)
=\frac{wp-\bar K^0_{G,w}(p)}{p(1-p)}.}                        \tag{4.1}
\]

Therefore at the charge root

\[
\boxed{
\Theta'_w(p_w^{ch})
=\frac{w-\bar K^0_{4,w}(p_w^{ch})
-\bar K^0_{8,w}(1-p_w^{ch})}
{p_w^{ch}(1-p_w^{ch})}.}                                     \tag{4.2}
\]

The transparent transfer oracle checks this identity internally through the Perron derivative.  This gives a direct physical interpretation of the eigenvalue-crossing slope as an occupation deficit of the two conditioned topological sectors.

## 5. Direct finite-size decomposition of the n^-4 shift

Using the repository reference `p_c=0.59274605079` only as a diagnostic,

\[
(p_c-p_w^{ch})w^4
=0.3402,0.3189,0.3093,0.3035,0.3002                         \tag{5.1}
\]

for `w=4,...,8`.  This is the same `Delta_1=4` behaviour reported by Jacobsen.

The new separation is

\[
\Theta_w(p_c)w^{17/4}
=1.1860,1.0996,1.0600,1.0358,1.0216,                         \tag{5.2}
\]

while

\[
\Theta'_w(p_w^{ch})w^{1/4}
=3.4851,3.4479,3.4262,3.4124,3.4031.                         \tag{5.3}
\]

Thus the pseudo-critical shift is visibly the ratio

\[
p_c-p_w^{ch}
\approx\frac{\Theta_w(p_c)}{\Theta'_w(p_w^{ch})},             \tag{5.4}
\]

with exponents `17/4 - 1/4 = 4`.

This motivates the separate Kac-(4,2) conjecture in `sector-odd-kac42-conjecture-20260914.md`; it is not needed for the exact Jacobsen/digital-Alexander equivalence.

## 6. Magnetic-gap cross-check

Jacobsen's CFT argument notes that open and closed sectors both determine the magnetic exponent

\[
x_m=5/48.                                                      \tag{6.1}
\]

For the square lattice, the per-row excitation cost therefore predicts

\[
wI^0_{G,w}(p_c)\to2\pi x_m=0.65449846949\ldots.             \tag{6.2}
\]

The safe site transfer gives

```text
w    G4             G8 complement
4    0.6808677452   0.6677642268
5    0.6702211440   0.6643381257
6    0.6651785592   0.6620431213
7    0.6622847725   0.6604282152
8    0.6604615126   0.6592750557
```

Both sectors approach the same magnetic gap while their difference is parametrically smaller.  This is an independent semantic check that the safe kernels are selecting the intended topological sectors.

## 7. Literature boundary

Primary comparison: J. L. Jacobsen, *Critical points of Potts and O(N) models from eigenvalue identities in periodic Temperley-Lieb algebras*, arXiv:1507.03027 / J. Phys. A 48 (2015) 454003.

That paper already contains the `n x infinity` root sequence, the open/closed eigenvalue criterion, the shared magnetic-exponent argument, and the empirical correction exponents `4,6,8,...`.  It explicitly remarks that equality of the leading magnetic exponent by itself does not derive the `n^-4` pseudo-critical shift.

Accordingly the present branch should claim only the new digital-Alexander/probability interpretation, the safe-site reference implementation, the Perron occupation-slope identity, and any separately audited CFT mechanism—not rediscovery of the threshold sequence.
