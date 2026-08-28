# Reliability signatures, threshold-rank tomography, and a Gaussian shape baseline

**Status:** exact finite combinatorics plus exploratory research synthesis, 2026-08-28.

This note proposes a new finite object for organizing several strands of Matching One:

- exact self-matching polynomials;
- Newman--Ziff / threshold-rank histograms;
- `kappa3`, `kappa5`, and the universal threshold profile;
- the exact two-sublattice tangent family;
- Boolean pivotal/influence and higher interaction observables.

The central observation is elementary but powerful:

> For any monotone torus wrapping event, the Newman--Ziff activation rank is exactly the activation analogue of a coherent-system reliability signature.

This gives a discrete object before any continuum fit. It also turns the N=26 failure of the frozen `Beta(s,s)` laws into a structural obstruction rather than merely a polynomial mismatch.

The reliability-signature literature represents coherent-system reliabilities as mixtures of `k-out-of-n` reliabilities / order-statistic distributions. See Samaniego's signature framework and Marichal--Mathonet--Waldhauser (2011), *On signature-based expressions of system reliability*.

## 1. From a wrapping event to a domination vector

Let

```text
f : {0,1}^N -> {0,1}
```

be an increasing event, for example rank-2/cross wrapping on a fixed finite torus. For every occupation number `k`, define

```text
a_k = number of successful k-subsets,
d_k = a_k / C(N,k).
```

The vector `d=(d_0,...,d_N)` is the layerwise success probability, also called a domination / tail-signature type object in reliability theory.

Under iid Bernoulli occupation with probability `p`,

\[
F(p)=\mathbb E_p[f]
=\sum_{k=0}^{N} d_k {N\choose k}p^k(1-p)^{N-k}.
\]

Thus the Bernstein coefficients used throughout this repository are exactly the layerwise reliability data of the monotone Boolean event.

## 2. Newman--Ziff threshold rank is an activation signature

Take a uniformly random permutation of the N sites and add sites in that order. Let

\[
K=\min\{k:f(S_k)=1\}
\]

be the first occupation rank at which the event occurs.

Conditioned on `|S_k|=k`, `S_k` is a uniform k-subset, so

\[
\Pr(K\le k)=d_k.
\]

Define

\[
q_k=\Pr(K=k)=d_k-d_{k-1},\qquad d_{-1}=0.
\]

Then `q=(q_1,...,q_N)` is the **activation signature** of the finite wrapping event.

This is not an analogy: it is exactly the histogram object estimated by an increasing-event Newman--Ziff threshold-rank run, up to the repository's declared `K_plus/K_minus` orientation and complement conventions.

The reliability curve can be reconstructed as the order-statistic mixture

\[
\boxed{
F(p)=\sum_{r=1}^{N}q_r I_p(r,N-r+1),
}
\]

where `I_p` is the regularized incomplete beta function. Equivalently, the continuum-looking beta mixtures already used by the threshold-rank scorer have a standard finite reliability-signature interpretation.

### Self-dual/self-matching symmetry

If a centered event satisfies

\[
F(p)+F(1-p)=1,
\]

then coefficient uniqueness gives

\[
d_k+d_{N-k}=1,
\]

and therefore

\[
\boxed{q_k=q_{N+1-k}.}
\]

So exact complement antisymmetry becomes mirror symmetry of the activation signature.

This is a useful exact regression target for every self-matching finite quotient.

## 3. Exact `Beta(s,s)` means a majority-core signature

The regularized beta CDF has the binomial-tail identity

\[
I_p(s,s)
=\Pr\{\mathrm{Bin}(2s-1,p)\ge s\}.
\]

Therefore an exact law

\[
F(p)=I_p(s,s)
\]

has a very concrete Boolean reliability interpretation:

> its layer-count / activation signature is exactly the same as majority on a fixed core of `2s-1` relevant variables, with the remaining `N-(2s-1)` variables acting as dummies.

This statement concerns the reliability signature. It does **not** imply the wrapping Boolean function is isomorphic to a literal majority gate; different coherent systems can share a reliability signature.

Degree elevation from the `2s-1` core to N sites gives the exact successful-layer counts

\[
\boxed{
a_k^{\rm Beta(s,s)}
=\sum_{j=s}^{2s-1}
{2s-1\choose j}{N-(2s-1)\choose k-j}.
}
\]

In particular, at the first possible successful layer,

\[
\boxed{a_s^{\rm Beta(s,s)}={2s-1\choose s}.}
\]

This is a one-coefficient obstruction: a finite geometry cannot have exact `Beta(s,s)` reliability unless the number of minimal-size successful configurations equals this majority-core count.

### N=10

The exact self-matching N=10 quotient has `F=Beta(3,3)`. The majority-signature condition requires

```text
s = 3
2s-1 = 5
a_3 = C(5,3) = 10.
```

This is a necessary signature identity behind the exact Beta law. It should be checked explicitly from the committed N=10 enumeration whenever the N=10 Beta result is used as a control.

### N=26

PR #152 exhaustively enumerates the primitive `(5,1)` N=26 C4 self-matching quotient. Its geometry-only frozen hypothesis was `Beta(5,5)`.

For an exact Beta(5,5) law, the majority-core signature requires

\[
a_5={9\choose5}=126.
\]

The exact target reveal instead has

```text
a_5 = 78
```

for the canonical `F=(1+M)/2` direction, with 156 raw five-site `either`-wrapping masks split 78/78 between quotient directions.

Thus the earliest N=26 coefficient discrepancy has a direct structural meaning:

\[
\boxed{78\ne126.}
\]

The geometry has only 78 critical five-site success masks where a Beta(5,5) / nine-variable-majority signature requires 126.

The `Beta(7,7)` antipodal-majority hypothesis has an even more basic obstruction: it has no successful layer before occupation 7, while the geometry has successful occupation-5 configurations. Its full exact score remains useful provenance, but the signature exposes why the first discrepancy must already occur at k=5.

This interpretation should replace attempts to rescue the N=26 result with a generalized beta fit.

## 4. The entire threshold curve is a signature, not just a root

The vector `q_k` is richer than the finite root or a few derivatives. It supports exact finite comparisons between geometries without referring to disputed infinite-volume `p_c` digits.

Useful operations include:

1. mirror-symmetry tests under self-matching;
2. comparison of signatures by stochastic order;
3. decomposition into `k-out-of-N` basis curves;
4. moments/cumulants of the activation rank;
5. comparison across torus modulus, orientation, or microscopic model;
6. information-optimal compression of a full threshold-rank histogram.

The important methodological point is that `q` is computed **before** fitting a continuum profile. If two proposed continuum mechanisms imply incompatible finite signatures, the finite combinatorics can reject them without a large-N exponent fit.

## 5. `kappa3` as a tilted-signature variance defect

The signature gives a new exact interpretation of the metric-free shape invariant.

Write

\[
F(p)=\sum_r q_r I_p(r,N-r+1)
\]

and define the beta-density value at the self-dual center

\[
h_r
=\left.\frac{d}{dp}I_p(r,N-r+1)\right|_{p=1/2}
=\frac{N{N-1\choose r-1}}{2^{N-1}}.
\]

Let

\[
A=F'(1/2)=\sum_r q_r h_r
\]

and define the slope-tilted activation signature

\[
\pi_r=\frac{q_r h_r}{A}.
\]

For the centered matching function

\[
M(p)=2F(p)-1,
\]

a direct derivative of the beta densities gives

\[
\frac{h_r''(1/2)}{h_r(1/2)}
=4\left[(2r-N-1)^2-(N-1)\right].
\]

Hence

\[
\boxed{
\kappa_3
=\frac{M'''(1/2)}{M'(1/2)^3}
=\frac{\mathbb E_{\pi}\left[(2K-N-1)^2\right]-(N-1)}{A^2}.
}
\]

So `kappa3` measures a **critical activation-rank variance defect** under a natural slope/pivotal tilt. This makes the bridge to Russo/pivotal mass (#100) and Potts thermal cumulants (#54) much more explicit.

The first derivative `A` is the critical density of the activation signature; for increasing events Russo's formula identifies the corresponding derivative with total pivotal mass. The third derivative measures how the activation-rank distribution departs from the binomial/order-statistic curvature at the center.

## 6. A principled Gaussian/majority baseline for `kappa3` and `kappa5`

The project currently treats the possible limit

\[
\kappa_3\approx-5/3
\]

as an aggressive candidate. A more informative null model than zero is supplied by the majority/Beta family itself.

For

\[
M_s(p)=2I_p(s,s)-1,
\]

let `d_s` be the `Beta(s,s)` density at 1/2. Exact differentiation gives

\[
\kappa_3(s)=-\frac{2(s-1)}{d_s^2},
\]

\[
\kappa_5(s)=\frac{12(s-1)(s-2)}{d_s^4}.
\]

Examples:

```text
s=3: kappa3 = -256/225                = -1.1377777778...
     kappa5 = 32768/16875             =  1.9418074074...

s=5: kappa3 = -131072/99225           = -1.3209574200...
     kappa5 = 4294967296/1093955625    =  3.9260891373...

s=7: kappa3 = -4194304/3006003        = -1.3953093194...
     kappa5 = 43980465111040/9036054036009
            =  4.8672202419...
```

As `s -> infinity`, the centered majority threshold has the usual Gaussian / error-function scaling limit. Since

\[
d_s\sim 2\sqrt{s/\pi},
\]

one obtains the parameter-free shape baseline

\[
\boxed{\kappa_3^{\rm Gaussian}=-\pi/2}
\]

and

\[
\boxed{\kappa_5^{\rm Gaussian}=3\pi^2/4}.
\]

Numerically,

```text
-pi/2       = -1.5707963267948966...
3 pi^2 / 4  =  7.4022033008170185...
```

This suggests a sharper scientific question for #16/#25/#54:

> Does the universal matching-threshold signature converge to the Gaussian/majority shape, or is its limiting profile genuinely non-Gaussian?

If a controlled multi-lattice limit excludes `-pi/2` and remains compatible with `-5/3`, that is much stronger evidence than noticing that `-1.67` is close to a simple rational number. `kappa5` should be used jointly: the Gaussian null predicts both numbers at once.

## 7. Anisotropic tangents are Boolean Fourier / Banzhaf tomography

The exact self-matching tangent in PR #148 supplies another connection.

For an event `f` and a sign assignment `sigma_i=+1/-1`, consider the biased product measure

\[
p_i=1/2+\sigma_i\lambda.
\]

Write the Walsh expansion at the unbiased measure as

\[
f(x)=\sum_S \widehat f(S)\chi_S(x).
\]

Then product independence gives

\[
\boxed{
\mathbb E_{\lambda}[f]
=\sum_S \widehat f(S)(2\lambda)^{|S|}\prod_{i\in S}\sigma_i.
}
\]

Therefore the coefficient of `lambda^r` is

\[
2^r\sum_{|S|=r}\sigma_S\widehat f(S).
\]

These signed level sums are Boolean Fourier interaction moments; in reliability/game-theory language they are closely related to Banzhaf influence and higher Banzhaf interaction indices. The first derivative is a signed pivotal/Birnbaum importance sum; higher odd derivatives probe coordinated multi-site interactions.

For the N=10 self-matching tangent, PR #148 finds exactly

\[
R_-(\lambda)=\frac54\lambda-4\lambda^5.
\]

Thus, for its declared even/odd sublattice sign convention, the signed Fourier level sums are

```text
level 1:  +5/8
level 3:   0
level 5:  -1/8
```

with the remaining allowed odd levels absent in that degree-5 polynomial slice.

This is a much more microscopic statement than the isotropic Beta(3,3) law. The two facts are compatible: a reliability signature fixes layer counts, not the Boolean function up to isomorphism. Two coherent systems can have the same isotropic signature and different anisotropic Fourier/Banzhaf spectra.

## 8. Cheap next exact experiment after the N=26 reveal

Do **not** fit another beta family to N=26.

Instead, reuse the exact N=26 engine to enumerate the bivariate occupation table

```text
(k_even, k_odd) -> successful configuration counts
```

for the C4 self-matching quotient. This costs essentially the same 2^26 traversal already demonstrated in PR #152.

From this one table reconstruct, exactly:

1. the isotropic domination vector `d_k` and activation signature `q_k`;
2. the two-sublattice tangent `R(t,lambda)`;
3. the signed odd interaction spectrum in `lambda` at `t=0`;
4. mixed thermal/anisotropy derivatives;
5. channel-by-channel versions before taking the exact matching-difference identities.

### Pre-analysis contract

The following are exact structural expectations and may be checked without model fitting:

- exchange/complement parity fixes the allowed parity of Taylor monomials;
- the isotropic signature must reproduce the already committed N=26 polynomial;
- self-matching implies mirror activation-signature symmetry;
- the Fourier/Banzhaf reconstruction from the bivariate table must agree with direct finite differences.

Do **not** assume the N=10 accidental sparsity `lambda + lambda^5` persists at N=26. Whether level 3 or higher interaction sums appear is the scientific result.

A useful outcome would be a size sequence of signed interaction spectra. It would tell us whether matching tangent-space simplicity is an actual RG structure or a tiny-quotient accident.

## 9. Consequences for the research program

### Exact self-matching work

The N=26 result should be reframed from “Beta family fails” to:

- exact self-duality gives a symmetric activation signature;
- N=10 happens to have the majority-core signature Beta(3,3);
- N=26 fails the majority-core signature at the first successful layer;
- the next exact object to compare is the full signature and anisotropic interaction spectrum, not another beta fit.

### Universal scaling function

The threshold-rank pipeline already estimates the finite activation signature. Preserve it as a primary scientific artifact rather than only a way to reconstruct roots and derivatives.

### `kappa3`

Add `-pi/2` and `3pi^2/4` as a coupled Gaussian/majority baseline for `kappa3/kappa5`. This is a mechanism-level competitor to a non-Gaussian rational candidate such as `-5/3`.

### Pivotal and operator bridges

The slope-tilted signature formula and anisotropic Fourier moments give a finite Boolean bridge to:

- Russo/pivotal measures (#100);
- four-arm anisotropy (#121);
- operator-mixing matrices (#125);
- Potts homology thermal cumulants (#54).

This bridge is attractive because every object is exactly defined before invoking an LCFT operator name.

## 10. Claim boundary

Exact statements in this note:

- threshold-rank activation probabilities are determined by the layerwise domination vector;
- the reliability curve is the activation-signature mixture of order-statistic / beta CDFs;
- self-dual complement symmetry gives mirror signature symmetry;
- exact `Beta(s,s)` reliability has the same domination/signature vector as majority on `2s-1` relevant variables plus dummies;
- the minimum-layer count obstruction `a_s=C(2s-1,s)`;
- the tilted-signature formula for `kappa3`;
- the Beta-family formulas and Gaussian limits for `kappa3/kappa5`;
- the product-measure/Fourier expansion underlying anisotropic tangent tomography.

Empirical inputs from open PRs, pending merge:

- PR #152: the exact N=26 `a_5=78` result and full polynomial;
- PR #148: the exact N=10 two-sublattice tangent polynomial.

Conjectural / research-level statements:

- that a limiting activation signature has a simple continuum characterization;
- that its non-Gaussianity is controlled by a particular LCFT module;
- that N=10 Fourier sparsity persists to larger self-matching quotients;
- that `kappa3=-5/3` is exact.

## References

- F. J. Samaniego, system-signature framework for coherent systems.
- J.-L. Marichal, P. Mathonet, T. Waldhauser, *On signature-based expressions of system reliability*, Journal of Multivariate Analysis 102 (2011) 1410--1416, arXiv:1010.0162.
- M. E. J. Newman and R. M. Ziff, *Efficient Monte Carlo algorithm and high-precision results for percolation*, Phys. Rev. Lett. 85 (2000), and the follow-up Newman--Ziff microcanonical algorithm paper.
- J.-L. Marichal et al., work on Banzhaf influence/interaction indices for pseudo-Boolean functions and their least-squares / difference interpretations.
- Existing Matching One Issues #16, #25, #54, #100, #118, #121, #125 and PRs #148/#152.
