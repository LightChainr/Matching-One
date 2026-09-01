# Exact p-biased Boolean noise-semigroup gate

Issue #227 asks whether wrapping and matching observables have a useful
Boolean spectral organization. This first gate makes only finite exact
statements on the existing `N=5` Gaussian and `N=10` C4 self-matching tori.

For `X_i~Bernoulli(p)`, use the orthonormal p-biased basis

\[
 \chi_S(X)=\prod_{i\in S}\frac{X_i-p}{\sqrt{p(1-p)}}.
\]

If `Y` keeps each bit of `X` with probability `rho` and otherwise resamples
it independently from `Bernoulli(p)`, centered observables obey

\[
 \operatorname{Cov}(f(X),g(Y))
 =\sum_{S\ne\varnothing}\rho^{|S|}\widehat f(S)\widehat g(S).
\]

The script uses `p=2/5`, so this is genuinely p-biased rather than a disguised
uniform Walsh calculation. It stores only the level sums

\[
 a_k(f,g)=\sum_{|S|=k}\widehat f(S)\widehat g(S),
\]

not a large individual-subset table.

## Independent exact cross-check

The Fourier computation enumerates single configurations and exact rational
p-biased moments. A separate path enumerates every noisy pair `(x,y)` and
groups it by counts of `(1,1)`, `(0,0)`, and unequal bits. It then evaluates
the exact joint weights

```text
P(1,1)=p^2+p(1-p)rho
P(0,0)=(1-p)^2+p(1-p)rho
P(1,0)=P(0,1)=p(1-p)(1-rho)
```

at `rho=0,1/4,1/2,3/4,1`. Every autocorrelation and declared cross-spectrum
agrees exactly between the two paths.

## Frozen observables

The common black field defines:

```text
orientation_difference = primal_direction0-primal_direction1
matching_odd_cross      = primal_cross-matching_cross(white complement)
```

The artifact retains autocorrelation generating functions for these and all
six Boolean indicator components. It also retains three signed cross-spectra:

1. primal direction-0 versus direction-1;
2. primal cross versus complement-matching cross;
3. orientation difference versus matching-odd cross.

The first two have nontrivial positive and negative Fourier-level sums at both
sizes. The third is exactly zero at every level on both tiny controls. This
is a useful finite orthogonality/null regression, not evidence that the two
sectors decouple asymptotically.

At N=10 the orientation-difference spectral support is exactly levels
`{2,3,4}`, while matching-odd cross has support `{1,2,3,4,5}`. These are
properties of this finite Boolean oracle, not fitted spectral exponents.

## Pivotal identity at rho=1

For every retained Boolean indicator,

\[
 C_f'(1)=\sum_S |S|\widehat f(S)^2
 =p(1-p)\sum_i\Pr(i\text{ is pivotal}).
\]

The right-hand pivotal mass is independently enumerated edge by edge in the
Boolean cube. The identity is exact for increasing primal events and for the
decreasing complement-matching events because the squared discrete derivative
is the unsigned pivotal indicator. Thus `C_f'(1)/(p(1-p))`, not the raw
derivative alone, equals total pivotal mass under this normalization.

## Boundary

The rho grid is one polynomial curve, not five independent tests. This tiny
gate does not infer a large-N spectral exponent, a continuum spectral sample,
or a Jordan/logarithmic mechanism. It also does not identify the spatial
Boolean spectrum with occupation-count Krawtchouk modes.

Reproduce with:

```bash
python3 scripts/exact_boolean_noise_semigroup.py \
  --output results/exact-boolean-noise-semigroup/oracle.json
python3 -m unittest discover -s tests -p 'test_exact_boolean_noise_semigroup.py'
```
