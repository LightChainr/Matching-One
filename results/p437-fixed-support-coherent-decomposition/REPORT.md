# The localized energy mostly resides above degree five

Archive-only decomposition of `386db0a`: **no new samples, support changes,
or noise streams**. The previous 14.97-SE B_S signal primarily measures
outside-dependent high-degree fluctuations, not one fixed degree-five
coefficient.

## Fourier interpretation

For the normalized derivative `D_S=2^-5 Delta_S` on the fixed five bonds,

`mu=E[D_SF]=Fhat(S)`.

This is exactly one degree-five Fourier coefficient. In contrast,

`B_S=|mu|² + sum_{T strictly contains S}|Fhat(T)|²`.

The second term is `Var(D_SF)` over the other 219 bonds and contains only
degree>=6 supports. The signed integer classes saved in the existing archive
identify both quantities without new simulation.

## Complex coefficient and covariance

`mu=(-7.109375e-5) + i(1.353165e-6)`.

Real/imaginary standard errors are `1.134741e-5` and `4.369614e-6`.
The full covariance of the complex mean is

```
[[ 1.2876362913e-10, -2.4690794464e-12],
 [-2.4690794464e-12,  1.9093523122e-11]]
```

The real component is -6.265 SE from zero; the imaginary component is .310 SE
from zero. The raw phase is 178.910 +/- 3.516 degrees (delta-method SE).

## Exact phase, not an inference from nonsignificance

Write a_i for the real degree-five coefficient of child i's Etop observer.
The lattice map `x -> -x` followed by translation `(9,0)` sends the period
lattice of child1 into child2 and preserves the five edge indices as a set.
The product measure and Etop are preserved, so **a1=a2 exactly**. Therefore

`Im mu=0`, while `Re mu=(a0-a1)/3` is unrestricted by this relation.

The finite certificate checks every translation and physical D4 map between
the three period lattices, also with geometric dual complement. Its only
nontrivial fixed-support relation is a1=a2; no anti-invariant fixed-support
map forces the real coefficient to zero. Allowed child coefficients span
`(1,0,0)` and `(0,1,1)`. Duality is geometric edge transport, **not naive
same-index bit complement**; an odd support does not by itself force zero.

This is a marginal coefficient relation, not a pointwise conjugation of the
common random stream. It was established by exact support transport, not by
the imaginary component's p value. The observed negative real sign remains
a measured result, not a symmetry prediction.

## Bias-corrected spectral energies

Let Z_b be the complex derivative mean in each of m=100 independent batches.
Use the cross-batch U-statistic

`U = sum_{b!=c} Re(conj(Z_b) Z_c)/(m(m-1))`

`  = |mean Z|² - tr Cov(mean Z)`.

It is unbiased for |mu|². The removed squared-mean bias is
`1.47857152e-10`. All energy uncertainties and their joint covariance are
computed by deleting one original batch at a time. Negative U estimates would
be retained and called unresolved, not interpreted as negative energy.

| Parameter | Estimate | Batch-jackknife SE |
|---|---:|---:|
| B_S | 3.23893229e-6 | 2.16390311e-7 |
| Coherent exact-degree5 weight | 4.90829519e-9 | 1.61618190e-9 |
| Outside-dependent degree>=6 weight | 3.23402400e-6 | 2.16409688e-7 |
| Coherent fraction | .00151541 | .00051128 |

Thus **99.8485% +/- .0511 percentage points** of the sampled energy is
assigned to the outside-dependent component. The fraction is a descriptive
ratio with jackknife uncertainty, not an exactly unbiased ratio or a claim
about all supports in the lattice.

## Stronger exact population bound from the same data

The degree-five term has multiplier `h5=9765/32768`; the residual terms have
degree at least six, so each multiplier is at least
`h6=615195/1048576=(63/32)h5`. Hence

`A_HP >= h5 |mu|² + h6 (B_S-|mu|²)`.

Using the U-statistic and full joint covariance, this new RHS parameter is

**`(1.89885057 +/- .12695997)e-6`**.

Its ratio to the preceding h5 B_S estimate is
**1.967282 +/- .000495**. This is an exact stronger population inequality
with an estimated RHS. The reported number is **not a statistically certain
lower bound**. It is an explicitly secondary reuse analysis, not a change to
the original frozen B_S primary.

## Provenance and scope

Input: the unchanged 20k block with SHA-256
`4e637b4d5d8388d3abfa2b1fca453da7ef39cd1d52dada01c3698d4a5a156063`.
Dependency group remains
`p437-N112-fixed-S5-lower-bound-fresh20k-20260831`.
The new JSON retains full complex-mean covariance, joint energy covariance,
the exact symmetry certificate, and the strengthened-bound propagation.

Reproduce:

```sh
python3 scripts/score_p437_coherent_decomposition.py \
  results/local-20260831/P437-N112-fixed-S5-20k \
  --output results/p437-fixed-support-coherent-decomposition/latest.json
```

Four focused tests cover the U-statistic identity, negative-estimate handling,
support symmetry, and exact h6/h5 ratio. No new MC, support scan, physical
field identification, PR, or merge is involved.
