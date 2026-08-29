# Norm-five chiral Hecke phase discriminator

This is an interest-driven mechanism test built on the exact character basis
of #239 and the selection rule of #244.  The selection rule removes unmarked
scalar detail responses, but it leaves a charged matrix element:

\[
S_\chi=\sum_{j,k}\zeta_5^k(X_{j,k}-p),\qquad
O_{\bar\chi}=\sum_{j,k}\zeta_5^{-k}o_{j,k},
\]

\[
R_m=\frac{E[O_{\bar\chi}S_\chi]}{p(1-p)}.
\]

Here \(k\in\mathbb Z/5\) is the exact fiber coordinate for the Gaussian
cover \(m=2\pm i\).  For production, \(o_{j,k}\) should be the existing
complement-odd local landing/pivotal \(H_4\) mark, transported to every fiber.
The tiny oracle below uses the pivotal derivative of matching-odd cross wrap,
because the five-site torus has no useful annular landing radius.

## Frozen chiral prediction

For a spin-\(s\) Hecke eigenfield the handed ratio is

\[
\frac{A_s(2+i)}{A_s(2-i)}
=\left(\frac{2+i}{2-i}\right)^s
=\left(\frac{3+4i}{5}\right)^s.
\]

The three targets are:

| hypothesis | exact ratio | phase |
|---|---|---:|
| \(H_4\) | \((-527-336i)/625\) | \(-147.480^\circ\) |
| \(H_8\) | \((164833+354144i)/390625\) | \(65.041^\circ\) |
| \(H_{12}\) | \((32125393-242017776i)/244140625\) | \(-82.439^\circ\) |

All have unit modulus and are exactly distinct.  The smallest angular
separation is about \(65.04^\circ\), so a single handed pair is mathematically
enough to distinguish the three candidates once the reflection/character
transport is frozen.  It is not enough if each hand is allowed an unrelated
complex normalization.  The proposed paired measurement avoids that freedom:
reflection maps \((x,y)\mapsto(x,-y)\), preserves the cyclic deck label, and
anti-linearly transports \(\chi_+\) to \(\bar\chi_-\); the common complex
normalization then cancels in \(R_{2+i}/R_{2-i}\).

This is a sharp discriminator under the explicit mechanism hypothesis that
one spin eigenfield dominates the measured charged row.  A mixture can land
between the three phase targets; that would be a result rather than a protocol
failure.

## Tiny exact channel oracle

The parent of order one and its two norm-five children are exhausted at
\(p=2/5\).  For both \(2+i\) and \(2-i\), all \(2^5\) configurations give

\[
E[O_{\bar\chi}S_\chi]/[p(1-p)]=-46/25.
\]

The value is represented exactly in \(\mathbb Q(\zeta_5)\), agrees with an
independent symbolic derivative of the inhomogeneous Bernoulli polynomial,
and obeys reflection conjugacy.  Its reality is a special symmetry of this
minimal torus; it validates the charged channel and transport but does not
select \(H_4,H_8\), or \(H_{12}\).

Reproduce the prediction and exact oracle with:

```bash
python3 scripts/norm5_chiral_hecke_phase.py \
  --output predictions/norm5_chiral_hecke_phase_20260829.json
```

A production implementation should accumulate one shared-randomness batch
vector containing `Re/Im R_(2+i)` and `Re/Im R_(2-i)`, retain its full 4x4
covariance, and score the complex handed ratio against the three frozen points.
That is a measurement prescription, not an instruction to multiply marginal
z scores.

