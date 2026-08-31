# From a positive birth covariance to an explicit safe-attachment mechanism

The new result is a microscopic distinction that the immediate rank and
Euler increment both miss: **among safe next sites of equal contact degree,
more contractible-cycle closure and less component merging are associated
with later first/completion births and a longer rank-one lifetime**. It is
visible at both production sizes without a new suffix run or fitted
high-dimensional model.

## 1. A finite-source failure of rank/Euler-only response closure

Let a next site have e occupied contacts in c old components. In R0-safe
attachments, `loop=e-c` is the number of new contractible graph cycles and
`merge=max(c-1,0)` measures component mergers. Both U,V in the comparison
preserve R0 and have the same e, hence the same Euler increment `1-e`.

The null `E[Y|Z,u]=f_Z(e(u))`, even allowing an arbitrary function for each
ordered prefix Z, forces the masked loop/future-response covariance to zero.
The actual same-degree covariance is nonzero at both sizes. The resulting
pooled projection slopes are:

| Clock | N325 +/- batch SE | N425 +/- batch SE |
|---|---:|---:|
| K1 | +.59009 +/- .04619 | +.71019 +/- .06240 |
| K2 | +.88512 +/- .05995 | +.96502 +/- .09445 |
| W=K2-K1 | +.29502 +/- .05763 | +.25483 +/- .08056 |

The near-critical F1 slopes are almost the same in the two sizes,
`-.01405 +/- .00110` and `-.01401 +/- .00125`. This is a finite-size
fingerprint worth explaining, not yet an exponent or universal amplitude.
At equal e, `delta loop=-delta merge` identically, so the reverse merger
direction is the same result, not another evidence block.

There is also an exact intervention interpretation. Within each safe degree
class A_e, preserve its mass pi_e and tilt the next-site choice by
`exp(t*pi_e*loop)`. This leaves the full immediate rank/Euler distribution
unchanged for every t, yet its derivative is exactly the measured covariance
numerator. Under the equal mixture of the two orientation-specific rules,
E(p_ref) increases by `8.299e-5 +/- 2.237e-5` and
`1.02777e-4 +/- 1.64317e-5` per unit t, while its integral decreases by
`1.58838e-5 +/- 3.11394e-6` and `1.14748e-5 +/- 3.64557e-6`.
This is a concrete Euler-invisible spatial tangent with a sign-changing
thermal response, not just another untyped attachment score. Its integral
effect is the already measured lifetime shift, not independent evidence.
The [exact perturbation and eight-coordinate response](https://github.com/LightChainr/Matching-One/blob/9ce53a5ad751c8b2aa37e5383cb8584637e19530/notes/p334-euler-invisible-next-label-tangent.md)
give the complete formula and joint errors.

The response moves both the birth center and the rank-one lifetime. A
pathwise common additive shift of K1 and K2 would leave W unchanged, so it
cannot alone reproduce this tangent. This does not exclude a more general
nonlinear clock reparameterization.

The [complete thermal curve](https://github.com/LightChainr/Matching-One/blob/7c60b8a728c909bbc87d359dc64422187bec4b01/notes/p334-euler-invisible-thermal-shape.md)
now locates the exchange. Main crossings occur at
**.616422 +/- .005387** and **.619771 +/- .004157** (local-delta batch SE;
all twenty delete-one branches retain the main root). The early positive
peaks are near .5925/.5965; stronger negative valleys lie near .6624/.6581.
The main negative areas are 5.34/3.64 times the early positive areas.

At both principal extrema, both F1 and F2 responses are negative. Early,
first birth loses more cumulative mass; later, completion loses more.
Thus E reverses sign through an exchange in **which delay dominates**,
not through early acceleration followed by late delay. The later tiny
positive tail lobes are below one batch SE and do not establish a third
physical phase. Complete integer birth histograms, all-p curves and
original-batch errors are published with that result.

[Complete contact-response result](https://github.com/LightChainr/Matching-One/blob/b9f79bfb6e1ba4177ff245f74f7b2e51c3bd2fdc/notes/p334-safe-contact-response-result.md)
provides raw covariance numerators, canonical and integrated responses,
the fixed-degree discriminator and all 87 original-batch coordinates.

## 2. What the original positive Gamma was comparing

The [three-mask analysis](https://github.com/LightChainr/Matching-One/blob/30c7ddb0b68fdcbfdb2e846bef6f3243290de4e3/notes/p334-rank-preserving-next-response.md)
shows that safe-safe position pairs carry only 6.82% +/- 1.81 pp and
7.82% +/- 2.47 pp of canonical Gamma, despite comprising 80.32%/83.15%
of all label quartets. About 78% lies in pairs comparing a rank-preserving
next site with a site changing at least one orientation's rank.

The safe-safe canonical residual is nevertheless positive:
`9.41865e-5 +/- 2.68289e-5` at N325 and
`7.79850e-5 +/- 2.69450e-5` at N425. It motivates the finer contact analysis,
without making it the principal source of all Gamma.

The mixed mask is **not** pure between-type covariance. Its exact
conditional expectation is

```
pi_safe*pi_birth * [B_safe+B_birth
                   +(mu_safe-mu_birth)(mu_safe-mu_birth)'].
```

The new exact census of every vacant label at every original checkpoint
supplies the true prefix-specific safety fractions, so the within/between
terms can be read separately rather than inferred from the mixed percentage.

The census changes that interpretation substantially. In the binary
safe/any-birth partition, the canonical between-type Gamma is only
`.00017343 +/- .00008245` and `.00005569 +/- .00002945`; within 01+10 it is
negative, `-.00028520 +/- .00007405` and `-.00021262 +/- .00003101`.
The large mixed mask therefore cannot be called an effect of merely
"a birth happens".

Resolve instead the two existing bits: does the next site change the first
orientation's rank, and does it change the second's? The four resulting
types 00/01/10/11 have exact prefix probabilities from the same census.
Their canonical between-type Gamma is
`.001297665 +/- .000047892` and `.000873098 +/- .000037316`, accounting
for about 93.94%/87.54% of the original total as point proportions.
Within 01+10 the corresponding point proportions are 96.15%/86.32%.
Thus **which orientation undergoes a birth**, not simply whether either
does, resolves most of the common-response covariance. The remaining
within-type covariance is smaller and is not set to zero.

The [final shared-covariance analysis](https://github.com/LightChainr/Matching-One/blob/56b383327c834236114be13f2a34f52688803e8c/notes/p334-trigger-contact-joint.md)
gives canonical four-type shares **93.94% +/- 2.78 pp** and
**87.54% +/- 3.36 pp**, integrated shares **89.21% +/- 6.52 pp** and
**76.99% +/- 6.34 pp**. Binary canonical between shares are only
12.56% +/- 5.87 pp and 5.58% +/- 2.87 pp. This same covariance includes the
contact-response coordinates above; the two findings are not treated as
independent blocks.

In the unmasked canonical Gamma, within-orientation terms contribute
`.00144939 +/- .00002913` and `.00101274 +/- .00002612`; the shared-label
cross-orientation corrections are only
`-.00006806 +/- .00002902` and `-.00001537 +/- .00003207`. The overall
positive direction is not chiefly supplied by the CRN-cross correction.
Those corrections can still be large and cancel across the three masks.

## 3. Actual contact architectures and the exact local mechanism

[Checkpoint contact coordinates](https://github.com/LightChainr/Matching-One/tree/959a7fa26677c416b874d272f1ba66523fb38f73/results/p334-next-label-contact-coordinates)
were attached to all 640,000 previously sampled next-label draws. About
63.44–63.66% of the observed safe R0 sites close at least one contractible
cycle across the four size/orientation combinations. Such a cycle is a
common safe attachment, not an exceptionally selected motif.

All 197 sampled R0->R2 draw records, representing 187 unique
prefix-label-orientation events, have a single old component supplying two
independent winding directions. The alternative two-component 2+2
architecture has zero records in this sample. It remains mathematically
possible; this sample count is not a theorem excluding it or 197 independent
prefix replications.

The [exact safe-role witness](https://github.com/LightChainr/Matching-One/blob/21bdb7b0e59155639452e26f3e75833234bfdaa5/notes/p334-safe-role-innovations-and-contact-mechanism.md)
constructs a genuine 5x5 torus prefix where a safe extension gives the next
first/completion probabilities `(4/16,1/16)`, while a local-loop insertion
gives `(2/16,0)`. Both responses move together, giving positive Gamma even
though the loop role has the lower response. This demonstrates why positive
cooperation need not mean that loop closure promotes birth. That example
does not fix e; the equal-degree population result supplies the stronger
Euler separation.

## Source boundaries and next scientific step

These are conditional-response and exact-contact analyses on the original
N325/N425 paired prefixes, using the completed e32a8593 auxiliary suffixes.
The own-orientation R0-safe contact comparison is distinct from the
both-orientations-safe Gamma mask and from an H4-projected contact response.
Projection slopes are not causal single-feature edits. The shared original
twenty batches per size carry all error propagation; no independent-source
claim, continuum field assignment or total-E explanation is made.

The new physical question is whether the attachment-partition coordinate
loads the **paired H4/E direction after orientation subtraction**, or mainly
describes an orientation-even common background. This requires applying the
fixed paired contrast to the same contact-source vectors, not acquiring
another batch or expanding the already solved R1 clock family.
