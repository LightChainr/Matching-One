# Independent finite intervention crosses the frozen mechanism threshold

The new one-million-prefix block reproduces the positive source-normal
birth-center response and crosses the predeclared practical-zero threshold.
The fixed primary mean is

```
T = 3.0852005663e-8 +/- 3.9187384066e-9 (new-batch SE)
T +/-3SE = [1.9095790443e-8, 4.2608220883e-8]
delta = 1e-8
frozen archive forecast = 3.6565e-8
```

The lower endpoint exceeds delta. The frozen rule therefore returns
`stop_complete_two_score_label_closure`. This is the planned finite-size
decision, not a new cutoff chosen after looking at this block. The archive
forecast lies inside the new3SE interval; its amplitude was not refitted.

## The two candidate explanations and what stops

The tested strong closure writes the complete conditional label mean as
`m_C(Z,u)=c_a(Z)+b_f(Z)s_f(Z,u)+b_s(Z)s_s(Z,u)`, allowing class-specific
intercepts but requiring a common pair of slopes across safe classes at
each prefix. It predicts exactly zero for the intervention below. We stop
promoting this as a complete state description for C at these sizes.

The positive discovery forecast transports to a new prefix population
under the specified finite policy. A relation limited to the first
Jacobian, such as J=B G, can still work: it does not claim completeness
of the conditional label mean. Class-dependent susceptibilities and higher
within-class structure also remain possible. We do not introduce or fit
additional descriptors to rescue the rejected full closure in this block.

## One new population, four predeclared coordinates

Each size has20 independent new batches of25000 uniform fresh prefixes.
All non00 prefixes are retained as zero contributions, never removed from
the denominator. There are36938 active00 prefixes atN325 and38876 atN425.
The two sizes use separate RNG domains. Within each size the two receiving
geometries and all24 recorded response coordinates retain paired covariance.

| N / own source and receiver | frozen C forecast (10^-8) | fresh C (10^-8), +/-SE | fresh A_ref (10^-7), +/-SE |
|---|---:|---:|---:|
|325 first|4.116|3.2314 +/-0.6151|-3.7989 +/-1.0875|
|325 second|3.233|2.8101 +/-0.8846|-3.0637 +/-1.4743|
|425 first|3.300|3.3661 +/-0.8345|-3.7114 +/-1.4527|
|425 second|3.977|2.9332 +/-0.7381|-3.4541 +/-1.2826|

Only the equal four-C average is primary. The negative A signs are named
secondary predictions, correlated with C, not four more independent tests.
The own W means remain weak. All cross receivers, K1/K2 and E/W are in
the saved readout, without selecting an attractive post-hoc coordinate.

## What was actually manipulated

At each prefix a full next-label census fixes the exact original scores
s and own second density score T_oo. Projecting T_oo off the exact source
Gram range gives phi_oo, with zero mean within every joint-safe degree
class and zero covariance with both original scores. B=max|phi_oo|.

The two arms are sampled **directly** from
`q_plus/minus(u|Z)=[1 +/- phi_oo(u)/B]/d`. They preserve each safe-class
mass, the outside-label law, both instantaneous rank/Euler distributions,
and both original score means. Labels use an exact rational common-CDF
coupling; ordinary uniform suffixes share remainder priorities.

The observed weighted contrast `B*(F_plus-F_minus)/2` equals the archived
normal derivative in expectation exactly. It has no finite-step Taylor
error. The B weight stays inside the population mean: this target is not
the unweighted difference of all plus and minus outcomes.

Eight N325 prefixes have rank1 source Grams; none were discarded. The
pre-frozen exact range projection handles them. All other active Grams
are rank2. All class-centering and score-orthogonality identities are
integer checks in the producer. No new model was selected from the output.

This intervention does not identify the origin of the unperturbed global
Matching-One anomaly or count continuum fields. Its baseline law is the
exact mixture of the two arms. The next scientific question is the map
from this finite susceptibility to a named global observer, not another
description of the same discovery residual.

## Execution and lifecycle

- Protocol/decision frozen and pushed at43079652/bc0a18c2;
  [Issue334 pre-production receipt](https://github.com/LightChainr/Matching-One/issues/334#issuecomment-5477517774).
  Exact dispatch manifest7fc6fcbb; producer513552c7; reader/driver4c533e0d.
- [Finite-source theory](p334-affine-normal-intervention.md)3532d5aa;
  [source implementation](p334-independent-normal-intervention-producer.md).
  Local technical smoke used a different seed and was excluded.
- Huawei NePnUn generatedN325 and551oUR generatedN425, each with14
  single-thread workers under measured14.5CPU/25GiB limits. Complete
  production drivers, including compilation, took12.578s and14.739s;
  these times exclude provisioning and transfer.
- Exactly500000 new prefixes per size;8 paired draws per active prefix
  and own axis. Unequal-label arms evaluated859754/923118 full suffix
  permutations respectively, each on both receivers. Equal-label pairs
  are exact-zero responses and skip duplicate tail evaluation.
- Every one of40 batches exited0. All162 raw/metadata/log artifacts were
  retrieved and committed atf1b36436 **before** the frozen scorer ran.
  [score.json](../results/p334-independent-normal-intervention/score.json)
  was committed at1164ba91. It retains every new-batch vector and factor;
  the old20 batches were not pooled into the new inference.
- No additional sampling after the primary result; the fixed block is
  complete. Full replication and broader validation remain with the user's
  team. The3SE rule is a fixed-batch diagnostic, not an anytime guarantee.

Lifecycle: independently sampled finite intervention / new random-prefix
population / own C primary, A secondary / N325,N425 paired geometries /
new dependency group separate from archive discovery / completed and
branch-delivered, not merged and not a closure of the general #334 issue.
