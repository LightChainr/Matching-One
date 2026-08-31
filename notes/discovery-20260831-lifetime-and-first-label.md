# The lifetime direction is already in complete E; first labels expose spatial information

This round moves the N425 W-squared observation into complete A/E observables,
separates its mean and spread mechanisms, and identifies how first-position
information can differ between mechanisms with identical complete clocks.
The full nine-cell A/E readout requested in Issue334 is completed on the
original archive, without new Monte Carlo, DP or path replay.

## 1. N425 has a named complete E lifetime direction

For `E(p)=1-F1(p)+F2(p)`,

```
integral E(p) dp = 1-E[W]/(N+1),  W=K2-K1.
```

Thus the mean-lifetime contrast is already a complete observer, not just a
piece hidden in the first thermal moment. At N425, first/second mean lifetimes
are18.29220 and18.67795 ranks. Their raw paired difference is
`-.38575 +/- .13114`. The common H4 normalization gives the following
full-observer integrals:

| N425 E integral | H4-normalized mean +/- original-batch SE |
|---|---:|
| Original complete paired paths | -.001014130 +/- .000344771 |
| Complete prefix-safe conditional readout | -.000833632 +/- .000348650 |

The two estimates are correlated readouts of the same source. The latter
uses the unchanged full-global gate: both orientations have reached rank1,
and every required R1 clock has the original whole-pair exact acceptance;
all other pairs remain original. It is not another independently significant
result or a fitted selection of the better estimate.

The additional conditional information changes A and E by the **same** amount
on every path, because it changes F2 alone. The E-A first-birth direction is
therefore exactly unchanged. The named endpoint decomposition is

```
H4 Delta K1 = -.37842998 +/- .20600989,
H4 Delta K2 = +.05358945 +/- .19640307.
```

Pointwise the lifetime change leans toward the first birth, but neither
endpoint contrast is resolved separately. Their correlated difference W is
more precise. We do not attribute the effect solely to an earlier first birth.

## 2. Nine checkpoint cells locate the complete contribution

The joint state is `(rank_first,rank_second)` in `{0,1,2}^2`. Each cell uses
one indicator shared by the paired orientations, so the constants in A and
E cancel **within every cell**. Full means are the sums of these contributions,
with all cross-cell covariances retained.

For the N425 safe E integral, the orientation-exchange grouped contributions
are:

| Joint rank group | H4 mean +/- SE |
|---|---:|
| 01+10 | -.000434177 +/- .000194806 |
| 12+21 | -.000136814 +/- .000183270 |
| 11 | -.000151698 +/- .000138928 |
| 00 | -.000072166 +/- .000072658 |
| 02+20 | -.000050214 +/- .000097256 |
| 22 | +.000011436 +/- .000082599 |

All nine individual-cell estimates and SEs are retained in the joint result;
this small table emphasizes the largest terms. The01/10 group carries about
52% of the total **point** contrast, not an independently precise fraction.
There is strong cancellation *inside* exchange pairs:01 and10 contribute
-.00554991 and+.00511574;12 and21 contribute+.00519074 and-.00532756.
After these cancellations, the main group sums point in the same direction.
This differs from claiming that one unpaired R1 stratum explains the global
observable or its variance.

The leading01/10 pairs contain R0, so the current global policy leaves them
original. More R1-only clock calculations cannot directly condition away this
leading contribution. That points the next physical continuation toward the
first-birth side of R0/R1 pairs and common-uniform spatial interventions,
rather than another enlargement of the solved R1 suffix family.

## 3. Lengthening and lifetime broadening reinforce, without resolving total spread

The exact pooled identity `E W^2=(E W)^2+Var(W)` gives N425 H4 contrasts

```
Delta (E W)^2 = 15.97182 +/- 5.44701,
Delta Var(W)  = 10.53423 +/- 4.39429,
Delta E W^2   = 26.50606 +/- 8.83952.
```

The approximately60/40 point allocation is not a precisely estimated ratio.
Both terms reinforce the previously reported J1-width direction. N325 does
not resolve corresponding orientation changes at this precision.

The [positive shape identity](p334-centered-birth-shape-energy.md) identifies
`q=A'/2` as the mixture of two normalized Beta densities and removes its known
binomial smoothing. For Y equally choosing K1/(N+1) or K2/(N+1), the fixed
reference shape energy is

```
R_ref = E[(Y-p_ref)^2]
      = E[(C-(N+1)*p_ref)^2+W^2/4]/(N+1)^2.
```

N425 H4 contrasts are `R_ref=3.94822e-5 +/- 2.61371e-5` and canonical
`Q_ref=3.95554e-5 +/- 2.61598e-5`. The corresponding centered ensemble
variance is `Delta Var(Y)=3.89365e-5 +/- 2.58576e-5`. All are less resolved
than the lifetime coordinate. The center and lifetime point contributions
reinforce; the lower precision of total shape does not arise from opposing
point signs. Center-shape uncertainty is the limiting term in this source.

Mean A/E curves depend only on birth marginals. The connected birth coordinate
`Cov(K1,K2)=Var(C)-Var(W)/4` is not identified by those mean curves and has
N425 contrast `-2.19403 +/- 3.79450`, unresolved here. No copula change or
continuum exponent is inferred from the W-squared signal.

## 4. Marked source reversal comes from marginal allocation, not extra memory

The exact accepted-R1 source formula decomposes as

```
A_D_integral = completion_D - marginal_debt_D - connected_debt_D,
connected_debt_D = r*Cov(K1,pi_D | accepted R1)/(N+1).
```

The collective connected debt is exactly its negative. Hence connected
past/winner coupling redistributes D/G loading but cancels from their sum in
every population and every delete-one replicate.

All four earlier source point-direction reversals remain after dropping the
connected term. Restoring it weakens their magnitudes. At N425 the direct
marginal debt's point decomposition is prevalence+.0005485904, conditional
K1−.0001092217, and winner share−.0010550887. The winner-allocation difference
outweighs the opposing prevalence contribution. This is an exact descriptive
accounting of the fixed accepted population, not a causal independence
intervention, non-Markovian memory claim or established population source sign.

## 5. First labels reveal an exact response-rank difference

The fixed five-site double star and C4-plus-inert mechanism have identical
full unmarked clocks and no initial direct absorber. Nevertheless the full
first-label conditional survival covariance has rank2 versus rank1. First
label identity explains89/504=17.6587% versus8/63=12.6984% of clock variance;
the binary direct/safe first-step information is zero in both.

The [first/last reciprocity](p334-first-last-label-reciprocity.md) shows how to
read these profiles from existing time-resolved final-birth marks. With safe
counts I_k and triggering-set counts b_v(k-1),

```
c_v(k)=I_(k-1)-c_v(k-1)-b_v(k-1),  c_v(0)=0,
Pr(T>k | first label=v)=c_v(k)/choose(d-1,k-1).
```

The site-centered recurrence is an invertible temporal transform: complete
first-label and time-resolved final-label profiles have the same response
rank. This makes the new rank interpretation computationally reusable, while
preventing it from being counted as independent evidence twice. A marginal
winner collision alone lacks the necessary time information. The constructed
rank2/rank1 example is not a claim of two versus one continuum fields.

### The same transform has now been applied to the two real prefixes

The old173-site pivotal integer arrays and full safe counts suffice. A thin
0.24-second recurrence, with no child solve, yields:

| Conditional first-label readout | A43042514269 | B43042505280 |
|---|---:|---:|
| Original direct sites | 0 | 0 |
| Fraction of clock variance explained by full first label | 3.12935% | 3.68332% |
| First numerical Gram mode trace fraction | 98.5365% | 96.2830% |
| Second numerical Gram mode trace fraction | 1.41095% | 3.64024% |
| Fraction of first-label innovation recovered by four old network roles | 35.2400% | 21.8961% |

The first two modes carry99.9475%/99.9232% of trace. This is a numerical
spectrum, not a certified rank bound. Both binary direct/safe floors are zero,
but the safe-label identities carry information, much of it within the old
port/interior/outside-core roles. A final-winner-inert label can still delay
the first-label clock by spending the first insertion without advancing the
event. These two real prefixes have different unmarked clocks; they are not
a new physical isoclock pair and do not represent population H4 inference.

## Sources and lifecycle

- Original population: e81dd59f, full births9c495ab1, conditional clocks0d1e586d.
- Root fixed-reference readout: da0080ec; complete nine-cell A/E vectors: bb79fd47.
- [Mean/variance and connected-birth result:be31a113](https://github.com/LightChainr/Matching-One/blob/be31a113ed9e7d1ed369261a1ab674f4af69062d/notes/p334-lifetime-square-mechanism.md).
- [Complete nine-cell, shape and source covariance:528793af](https://github.com/LightChainr/Matching-One/blob/528793af816e56e30a8d3a045e5bcad1ad3dca47/notes/p334-complete-shape-source-joint.md).
- [Marked source marginal/connected decomposition:9ed1e508](https://github.com/LightChainr/Matching-One/blob/9ed1e5082ac114d818da03d07e6cf2a315d75023/notes/p334-first-birth-winner-connected-coupling.md).
- [Exact first-label counterexample:d95a6045](https://github.com/LightChainr/Matching-One/blob/d95a6045325a4f89dea3322d0b6128afd1b0dc4d/notes/p334-isoclock-first-label-innovation.md); reciprocity:31c17d48.
- [Completed real-prefix first-label readout:65e4c677](https://github.com/LightChainr/Matching-One/blob/65e4c677097d089f564b9c5e2bd698489afe2cc0/notes/p334-physical-first-label-reciprocity.md).

All population effects use the same original20 paired batches per size, with
one common covariance coordinator. The exact constructed graphs and selected
physical prefixes are kept separate from that population inference. No new
remote job, GPU, validation suite, PR merge, issue closure or history rewrite
was used. The broader research goal remains open; this is a completed mechanism
increment rather than a declaration that every repository task is finished.
