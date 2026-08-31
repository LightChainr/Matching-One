# N400 separates a weak transport defect from a resolved change of clock shape

The homothetic N400 experiment is complete. It contributes genuinely new
data: eight million shared-counter replicas for each of three modulus pairs,
with 400 aligned batches. All period entries are doubled from N100. The
three shapes share one N400 block; that block is independent of N100.
The acquisition and main six-moment readout were frozen at `894b3d8`, and
the complete histograms are committed at `3e01b49`.

**Two things happen at once:** the common-coordinate defect is no longer
resolved at this precision, while the area-normalized odd source profile
changes its thermal width very clearly. These are different observations.
An unresolved defect is not evidence that a model has become true.

## 1. The new block does not resolve the N100 transport violation

Write D=Y(4i)-Y(2i), U=Y(1/2+i)-Y(2i). The frozen class allows

`U_j(p)=r_j phi'(p) D_j(phi(p)), j=A,E`,

with one arbitrary increasing endpoint-fixed map and separate amplitudes.
Its necessary cumulative cross-moment remainders, m=1,...,6, give:

| readout | N100 | N400 |
|---|---:|---:|
| common density-coordinate remainder chi-square / df | 53.91436 / 6 | 3.90086 / 6 |
| nominal Gaussian-reference p | 7.68e-10 | .69009 |
| fixed-p free-common-secant chi-square / df | 19.21808 / 3 | 1.55791 / 3 |
| fixed-p affine E4 chi-square / df | 58.26869 / 4 | 6.69561 / 4 |

The N400 first cumulative remainder is `2.45776e-5 +/- 4.21080e-5`.
The source area ratios are `r_A=-.2109552` and `r_E=-.4015635`.
The exact one-subdivision empirical positivity certificate used at N100
does not certify both N400 curves; this neither proves a sign change nor
invalidates the necessary moment condition. No N400 quantile-map existence
claim is made.

Allowing one signed amplitude between the two six-dimensional remainder
vectors gives N400/N100 `g=-.03026`, chi-square `3.79172/5`. The approximate
95% profile interval is `[-.20682,.18347]`; zero costs only delta-chi-square
`.10914`. This is a **raw remainder amplitude**, not a dimensionless field
coupling. Its fitted negative sign is not a measured reversal. Both source
and target covariance are retained. Neither p-values nor the oriented-area
sub-readout are pooled as additional evidence.

## 2. A positive observation: the odd clock broadens in the chosen z frame

Use the declared finite-size coordinate `z=N^(3/8)(p-p_ref)` and normalize
D_A by its signed full-p area. No exponent has been fitted by this choice.
The complete histograms give:

| D_A readout | N100 | N400 |
|---|---:|---:|
| signed area | .004084910 +/- .0000686 | .000526490 +/- .0000261 |
| mean z | -.422058 +/- .0140 | -.284481 +/- .0399 |
| centered second moment in z | 1.439994 +/- .0156 | 2.123381 +/- .0738 |
| signed area inside abs(z)<=1 | .358959 +/- .00802 | .255364 +/- .0267 |
| signed area inside abs(z)<=2 | .952717 +/- .000875 | .868177 +/- .00615 |

The centered second moment increases by about `.68339 +/- .07546`.
The abs(z)<=2 fraction decreases by about `.08454 +/- .00621`. These two
summaries reuse the same underlying profile and are not independent votes.
The window half-widths .5, 1 and 2 are post-reveal descriptive choices;
full/window A/E residual moments and their complete joint covariance are
saved in the scale-transport JSON.

Thus the whole normalized clock is **not yet a size-independent curve in
this z frame**. Its width in p still shrinks; it shrinks more slowly than
the chosen critical-width coordinate over these two sizes. The U profile
is much weaker and does not support the same precision. This is a finite
two-size result, not an asymptotic exponent or a rejection of eventual
critical collapse with corrections.

## 3. Scientific consequence

The useful next question is no longer whether more N100 precision makes
the already-visible defect significant. It is which parts of the complete
birth-clock profile carry the changing width, and whether canonical
binomial smoothing or the underlying rank clock creates it. The retained
K1/K2 histograms permit that decomposition without another simulation.

The ordinary scalar class `U=a D composed with phi`, with **no Jacobian**,
is a different problem. Its ordered-extrema follow-up is being reported
separately; the density result cannot settle it.

## Lifecycle and reproduction

- Observer/sector: ordinary normalized P4[A_top,E_top], and integrated C/W.
- Geometry: tau=2i,4i,1/2+i; common rational rotation; Smith pairs
  (1,100)/(5,20) at N100 and (2,200)/(10,40) at N400.
- Source/dependency: N400 seed 20260831134001, offset 267400000000;
  one aligned block, independent of the source N100 block.
- Production: 24 million geometry-pair evaluations, not 24 million
  independent counter units; three local single-threaded processes,
  about 363--429 seconds each. No Huawei job or GPU rental was used.
- Changed mechanism space: strong N100 noncollapse cannot simply be
  promoted to a resolved larger-scale defect; whole-clock width evolution
  remains directly observable.
- Not established: model recovery, a field count, a continuum identity,
  a measured sign reversal, or a two-size asymptotic law.

```sh
python scripts/etop_finite_transport_invariants.py \
  --source results/etop-n400-three-modulus \
  --output results/etop-n400-finite-transport-invariants \
  --source-commit 3e01b495b5b637b0070705e37b4137a9a0ef0d8b
python scripts/etop_n100_n400_scale_transport.py
```

These commands analyze the saved production. They do not rerun Monte Carlo.
