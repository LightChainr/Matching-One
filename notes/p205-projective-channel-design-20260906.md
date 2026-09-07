# Gate 2: the ranking transports, the price does not — and pricing costs more than running

**Date:** 2026-09-06
**Issues:** #595 v2, #596 Gate 2; blocks the #583 purchase decision
**Evidence type:** design scan over archived #205 histograms. Angular-resolution and sample-cost quantities only. **No Smith verdict is read from this block in any channel** (GOVERNANCE §2E), and every primary number is computed under the design null `delta = 0`.

## The objective, corrected

#595 v1 maximized `|A4| / se(A4)`. That was the wrong target and the owner's correction on #595 is adopted here in full: the decision #583 waits on is whether the angular direction can be separated from the Smith nuisance well enough to bound the induced contamination, so the design object is the joint `(A4, delta)` estimator with its full 2×2 covariance, scored by the smallest future sample multiplier for which an α-level **Fieller** set for `rho` lies inside the declared window

```text
|A8_fit/A4_fit| = |0.324936 rho / (1 +- 0.771701 rho)| < epsilon,
window(0.20) = (-0.41729, +0.41729),   window(0.05) = (-0.13754, +0.13754).
```

Both windows lie strictly inside the nuisance poles at `rho = ∓1.29584`, so pole exclusion is *implied by* the band requirement, not additional to it.

The correction has a clean form worth stating, because it is the entire difference from v1. The Fieller half-width at large `n` is `z sqrt(var(delta)/n) / |A4|`, so the figure of merit is

```text
A4^2 / var(delta)        not     A4^2 / var(A4).
```

The optimal linear combination is therefore a matched filter against the **offset** covariance, `v* ∝ S^-1 a` — a different filter from v1's, which whitened by the `A4` covariance. The code uses the exact Fieller form rather than this asymptotic shortcut; the shortcut understates the requirement whenever `A4` is weakly resolved, which is the entire regime this block sits in.

## Full-sample result

Cost is the sample multiplier on the archived 10M/node, band 0.20, each readout at its own best `p` over the predeclared grid:

```text
N=325   baseline M @ p_ref = x66.92          N=425   baseline M @ p_ref = x42.43
  combination  x  4.07   A4 snr 4.26           combination  x 11.40   A4 snr 2.26
  Sp           x  5.44   A4 snr 3.54           Sp           x 12.73   A4 snr 2.18
  S            x 18.59   A4 snr 1.91           S            x 19.30   A4 snr 1.64
  M            x 25.53   A4 snr 1.62           M            x 35.11   A4 snr 1.20
  Dp           x 91.59   A4 snr 0.90           Dp           x 99.34   A4 snr 0.73
```

**The ranking is identical at both sizes.** N=325 and N=425 use independent seeds and different geometries — they are independent by the frozen #205 contract — so this is a replication of the ordering, not a re-reading of one sample. (I looked at it *after* the fold cross-fit came back ambiguous; it is a check, not a selection.)

Against the historical choice (`M` at the frozen `p_ref`), the best readout is a **16.4× saving at N=325 and 3.7× at N=425**. So the channel work bought something real: the status quo would need ~67× the archived block to calibrate the nuisance; a chosen readout needs a few.

## Where it fails: the price is not estimable

Nested three-fold cross-fitting, choose on two folds and score on the held-out third:

```text
raw          all rounds beat baseline: False   fold-to-fold cost spread up to x527
stabilized   all rounds beat baseline: True    fold-to-fold cost spread up to x73
```

The two variants give different Gate-2 branches, so the script refuses a branch letter and returns `GATE2_NOT_ESTABLISHED`. That refusal is the honest output, and the reason is worth naming precisely.

**The validation statistic is itself a weak-denominator ratio.** The cost is proportional to `var(delta)/A4^2`, and on a 33-batch fold the `A4^2` denominator is barely resolved — the same disease this whole exercise exists to avoid, reappearing one level up. The stabilized variant shares the full-sample `A4` across discovery and validation and cross-fits only the covariance; it was introduced after seeing the raw variant's instability and is labelled as such in the code.

Even stabilized, the fold costs swing by x73, and two N=425 rounds return multipliers *below* 1 — i.e. "the experiment is already done", which it is not. That is diagnosable: `v* = S^-1 a` on 33 batches inverts a covariance estimated from 33 samples, so the matched filter can land on a poorly-estimated direction with spuriously small apparent `var(delta)`. **The combination is the unstable component**; the single channels do not have this problem, and the eigenvalue floor I used (`1e-8` of the largest) is far too permissive for a 33-sample covariance. A shrunk or regularized covariance is needed before the combination can be quoted.

## The finding that decides the gate: there is no cheap pilot

The obvious response to an unpriced experiment is to buy a small pilot and price it. That does not work here, and the reason is structural.

Both the price and the answer depend on resolving the same quantity, `A4`. Since cost `∝ 1/A4^2`, knowing the cost to ±10% needs `A4` at about 20σ:

```text
N=325, Sp:   running the band-0.20 calibration needs   x5.4
             pricing it to +-10% needs                 x31.9

N=425, Sp:   running                                   x12.7
             pricing                                   x84.5
```

**Pricing the experiment costs five to seven times more than running it.** So the choice is not "pilot then decide"; it is buy or don't.

The conservative purchase multiplier makes the same point from the other side. Taking the lower confidence edge of `|A4|`:

```text
                    point    1-sigma conservative    2-sigma conservative
N=325 combination   x 4.07          x  6.96                x   14.50
N=325 Sp            x 5.44          x 10.56                x   28.76
N=425 combination   x11.40          x 36.64                x  851.41
N=425 Sp            x12.73          x 43.59                x 1944.27
```

At N=425 the 2σ column diverges because `A4` there is only 2.2σ from zero, so its 2σ lower edge is nearly zero. **A purchase order cannot be written against that.**

## What this means for #583 and #589

Read against #596's Gate-2 branches:

- **Not branch B.** The improvement transports — every stabilized fold beats the baseline, and more convincingly the ranking replicates across two independent sizes.
- **Not branch A.** A channel cannot be "used prospectively" at a cost that is uncertain by two orders of magnitude at the conservative end.
- **Effectively branch C, with the reason sharpened.** It is not that the best honest channel leaves the nuisance unbounded; it is that **the archived block can rank readouts but cannot price them, and no affordable measurement inside this block can fix that.**

So the recommendation to #583 is unchanged in direction and firmer in basis: **do not buy N=650 as a pure `A8/A4` experiment.** The live options remain #589's 3 and 4 — the all-cyclic N=2210 design or a crossed design — neither of which needs this calibration at all. If the calibration route is taken anyway, it must be bought at a declared multiplier in the `Sp` channel with the uncertainty stated, not priced first.

`Sp` rather than the combination: `Sp` is within 10–35% of the combination's cost at both sizes, replicates its ranking across sizes, and does not depend on a 4×4 inverse covariance estimated from a third of the block. Its `p` optimum is also broad (8/21 and 11/21 grid points within 2× of the minimum) where the combination's is narrow at N=325 (6/21).

## Open, and deliberately not fixed here

- The `p` optimum sits **on the predeclared grid boundary** for several readouts, so these costs are lower bounds. A wider grid is exploratory only and the script refuses it a verdict; a finer/wider scan is specified in **#600** rather than done here.
- The combination needs a shrunk covariance. Also #600.
- `R''` and higher tail derivatives are unscanned. Also #600.

## Reproduction

```bash
D=results/server-20260829/P205-norm5-conjugate-coalescence/raw
python3 scripts/p205_projective_channel_design.py \
  --pair 325:A:$D/n325_C_A_10m.hist.csv:$D/n325_C_A_10m.moments.csv:$D/n325_C_A_10m.metadata.json \
  --pair 325:B:$D/n325_C_B_10m.hist.csv:$D/n325_C_B_10m.moments.csv:$D/n325_C_B_10m.metadata.json \
  --pair 425:A:$D/n425_C_A_10m.hist.csv:$D/n425_C_A_10m.moments.csv:$D/n425_C_A_10m.metadata.json \
  --pair 425:B:$D/n425_C_B_10m.hist.csv:$D/n425_C_B_10m.moments.csv:$D/n425_C_B_10m.metadata.json
```

~30 min. Tests: `python3 -m unittest tests.test_p205_projective_channel_design` (13).
