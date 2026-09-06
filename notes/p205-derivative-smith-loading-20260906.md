# The #205 block cannot calibrate the Smith nuisance — in any channel

**Date:** 2026-09-06
**Issues:** #589 (option 2), #591; blocks the interpretation of #583
**Evidence type:** re-view of the frozen #205 evidence block in additional observable channels. **Zero new sampling, zero refit, no harmonic selected.** This is the same evidence viewed differently, not a new independent vote (GOVERNANCE §2E), and it cannot move #205's H4 verdict in either direction.

## What was asked

#589 found that the N=650 three-orientation design changes the lattice angle and the
Smith/quotient class at the same time, so its third fitted coefficient is not a clean
`A8`. #591 made that exact: a unit offset `delta` on the one noncyclic row lands as

```text
Delta C  = +0.4593316585 delta      (both families)
Delta A4 = +0.7717014729 delta      (square)   /  -0.7717014729 delta (rectangle)
Delta A8 = +0.3249364887 delta      (both families)
```

recomputed here from scratch in exact rational arithmetic (`matches_issue_591: true`).

#589 listed four ways to recover a pure angular reading and recommended trying the
cheapest first: **option 2**, calibrate the noncyclic loading from #205's archived
histograms in the derivative channel that #583 is actually motivated by, before
commissioning anything larger.

This note reports that attempt. **Option 2 fails on the existing data.**

## The measurement

For each size the frozen #205 integer H4 null is, exactly, the two-cyclic-row `cos(4θ)`
interpolation — verified, not assumed:

```text
N=325   frozen [5, -11, 6]     ==  M_C - (11/5) M_A + (6/5) M_B      defect 0
N=425   frozen [20, 13, -33]   ==  M_C + (13/20) M_A - (33/20) M_B   defect 0
```

so the published `residual` **is** `delta`, the noncyclic offset. Alongside it, the same
three correlated numbers give the angular amplitude from the two *cyclic* rows only,

```text
A4 = (x_A - x_B) / (cos4_A - cos4_B),
```

which deliberately gives the noncyclic row weight zero: the denominator must not absorb
the offset the numerator is measuring. The exported quantity is the dimensionless

```text
rho = delta / A4
```

with a delta-method error that keeps the `delta`/`A4` cross term (the two functionals
share the A and B rows, and the common priority field correlates all three; the measured
correlation is |0.76|–|0.90|, so dropping the cross term is not an option).

Same frozen `p_ref = 0.59274605079`, same aligned delete-one jackknife, same guards.

## Control: this is the published pipeline

The `M` channel run through this code path reproduces #205's archived `analysis/score.json`:

```text
H4   chi2_2 = 2.722843  p = 0.25629   z = +1.1793974 / -1.154064
H12  chi2_2 = 0.711888  p = 0.70051   z = -0.032126 / +0.8431228
H8   chi2_2 = 0.469946  p = 0.79059   z = +0.4996773 / -0.469328
```

A second, free control: `obs()` defines `D = (rg-rh)/2` and `M = rg-rh`, so `D` is exactly
`M/2` and `rho` — being a ratio — must be bit-identical. It is. That is a scale-invariance
check on the ratio machinery, and it also means `D` carries no information beyond `M`.

## Result

```text
  ch     N   delta_z      A4_z         rho          +-  |A8/A4|2s  pole?  x samples
  M    325 1.1793974 1.0236792 1.616698183 0.737686559    0.72483    YES    12.5002
  M    425 -1.154064 1.0826992 -1.23524991 0.764565799    0.79261    YES    13.4277
  S    325 0.7736940 0.5921551 1.718020674 1.386120087    0.59187    YES    44.1341
  S    425  1.286229 -1.230702 -1.27645747 0.702033425    0.81511    YES    11.3211
  D    325 1.1793974 1.0236792 1.616698183 0.737686559    0.72483    YES    12.5002
  D    425 -1.154064 1.0826992 -1.23524991 0.764565799    0.79261    YES    13.4277
  Sp   325 1.8260857 2.1824433 1.269358388 0.349250922      1.233    YES    2.80187
  Sp   425 -0.833353  1.602811 -0.68650016 0.529980996     1.6319    YES    6.45199
  Dp   325 0.4217717 0.1165169 5.103477982 33.38451791     0.4288    YES    25601.4
  Dp   425 -0.667873 0.0581840 -15.7677861 253.5042081    0.42211    YES 1.47619e+6
```

`x samples` is the sample-count multiple that would shrink the 2σ window into the
`|A8/A4| < 0.20` band. Pre-registered bands (written before any derivative number was
computed; the only #205 numbers already seen at that moment were the published `M`
scores, which is why `M` and its exact rescaling `D` are excluded from the verdict):

```text
worst induced |A8/A4| at 2 sigma  <  0.05   ->  N=650's third coefficient readable as A8
                                  <  0.20   ->  readable with a stated contamination band
                                  otherwise ->  needs an independent quotient control
```

Verdict, read from `Sp`/`Dp`:

```text
DERIVATIVE_SMITH_LOADING_UNCONSTRAINED__N650_NEEDS_AN_INDEPENDENT_QUOTIENT_CONTROL
```

## Read this correctly: it is a power failure, not a detection

The noncyclic offset is **not detected anywhere**. The largest `|delta_z|` across all ten
channel×size cells is 1.83 (`Sp`, N=325); every cell is consistent with zero offset. That
is the same reassurance #205 already gave, now extended to the derivative channels.

What is missing is the *denominator*. The angular amplitude `A4` is itself only
0.06σ–2.18σ from zero at these sizes. A ratio of two things that are each consistent with
zero is unconstrained, and that ratio is precisely what has to transfer to N=650.

The sharpest form of the failure is structural rather than numerical. The induced ratio

```text
A8_fitted / A4_fitted  =  0.32494 rho / (1 +- 0.77170 rho)
```

has a pole at `rho = ∓1.29584`, where a quotient response cancels the true angular
amplitude exactly and N=650's fitted `A4` is *entirely* nuisance. **Every one of the ten
cells has a 2σ interval that reaches that pole**, in at least one modulus family. So the
existing block does not merely fail to bound the contamination; it fails to exclude the
case in which the N=650 fit has no angular content to report at all.

This holds in the published `M` channel too. #205's "no gross cyclic-to-noncyclic
breakdown" was and remains a correct statement about an *offset*; it was never a bound on
the offset *relative to the angular amplitude*, and #583 needs the latter.

## What it costs to fix by this route

Assuming the true offset is zero — load-bearing, and stated as an assumption: if the
offset is real but unresolved, more samples make it *measurable*, not small — the error
shrinks as `1/sqrt(n)` and the required multiples are:

```text
channel  N     x samples for banded (0.20)   x samples for clean (0.05)
Sp      325             2.8                          25.8
Sp      425             6.5                          59.4
S       425            11.3                         104.2
M       325            12.5                         115.1
Dp      325        25,601                        235,652
Dp      425     1,476,190                     13,587,800
```

Two things follow.

1. **`Sp` is the only viable channel.** At 10M samples per node it is already the tightest
   cell; ~30M–65M per node would reach the banded window, ~260M–600M the clean one.
   `Dp` is out of reach by four to seven orders of magnitude — its `A4` is 0.06σ, so the
   derivative *difference* channel has essentially no angular leverage at N≈325–425.
2. **The cheap route is no longer clearly the cheap one.** #589 recommended option 2
   because it required no new sampling. It does require new sampling, of order 3–6× the
   existing #205 production, and only for a *banded* reading.

## Consequence for #583 and #589

- #589's option 2 is **attempted and closed** on existing data. It is not refuted as a
  method; it is unaffordable at #205's sample size, and its cost is now quantified.
- **#583 must not be interpreted as measuring `A8/A4` at N=650.** Its decision-rule
  amendment stands, and the branch "large fitted third coefficient + derivative-channel
  Smith control small" is currently unreachable — the control is not small, it is absent.
- The live options are #589's 3 and 4: the all-cyclic N=2210 design (four primitive
  directions in both families, worst amplification 0.611 vs N=650's 0.893), or a crossed
  design that estimates the quotient loading while keeping angular leverage.
- Before either, one cheap prospective question is worth asking, because it changes the
  cost by a factor of ten: **is there a channel with more angular leverage at fixed
  sample count?** `Sp` beats `M` by 4.5× and `Dp` by 9,000× on exactly this figure of
  merit, and nothing in the repository has yet optimized the channel *for* `A4`
  resolution. That is a design question, answerable from the same archived histograms,
  and it is opened as a separate item rather than folded in here.

## Reproduction

```bash
D=results/server-20260829/P205-norm5-conjugate-coalescence/raw
python3 scripts/score_p205_derivative_smith_loading.py \
  --pair 325:A:$D/n325_C_A_10m.hist.csv:$D/n325_C_A_10m.moments.csv:$D/n325_C_A_10m.metadata.json \
  --pair 325:B:$D/n325_C_B_10m.hist.csv:$D/n325_C_B_10m.moments.csv:$D/n325_C_B_10m.metadata.json \
  --pair 425:A:$D/n425_C_A_10m.hist.csv:$D/n425_C_A_10m.moments.csv:$D/n425_C_A_10m.metadata.json \
  --pair 425:B:$D/n425_C_B_10m.hist.csv:$D/n425_C_B_10m.moments.csv:$D/n425_C_B_10m.metadata.json
```

writes `results/p205-derivative-smith-loading/latest.json` (~85 s). Tests:
`python3 -m unittest tests.test_p205_derivative_smith_loading`.

The raw #205 production block was brought into the tree from
`origin/results/p205-norm5-conjugate-coalescence-20260829` (an unrelated-history delivery
branch); all 24 SHA-256 checksums verify from the repository root.
