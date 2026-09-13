# Erratum and second-round locks — span-spectrum diagnostic (2026-09-13)

**Scope.** Corrections to `span-spectrum-diagnostic-20260913.md` (commit `7226a2c`) and to
`scripts/span_spectrum_solve.py`, prompted by the owner's return read at the same commit
(PR #739, comment `5653178247`). No result is overwritten, no claim is upgraded; three claims
are downgraded to what was actually run, two source defects are fixed, one missing artifact is
restored.

## 1. Coverage claim corrected: the measured widths are 2–7, not 2–8

The note's §1/§2 say "widths 2-8" and "38 configurations"; the artifact that carries the data
(`results/geometric-consistency/span-spectrum-20260913.json`) holds **30 configurations,
widths 2–7**. Every w=8 statement in the note is about the *builder's capability and measured
build cost* (469 s, 2 505 625 frontier states, 604M transitions, 15.4 GB intermediate), not
about a solved spectrum: the w=8 solve never produced a result, and the `light` flag was
introduced for it (see §3). The width-8 row therefore must not be read as a measurement.
Tables §3/§4 are unaffected — they stop at w=7 by construction.

The results JSON itself was honest about this — its `honesty` list says "w=8 runs were still
streaming at delivery time and are NOT included". The overclaim is in the note's prose only.
Two further prose fields inside the JSON are wrong at claim level and are corrected in §6.1.

## 2. Tail-negligibility sentence corrected: 4 of 30 configurations are censored

§2 says "Tail fractions are negligible everywhere (<= 2e-6, mostly < 1e-20)". That is false.
With `D_MAX = 48` (w <= 6) the tail bin (`span > D_MAX`, excluded from all moments) carries:

| configuration | tail fraction | E[L] as printed |
|---|---|---|
| matching p=1/2, w=4 | **2.411e-02** | 14.253297 (censored) |
| matching p=1/2, w=3 | 2.113e-03 | 9.578554 (censored) |
| NN p=1/2, w=4 | 1.569e-06 | 7.087711 (censored) |
| matching p=1/2, w=2 | 1.294e-06 | 4.761844 (censored) |

**No published number is materially affected.** The two matching p=1/2 rows were computed but
never printed in §3/§4; the only published censored row is NN p=1/2, w=4, where a 1.6e-6 tail
shifts E[L] by ~1.2e-4 — below the two decimals quoted. All subcritical-p rows (p = 1/8, 1/4,
1/16) are complete to better than 1e-19. The solver now returns `moments.censored` and
`moments.tail_fraction` explicitly, and `tests/test_span_spectrum.py` locks the count at 4.

## 3. Source defect 1 — `light=True` skipped the `p_den**W` normalisation (fixed)

`light=True` handed the raw integer-weight float32 chain to `stationary()`, while the normal
branch divides by `P_W = p_den**W`. The stationary *distribution* is scale invariant; neither
solver is. Measured, on the committed tables, with this worktree's fixed copy sidelined:

| table | normal branch | pre-erratum `light` branch |
|---|---|---|
| `w2_nn.bin` (exact control) | exact rational | **RuntimeError: Factor is exactly singular** |
| `w4_nn.bin`, p=1/4 (n=11245, splu) | nu total 4.429301e-03, res 4.5e-15 | nu total 3.739144e-01 (**x84 wrong**), res **255** |
| `w5_nn.bin`, p=1/4 (n=52061, power) | nu total 1.344317e-03, res 2.8e-14 | **nan** (overflow to inf per sweep) |

With `(K^T - I)` nonsingular for the unscaled `K`, the `splu` path returned a non-stationary
vector with no error raised — a silent-wrong branch, not a crash. The fix divides exactly like
the normal branch; `w2_nn.bin` light and exact now agree to 0 difference (residual 1.4e-16),
locked by `test_light_chain_is_normalised`.

**No reported number used this path.** The 30 delivered runs carry modes: 6 `exact-rational`,
12 `float64(splu)`, 8 `float64(power)`, 4 `float64(power) + certificate`. Zero `light`. The
flag was dead code for every published value.

## 4. Source defect 2 — the certificate bounded a different vector than the one printed

The int64 residual certificate bounds `|observable(a/2^k) - observable(pi_exact)|` for the
**rational candidate** `a/2^k = round(pi * 2^k)/2^k`, while the printed observables come from the
float `pi`. The gap was implicit. The solver now emits both sides of it:

* `cert_candidate` — the candidate's spectrum, moments, and its `max_abs_diff_dh` and
  `rel_diff_sum_dh` against the float solve;
* `observable_bound_cert` — the certificate, unchanged in meaning;
* `observable_bound_float` — the same lemma applied to the float solve through its *reported*
  residual `float_residual_l1` (same `ginf/delta` factor);
* `ginf_num`, now recorded (it was not, which is why this could not be recomputed post hoc).

For the four certified runs, the candidate/float split is bounded by
`observable_bound_cert + observable_bound_float`, with `k >= 20`: the rounding term alone is
`<= 2^-20` relative, and the float residuals are 1.2e-13 … 1.7e-13, so the sum is far below the
six significant digits printed. The exact numbers are produced by the pending rerun in §6.

## 5. Mechanism reading withdrawn (owner's correction accepted)

§4 explained the declining ratio by "its span is the sum of ~w/mu piece heights". That
reasoning is wrong and is withdrawn. The span is
`max_i(y_i + upper_i) - min_i(y_i + lower_i) + 1` of the **whole** component; the engine's rule
was already that (per-component min-row offset, merge `max`, retire `span = offset`, verified by
`d_1 = p^w (1-p)^(2w)` and `sum_{h<=3} d_h = 9087/1048576`), but the one-line *explanation* was
not. If decorations are logarithmically Hausdorff-close to the skeleton — CIV §1.3.3, eq (1.10),
reprint p.11, then Theorem C's Brownian bridge for its bond/open-connection model — then
decorations do not force an extensive span and a range limit transfers to the full component.
What this delivery establishes is therefore the **measurement** only: no plateau at the
largest accessible width, ratio 2.7x the Brownian-bridge target and still decreasing. It does
not establish the additive mechanism, and it does not settle the asymptotics — the tagged
resolvent's conjecture (R) remains the right place for that.

Two related wordings are also narrowed: the closure `sum_h d_h + tail = nu_w` is a *theorem*
(the depth-clamped chain is an exact lumping) but is *verified* only where a certified `nu_w`
reference exists — w = 2, 3, 4 at the published p (plus #741's w=8 NN p=1/4 and matching p=1/8
rationals, unused here because no w=8 spectrum exists). §1's "at every width" means that.

Checked invariant, previously unstated: the 8-slot retire record cannot overflow. `nret` counts
distinct span bins among *retiring components*, each component contributes one bin, and at most
`W <= 8` components retire in a transition (`span_spectrum_build.cpp` line 182 caps W at 8), so
the `nret < 8` guard is unreachable and no retirement can be silently dropped.

## 6. Restored artifact and independent corroboration

`results/geometric-consistency/span-spectrum-20260913.json` (the 30-run table behind §3/§4) was
missing from `7226a2c` and is committed here, byte-identical to the artifact the tables were
computed from. Its headline moments regenerate from committed code via
`spectrum_moments(...)` (`test_published_moments_regenerate`), which they could not before: the
driver that produced §3/§4 was never committed.

The owner's one-lineage tagged resolvent — all heights, no bins, 3 963 reachable states at w=8
against this delivery's 2 505 625 — reports

| family | E L (this delivery, w=6 → 7) | E L (tagged, w=8) | CV^2 (w=6 → 7) | CV^2 (tagged, w=8) |
|---|---|---|---|---|
| NN p=1/4 | 3.93802 → 4.26752 | 4.56334396378038 | 0.15706 → 0.14000 | 0.127163500068044 |
| matching p=1/8 | 4.06011 → 4.40122 | 4.711805413158092 | 0.13990 → 0.12660 | 0.116004934830910 |

Both columns continue the same trend within the quoted digits — two independent chains agreeing
on the measured object, and a route that removes the width-8 age-state computation this
delivery could not finish.

## 6.1 Two claim-level corrections inside the results JSON

The JSON is committed with every number byte-identical (the `runs` array is unchanged, verified
by a sorted re-serialisation equality check before writing). Two prose fields were false and are
replaced, so that the machine-readable artifact cannot be read as claiming a validation it never
ran:

* `honesty[2]`, was: *"the exact int64 certificate is reported where the scale fits; w=7,8 use
  float64/float32 power iteration validated by exact closure against certified nu_w"* — wrong
  twice: no run ever used float32, and no closure against a certified `nu_w` was run at w >= 5
  (`gap_float` is `None` for all 11 of those runs). Now: the certificate count is 4 of 30, the
  other 26 are float64 solves with the float residual reported, and closure is verified at
  w = 2, 3, 4 only.
* `findings.interpretation` carried the additive-span mechanism. Now marked WITHDRAWN with the
  corrected reading of §5.

## 7. Pending, and not claimed here

1. Width-8 spectrum (needs the fleet; the `light` path is fixed but the state count is 2.5M).
2. Rerun of the four certified configurations to print the certified candidate of §4 exactly.
3. Closure at w >= 5 against a certified `nu_w` at those widths (or the tagged resolvent's
   total), the one validation this delivery could not run.
