# N=580 aspect ladder — pilot, 2026-09-05

**Pilot only. The standard errors are the output. The central values are not, and
must not be read as a measurement or pooled into the scoring run.**

Artifact: `results/aspect-ladder-n580-pilot/latest.json`.
Design: `predictions/aspect_ladder_n580_20260905.yaml`.
Ticket for the production run: #567.

## What was run

10M samples per rung, 20 batches, seed `20260905`, replica offset `800000000`,
three rungs of the frozen ladder at 580 sites.

| rung | `A4` | jackknife SE | relative | per-difference noise |
|---|---:|---:|---:|---:|
| r=1 | +0.003119 | 0.001024 | 33 % | 0.001963 |
| r=2 | +0.000090 | 0.001277 | 1416 % | 0.002448 |
| r=4 | +0.004158 | 0.000929 | 22 % | 0.001781 |

Every rung shares `Δcos4 = 8064/4205`, so the three noise figures are directly
comparable and none of them is the bottleneck by more than 40 %.

## Why the central values may not be read

A 200,000-sample pilot of this same channel once returned an amplitude **five
times** its 20M value, and a sample count projected from it was wrong by about
forty times. Small-sample estimates of this projector run high.

Standard errors are a different matter: they extrapolate as `1/√n` and were
checked to. A 1M pilot at N=290 predicted the committed 200M standard error to
within 22 %, which is inside the spread you expect from a 10- to 20-batch
jackknife. That is the whole reason a pilot is worth running at all here.

## Budget

Throughput is about **33 s per million samples** at 580 sites, measured.

| target relative precision on `A4(i)` | samples/rung | hours/rung |
|---|---:|---:|
| 30 % | 86 M | 0.8 |
| 20 % | 193 M | 1.8 |
| 15 % | 344 M | 3.2 |

**200 M per rung, three rungs, about 5.5 hours.** That puts ~20 % on the
denominator of both score entries, which separates `4.00` from `10.99` at r=4
with room to spare.

## What this pilot replaced

The design was at N=1300 with three orientations per rung until a 1M pilot there
returned two results, neither of them the one it was run for.

**It found a bug.** The analysis path returned zeros at N=1300 — the binomial
tail's recurrence starts at `(1−p)^N`, which underflows to exactly zero near
**790 sites** at the percolation threshold and then stays zero, silently. That
bound had never been noticed and had capped every large-`N` plan in the
repository. Fixed by anchoring the recurrence at the mode; `N ≈ 4000` now
analyzes cleanly. Follow-up audit: #570.

**It priced the design out.** Per-difference noise at N=1300 was `0.0131` per 1M
against `0.0065` at N=290 — a factor 2.0 for a 4.5× larger torus — while the
amplitude falls roughly as `N^-5/4`. A decisive ratio there is about three orders
of magnitude beyond what we can spend.

N=580 needs no three-orientation fit, which is what N=1300 was bought for: its
r=1 and r=4 rungs carry the same spin-8 leakage, so the bias cancels to leading
order in the ratio that discriminates. Three runs instead of six, at 580 sites
instead of 1300, and a systematic handled by geometry rather than by fitting.

## The honest summary

Two pilots, about twenty minutes of compute between them, turned a design we
could not afford and could not have analyzed into one that runs in an afternoon.
Neither pilot measured anything about the lattice, and that was not the point.
