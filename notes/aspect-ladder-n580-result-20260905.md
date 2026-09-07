# N=580 aspect ladder — result, 2026-09-05

**The ladder came back underpowered. That is a legitimate result and is
reported as such, not worked around.**

Artifact: `results/aspect-ladder-n580/latest.json`.
Design (frozen before any scoring run): `predictions/aspect_ladder_n580_20260905.yaml`.
Ticket: #567 (run), #575 (deterministic replay that measured the covariance).
Pilot (standard errors only, may not be pooled): `results/aspect-ladder-n580-pilot/latest.json`.

## What was run

200M samples per rung, 100 batches, seed `20260906`, replica offset
`900000000` (disjoint from the pilot's `800000000`). Three rungs of the frozen
ladder at 580 sites, one rung per machine on three ARM 16-vCPU boxes, each at
`--threads 16`. Thread count does not change the output — it was verified before
launch that the same seed/offset/samples/batches produce a bit-identical
histogram at one thread and sixteen — so the three rungs are one experiment
rather than three.

**#575 re-ran this to recover the delete-one replicates.** The first pass threw
away the per-batch delete-one values (a GOVERNANCE §2 B violation), so the
cross-rung covariance could only be reconstructed. The replay ran the exact same
parameters and every rung reproduced its committed `histogram_sha256`
bit-for-bit (`r1` `f71c0faf…`, `r2` `68b6253a…`, `r4` `5d00fe00…`), which is the
independent confirmation that the engine is deterministic in the way the
"threads do not change the output" claim assumes. The per-rung shards, carrying
the delete-one replicates, are committed at `results/aspect-ladder-n580/shards/`.

Per-rung `P4_S_prime` (the frozen spin-4 projector of the matching-odd slope):

| rung | modulus | `A4` | jackknife SE | relative |
|---|---|---:|---:|---:|
| r=1 | `i` | +0.0009016 | 0.0002491 | 28 % |
| r=2 | `2i` | +0.0029110 | 0.0002591 | 9 % |
| r=4 | `4i` | +0.0041318 | 0.0001951 | 5 % |

The standard errors extrapolate as `1/√n` from the pilot exactly as expected
(pilot r=1 SE 0.001024 at 10M → 0.000229 at 200M predicted; measured 0.000249).

## The statistic, and the deviation it recorded

The ratio `A4(4i)/A4(i)` is `4.5825`, but its denominator `A4(i)` is only
3.6σ from zero. Dividing first and forming `(Rhat − R0)/SE(Rhat)` assumes the
ratio is roughly normal, which a weak denominator breaks — so the committed
runner's ratio z-test was replaced by the **Fieller contrast** on
`Y − R0·X`, where `X` and `Y` are batch means and the CLT applies to them
directly. That is a deviation (the ratio z was effectively frozen by being in
the committed runner), and it is recorded in
`predictions/aspect_ladder_n580_20260905.yaml` rather than applied quietly. The
pre-registered numbers stay in the artifact under `ratio_z_not_used`.

## The covariance, now measured

The pair `(A4(i), A4(4i))` is measured on the same batches, so its covariance is
a delete-one jackknife quantity. The first pass dropped the replicates and
backed ρ out of the ratio standard error, giving **ρ = −0.1526**. The replay
measured it directly:

| quantity | reconstructed | measured (#575) |
|---|---:|---:|
| `cov(A4(i), A4(4i))` | −7.45e−9 | −8.01e−9 |
| **ρ** | **−0.1526** | **−0.1648** |

The measured ρ is close to the reconstruction, and rescoring with it moves no
verdict at 3σ — the conclusions do not rest on the correlation (there is a test
pinning that over ρ ∈ [−0.5, +0.5]).

## The score

`r2_over_r1` keeps the spin-8 systematic (the leakage does not cancel there), so
it is reported and not used to exclude. Exclusion is decided on `r4_over_r1`
alone. The 3σ Fieller interval on `A4(4i)/A4(i)` is **`[2.40, 27.47]`** — bounded
only just, because the denominator clears 3σ by a hair.

| law | r=4 prediction | Fieller z | ratio z (not used) | outcome |
|---|---:|---:|---:|---|
| weight-4 modular shape | 10.9908 | −2.08 | −4.87 | **survives** |
| bare aspect ratio | 4.0000 | +0.50 | +0.44 | **survives** |
| plain area scaling | 16.0 | −2.56 | −8.67 | **survives** |
| no modulus dependence | 1.0000 | **+9.48** | +2.72 | excluded |
| weight12 `δ` | 0.00003 | +21.2 | +3.48 | excluded |
| weight8 `E8` | 120.80 | −3.48 | −88.3 | excluded |

The statistic change moves verdicts in **both** directions, which is why it is a
fix rather than a preference: `no_modulus_dependence` goes from surviving at
+2.7 (ratio) to excluded at +9.5 (Fieller), while the weight-4 shape and plain
area scaling go from excluded at −4.9 / −8.7 to surviving at −2.1 / −2.6.

## What it means

**Three predictions survive at r=4** under the Fieller contrast — the bare
aspect ratio (`4.00`, z = +0.50), the weight-4 modular shape (`10.99`,
z = −2.08) and plain area scaling (`16`, z = −2.56) — so the run is
**underpowered**: it does not identify a single amplitude law. `no_modulus_dependence`
(`1.00`) is the one clean exclusion at +9.5σ.

The weight-4 modular shape is no longer excluded here. Under the ratio z-test it
was dead at −4.9σ; the Fieller contrast, which does not punish the weak
denominator the way the ratio test does, brings it back inside 3σ at −2.1σ. The
r=2 entry, at `3.23 ± 0.93`, is compatible with `2.75` (z = +0.51) and `2.00`
(z = +1.32) and cannot split them.

There is no optional stopping. The frozen rule is that a further run is a new
frozen design with a new file. A law that fits three points is a law that fits
three points; the run says the N=580 ladder, as sized, cannot adjudicate between
the three survivors.

## Where this leaves the claims

The weight-4 modular shape survives as the paper's own prediction only if the
normalization that removes the same thermal-primary block holds — the additive
shape `A~(τ)` is not fixed by the Jordan relation
(`docs/astra/Q2-additive-shape-ambiguity.md`). This run tests the amplitude law,
not the module identification, and under the correct statistic it neither
excludes nor confirms the weight-4 amplitude law.

The bare aspect ratio earned no unique prospective standing here — it is one of
three survivors, not the survivor. `docs/astra/Q2-additive-shape-ambiguity.md`
still says the bare aspect ratio has no prospective standing; whether this run
upgrades that sentence is a decision for the manuscript, not this note.
