# N=580 aspect ladder — result, 2026-09-05

**The ladder came back underpowered. That is a legitimate result and is
reported as such, not worked around.**

Artifact: `results/aspect-ladder-n580/latest.json`.
Design (frozen before any scoring run): `predictions/aspect_ladder_n580_20260905.yaml`.
Ticket: #567. Pilot (standard errors only, may not be pooled): `results/aspect-ladder-n580-pilot/latest.json`.

## What was run

200M samples per rung, 100 batches, seed `20260906`, replica offset
`900000000` (disjoint from the pilot's `800000000`). Three rungs of the frozen
ladder at 580 sites, one rung per machine on three ARM 16-vCPU boxes, each at
`--threads 16`. Thread count does not change the output — it was verified before
launch that the same seed/offset/samples/batches produce a bit-identical
histogram at one thread and sixteen — so the three rungs are one experiment
rather than three.

Per-rung `P4_S_prime` (the frozen spin-4 projector of the matching-odd slope):

| rung | modulus | `A4` | jackknife SE | relative |
|---|---|---:|---:|---:|
| r=1 | `i` | +0.0009016 | 0.0002491 | 28 % |
| r=2 | `2i` | +0.0029110 | 0.0002591 | 9 % |
| r=4 | `4i` | +0.0041318 | 0.0001951 | 5 % |

The standard errors extrapolate as `1/√n` from the pilot exactly as expected
(pilot r=1 SE 0.001024 at 10M → 0.000229 at 200M predicted; measured 0.000249).

## The score

| entry | value | SE |
|---|---:|---:|
| `A4(2i)/A4(i)` | 3.2285 | 0.9306 |
| `A4(4i)/A4(i)` | 4.5825 | 1.3167 |

The `r2_over_r1` entry keeps the spin-8 systematic (the leakage does not cancel
there), so it is reported and not used to exclude. Exclusion is decided on
`r4_over_r1` alone, which is the clean entry.

| law | r=4 prediction | z | outcome |
|---|---:|---:|---|
| weight-4 modular shape | 10.9908 | **−4.87** | excluded |
| bare aspect ratio | 4.0000 | +0.44 | **survives** |
| plain area scaling | 16.0 | −8.67 | excluded |
| no modulus dependence | 1.0000 | +2.72 | **survives** |
| weight12 `δ` | 0.00003 | +3.48 | excluded |
| weight8 `E8` | 120.80 | −88.3 | excluded |

## What it means

**The weight-4 modular shape is dead at r=4.** It predicted `10.99` and the
measurement is `4.58 ± 1.32`; it is excluded at 4.9σ. Plain area scaling (`16`)
dies harder, at 8.7σ. This is the same direction the N=290 fingerprint pointed
(`1.880 ± 0.177` against `11/4`), now at a different site count and in a
geometry where the leakage enters the discriminating ratio once rather than
twice.

**But the run is underpowered.** Two predictions survive at r=4: the bare aspect
ratio (`4.00`, z = +0.44) *and* no modulus dependence (`1.00`, z = +2.72). The
bare aspect ratio is the N=290 reading's post-hoc candidate, named in the frozen
file before any block ran, and it did **not** lose its one chance — but it also
did not uniquely win, because `1.00` sits only 2.72σ below the measurement, just
inside the 3σ exclusion line.

There is no optional stopping. The frozen rule is that a further run is a new
frozen design with a new file. The r=2 entry, at `3.23 ± 0.93`, is compatible
with `2.75` (z = +0.51) and `2.00` (z = +1.32) and cannot split them — which is
precisely why the ladder went to r=4 in the first place.

## Where this leaves the claims

The weight-4 modular shape survives as the paper's own prediction only if the
normalization that removes the same thermal-primary block holds — the additive
shape `A~(τ)` is not fixed by the Jordan relation (`docs/astra/Q2-additive-shape-ambiguity.md`).
This run tests the amplitude law, not the module identification, and it
excludes the amplitude law at 4.9σ at r=4.

The bare aspect ratio earned no unique prospective standing here — it is one of
two survivors, not the survivor. `docs/astra/Q2-additive-shape-ambiguity.md`
still says the bare aspect ratio has no prospective standing; whether a 4.9σ
death of the weight-4 shape plus a bare-aspect-ratio survival at z = +0.44
upgrades that sentence is a decision for the manuscript, not this note.
