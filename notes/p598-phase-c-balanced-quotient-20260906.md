# The balanced realization factors through the reflection quotient, and the odd-width break proves the comparison had teeth

**Date:** 2026-09-06
**Issues:** #600 (Phase C); closes the compute half of #598 and #596's Gate 2 prelude
**Evidence type:** numpy finite-horizon Gramians + Hankel SVD on the frozen #588 construction, A-vs-B per width. No sampling, no new width, no fit.

## What was asked, and the one hard constraint

#598 Phases A and B established that P398's coarsest positive lumping is the orbit partition of the reflection `R = (i -> (w-1)-i)` and that the `R`-odd sector of the single-point tilt carries no response for `R`-even sources and readouts. Phase C was deliberately left unrun because the previous session had no numpy/scipy and the two Gramians plus Hankel SVDs on a 1430-state operator are the wrong tool in pure Python.

This session had a real linear-algebra stack, so Phase C was run directly, two ways:

```text
A. the microscopic chain on Catalan(w) states,
B. the exact reduction to the R-orbit quotient on r_reflect(w) states,
```

on the frozen declared sources (all four `R`-even) and readouts, at `w = 4..8`, and compared on (1) finite-horizon Hankel singular values, (2) balanced order at the frozen tolerance, (3) response matrices, (4) transport under `H_even` at `eta = 0, ±1/4`.

## The protected dictionary factors through the quotient — bit-for-bit in the resolvable spectrum

Restricted to `D0` (the `R`-even dictionary) and the even readouts, A and B agree on every score the frozen construction produces:

| w | states → orbits | max |A-B| Hankel spectrum | response rel-Frobenius |
|---|---|---|---:|---|
| 4 | 14 → 10 | 1.5e-16 | 9.7e-16 |
| 5 | 42 → 26 | 3.3e-16 | 1.3e-15 |
| 6 | 132 → 76 | 3.0e-16 | 2.6e-15 |
| 7 | 429 → 232 | 3.4e-16 | 4.2e-15 |
| 8 | 1430 → 750 | 2.7e-15 | 1.1e-14 |

The Hankel spectrum is identical to the bit in every resolvable direction, and the numerical rank is preserved at every width. The balanced order at the frozen tolerance is the same on both chains, and the reduced generator (the quotient form of the tilted operator) agrees with the reduced microscopic generator to ~1e-10, which is the float64 arithmetic tail on a 1430-state object, not a Hankel direction.

## The odd sector is exactly inert

For `R`-even sources and `R`-even readouts the odd sector (`dim = Catalan - r_reflect`, i.e. 4/16/56/197/680) contributes no Hankel energy in either direction:

```text
w    odd dim   reach odd energy   observe odd energy
4       4         1.7e-17            5.9e-18
5      16         3.8e-17            4.5e-18
6      56         4.9e-17            1.4e-17
7     197         5.3e-17            7.5e-18
8     680         5.7e-17            7.6e-18
```

Exactly uncontrollable and unobservable, as the parity selection rule requires. This is why the balanced realization "sees" only `r_reflect` directions: the reduction is not an approximation.

## The exposed dictionary breaks at exactly the two odd widths — the control that makes the first claim a test

`halves_linked` is `R`-even only at even widths, so including it in the dictionary must break the quotient at `w = 5, 7` and only there:

| w | max odd content of exposed dict | restriction max |A-B| | symmetrization max |A-B| |
|---|---|---|---|---:|---|
| 4 | 0.000 | 1.7e-16 | — (identical) |
| 5 | **0.306** | **1.8e-2** | 1.7e-16 |
| 6 | 0.000 | 3.3e-16 | — (identical) |
| 7 | **0.250** | **8.7e-3** | 3.2e-16 |
| 8 | 0.000 | 4.3e-15 | — (identical) |

Two readout conventions were compared, kept apart because they mean different things:

- **restriction** — the microscopic readout is restricted to a representative of each orbit. At the two odd widths the odd part of `halves_linked` carries real signal, so restriction must disagree (1.8e-2 and 8.7e-3). It does.
- **symmetrization** — the readout is symmetrized (`v -> (v + Rv)/2`) before restriction, which removes exactly the odd part. Under this convention the two chains agree again to 1e-16.

The disagreement appearing under restriction and *disappearing* under symmetrization is the proof that the break is the `R`-odd part of `halves_linked` and nothing else — a run where B reproduced A on every readout at odd width would have been a bug, not a stronger result.

## Continuity gate: the numpy port reproduces the committed pure-Python numbers

The committed #588 artifact stores ten significant digits and takes its singular values from a Jacobi solver on the *squared* Gramian (which squares the condition number). The gate re-runs A on the frozen intervention and compares:

```text
w    spectrum above tolerance      max score difference
4    2.1e-11 over 8 entries        1.2e-10
5    3.3e-11 over 10 entries       4.7e-12
6    4.8e-11 over 11 entries       8.7e-12
7    3.7e-11 over 11 entries       5.0e-12
8    3.6e-11 over 11 entries       4.5e-12
```

All differences sit at the committed artifact's own storage precision. The port carries the same numbers; no B result below is an artifact of the port.

## Decision

```text
protected_factors_through_the_quotient : True
odd_sector_is_inert                    : True
exposed_disagrees_at_odd_widths        : True
exposed_agrees_at_even_widths          : True
```

All three legs of the pre-declared claim hold. The balanced realization factors through the reflection quotient for the invariant task, and the exposed dictionary's odd-width break is the control that shows the agreement was a theorem about the input/output object rather than an identity.

The claim is deliberately bounded: this is a statement about the frozen finite-horizon input/output description of one exactly known process, not about the state space, and not a percolation result. The quotient reproduces the task; it does not claim the `Catalan - r_reflect` "extra" states were unphysical.

## Second job — the projective design scan was not run (the referenced script does not exist)

The ticket's second job points at `scripts/p205_projective_channel_design.py` (#595 v2 / #596 Gate 2). That file does not exist anywhere in the repository — not on any branch, not in the history, not in any PR (checked `git log --all`, `gh search prs`, `gh search code`). What exists is `scripts/p205_channel_leverage.py`, which is the *demoted* #595 v1 (`|A4|/se(A4)`, explicitly "not the objective the downstream decision needs") and does not compute the Fieller confidence set on `(A4, delta)` that the ticket describes.

Implementing the projective design scan is a from-scratch build (joint `(A4, delta)` covariance, the #579 denominator-free Fieller/projective objective, the contamination-band sample multiplier, plus the three requested extensions), not a port. It is left as a separate piece of work and flagged here rather than silently substituted.
