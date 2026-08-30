# P250/P333: the adaptive order signal survives medium Gaussian tori

Date: 2026-08-30. Status: **frozen two-geometry intervention pilot completed**.

The exact `L=3/4` adaptive cut/connector witness was transported without changing its rule to two medium Gaussian quotients.  This is the first actual-data answer to the narrow question “does the constructed order response leave the tiny exact witnesses?”  It does.

## Frozen design

Commit `ab28bac` froze the runner and preregistration before the Huawei runs:

| geometry | Gaussian period | matrix | samples / batches | counter interval |
|---|---|---|---:|---|
| N325 | `(17,6)` | `[[17,-6],[6,17]]` | 20,000 / 40 | `[10250000000,10250020000)` |
| N425 | `(16,13)` | `[[16,-13],[13,16]]` | 20,000 / 40 | `[10260000000,10260020000)` |

Each replica uses a counter-selected translation and cycles through four rotations and their four reflected partners.  The marked stencil is the exact-witness pattern `a_D=(0,0), a_J=(1,1), c=(0,1)`.  Undefined partial rectangles score zero in the primary unconditional response; their frequency is reported separately.

Every defined primary row is paired with the complement field, exchanged NN/matching hand, and swapped anchors.  Freezing the initially selected `D0/J0` sites supplies an exact order-null control.

## Results

| geometry | defined | definition probability | unconditional `R_minus` | `R_minus` given defined | unconditional `R_plus` |
|---|---:|---:|---:|---:|---:|
| N325 | 764 | `0.03820 +/- 0.00118` | `0.06785 +/- 0.00235` | `1.7706 +/- 0.0144` | `-0.00290 +/- 0.00088` |
| N425 | 739 | `0.03695 +/- 0.00125` | `0.06655 +/- 0.00226` | `1.8032 +/- 0.0141` | `-0.00190 +/- 0.00078` |

All 1,503 defined medium-torus rectangles have nonzero positive order response:

```text
N325: Rminus=1:171, Rminus=2:593
N425: Rminus=1:147, Rminus=2:592
```

The primary responses agree across geometries:

```text
N325 - N425 = 0.00130 +/- 0.00326  (z=0.40).
```

Definition probability also agrees (`z=0.73`).  Conditional amplitude differs by only `-1.62` combined standard errors.

## Controls

Across both runs:

- typed defined/undefined mismatches: `0`;
- typed support-site mismatches: `0`;
- maximum `Rminus-primary - Rminus-dual`: `0`;
- maximum `Rplus-primary + Rplus-dual`: `0`;
- frozen-support order-null failures: `0`.

The rotation `q1`, rotation `q2`, and reflection-odd leakage modes are all small.  The largest absolute leakage score is `0.84` standard errors at N325 and `1.07` at N425.  The analogous availability projections are also unresolved.  Thus the observed scalar response is not being carried by one chosen stencil orientation or reflection hand at this precision.

`Rplus` is much smaller than `Rminus` but mildly negative in both geometries.  It is secondary: the typed pairing makes primary/dual `Rplus` exactly opposite, and this pilot was not designed to identify its mechanism.

## Scientific reading

The useful conclusion is specific: the physical, state-dependent intervention introduced by `6fbbe5e` does not collapse when the quotient grows from the exact `L=3/4` witnesses to N325/N425.  Its unconditional magnitude and definition rate are stable across two independent Gaussian geometries, while every exact null and the orientation-leakage controls behave correctly.

This is not evidence that an unperturbed percolation trajectory possesses hidden memory.  Noncommutation is introduced by the frozen adaptive support rule itself.  Nor is `Rminus` identified with a CFT operator, a Jordan cell, a universal amplitude, or an asymptotic exponent.  The result is a finite intervention response and a working production observable, nothing more.

## Provenance

Environment: Huawei `DevEnvC_HZsCM6`, ID `033945d8bf8b47a7acf475c595169e07`, 16 Kunpeng cores, Python 3.9.9.

```text
N325 wall: 5.732 s
N425 wall: 7.461 s
N325 batches sha256: f5127be480375d9226e2b1f54b283d99d9af6a2ecbc91ab7ae2a168a42393d05
N425 batches sha256: 5f70875add2af083423ae7142e7941a6ffe9c20735f7a5b72a856a503d1624c1
```

Artifacts live under `results/server-20260830/P250-adaptive-order-pilot/`.  No production extension was launched.
