# Gaussian commuting-square root character: production power gate

Status: source-only planning analysis for Issue #158.  No N=650/N=850 target
data were generated or read.

## 1. The fixed-`p` and root characters are different objects

For norm `Q=10`, the exact angular factor is

```text
DeltaCos4_child / DeltaCos4_parent = 14/25.
```

If the fixed-probability matching contrast has the thermal H4 scaling
`DeltaM ~ DeltaCos4 N^(-13/8)`, while the mean thermal slope scales as
`Mbar' ~ N^(3/8)`, then

```text
r_M     = (14/25) 10^(-13/8) = 0.0132796927517053,
r_slope = 10^(3/8)             = 2.37137370566166,
r_root  = r_M/r_slope          = 7/1250 = 0.0056.
```

The root factor must not be applied to fixed-`p` `DeltaM`.  Conversely, the
fixed-`p` factor must be divided by the thermal-slope factor before it becomes
a root prediction.  The tiny mismatch between direct transported roots and
`-DeltaM/Mbar'` below is exactly the already measured P45 nonlinear closure
error (`7.8e-5` and `8.8e-5` relative), not a character discrepancy.

## 2. Source-implied leading targets

Using the revealed P45 roots, fixed-`p` gaps and slopes gives:

| lineage | target `DeltaM(p_ref)` | target `Mbar'` | target root gap |
|---|---:|---:|---:|
| 65 -> 650 | `+1.5091265e-5` | `19.87027` | `-7.5954907e-7` |
| 85 -> 850 | `+1.0715235e-5` | `21.95120` | `-4.8809610e-7` |

These are leading-character planning targets, not target observations and not
a fitted finite-size correction.

## 3. The direct production is severely underpowered

P57 supplies a much better large-`N` variance proxy than P45.  Its 500M
same-orientation-pair runs give the per-replica fixed-`p` variance coefficients

```text
N325 lineage: n Var(DeltaM) = 0.5143361
N425 lineage: n Var(DeltaM) = 0.5039902.
```

Transporting these coefficients to N650/N850 and retaining the archived P45
parent uncertainty, a two-sided `alpha=.01`, 80% known-direction test against
a zero child contrast needs:

| target | target samples | P57-scaled 8-thread time | expected z at 500M |
|---|---:|---:|---:|
| N650 | `27.21B` | `96.5 h` | `0.470` |
| N850 | `55.18B` | `251.6 h` | `0.337` |

Root-gap propagation gives `27.208B` and `55.189B`, respectively.  The near
identity is the slope conversion above: scoring roots instead of fixed-`p`
`DeltaM` does not create information.

A joint matched-direction score can reach the same detection criterion with
`17.88B` samples in *each* target stream (`35.76B` total), but that only detects
the two-lineage predicted direction.  It does not demonstrate separate
commuting-square closure in either lineage.

More importantly, detection is not precise measurement of the rational
character.  With the current parents, one-standard-error relative precision
requires:

| target | 20% relative SE | 10% relative SE | 5% relative SE |
|---|---:|---:|---:|
| N650 | `60.44B` | `306.60B` | impossible from archived-parent floor |
| N850 | `129.38B` | `1.118T` | impossible from archived-parent floor |

Thus a 500M or 1B run would be a runtime/variance pilot only.  It could not
fairly be presented as a test of `7/1250`.

## 4. A norm-10 covering CRN cannot rescue this design

The existing integer-period engine already uses the useful coupling: the two
same-`N` orientations share one HNF-label permutation.  A new parent-child
cover coupling attacks only the transported parent term, whose coefficient is
`0.01328` in fixed-`p` space (`0.0056` in root space).

Cauchy--Schwarz gives a coupling-independent no-go.  Even with a fresh parent
run at equal depth and unattainable absolute correlation one, the maximum
total variance reductions are only:

```text
65 -> 650: 2.15%
85 -> 850: 2.34%.
```

No implementation choice can exceed those bounds.  A symmetric ten-fiber
measure-preserving parent-rank map could be built as a coupling pilot, but it
cannot materially change the production decision and is therefore lower
value than target samples.

Also, the direct, norm-2-then-norm-5 and norm-5-then-norm-2 constructions are
not three stochastic observables.  Exact Gaussian multiplication produces the
same target period lattice up to the already audited canonical relabeling.
Re-running the relabelings would duplicate one target law, not add a second
commuting-square measurement.

## 5. Backend readiness and the minimal remaining software gap

`src/threshold_rank_integer_period_mc.cpp` already accepts the targets through
its arbitrary-matrix CLI:

```text
N650: ((23,-11),(11,23)) / ((19,-17),(17,19))
N850: ((29,-3),(3,29))   / ((27,-11),(11,27))
```

A local 100-replica smoke run accepted both, recorded Smith `(1,650)` and
`(1,850)`, and passed the engine self-test.  No new production runner is
needed.

The scorer is not yet production-ready for these sizes:
`score_angular_root_amplitude.py` hard-codes `[65,85]` and a common provenance
block.  Before any target launch, parameterize the size pair and allow
independent target counter groups while retaining full per-batch root,
`DeltaM`, slope and covariance reconstruction.

## Decision

Do not launch the nominal N650/N850 production at 500M-scale.  The useful next
gate is either:

1. accept a `50M x 50 batches` target-only pilot as a variance/runtime audit,
   with no scientific character claim; or
2. redesign the observable so that the child signal is not suppressed by
   `7/1250`--for example, score an amplitude-normalized H4 coordinate whose
   deterministic `N^2/DeltaCos4` rescaling is applied inside every batch.

The second option changes numerical conditioning, not Fisher information, if
it is only a deterministic rescaling.  Genuine improvement requires an
additional strongly correlated target observable/control variate, not a new
parent-child cover map.

