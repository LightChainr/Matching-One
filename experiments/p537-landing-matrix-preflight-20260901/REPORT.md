# P537 provisional clean-two-bridge landing-minor audit

## Decision

Under the explicit **provisional** landing contract below, the finite N25
rank-one route fails.  For each of the two nonzero C4 landing orbits and for
their sum, all six `2 x 2` minors of the four-readout source/thermal matrix are
strictly nonzero at the saved pooled root.  The minor remains nonzero after the
original P4 row projection and after the fixed-M Schur column elimination.
The complete root-conditioned mixed Hessian `T_p/M_p` is also strictly
nonzero.

This is an exact finite counterexample to the **provisional contract**.  It is
not yet licensed as a counterexample to the repository phrase “ordinary
four-arm/no-extra-branch”, because that phrase has no formal event definition
or transition-row contract in the checked-in files.

## Why the checked-in N13 aggregate is not enough

`thermal_gate.py` can regenerate every N13 state/site/pair record, and it was
used here to verify the Bell-key and rank semantics.  Under the same
provisional six-block source at `p=1/2`, the single-geometry `(q,E)`
source/thermal minor is already

```text
10775865/34359738368 = 0.00031361894798465073 != 0.
```

That is only a scalar finite warning.  N13 has no inequivalent same-norm
Gaussian representation, so it cannot supply the original axis/tilted P4
row.  The checked-in audit also exports only aggregate flip counters and one
witness, rather than all landing rows.  For the actual pooled-root/P4/Schur
test below, this audit therefore reused the existing exact N25 joint-source
producer and changed only its Bell-key grouping.

There is a second exact N13 warning that does not use the provisional landing
filter.  Grouping every physical record only by rank transition at `p=3/5`,
the `(kernel reconnection, readout pivotal) x (0->1,1->2)` matrix for `E`
has determinant

```text
-2193536256714011846008984563072
 /600408611717284657061100006103515625
= -3.6534057205476687e-6 != 0.
```

The `0->2` readout entry is exactly zero while its kernel entry is nonzero.
This kills the stronger claim that rank-transition aggregation is
automatically rank one, but it still mixes every landing type and precedes
C4/root-Schur projection.  It is therefore a guardrail rather than the #537
falsifier.

## Provisional landing contract

An eight-port Bell partition is retained exactly when:

1. it has six exterior blocks;
2. exactly two blocks occur once at each marked four-port group;
3. the other four blocks are singleton ports; and
4. the two shared landing positions are adjacent at each mark.

This is the smallest port-level meaning of “two bridges with no port branch”.
Separate C4 rotations at both marks plus mark exchange give four orbits among
the 72 six-block two-bridge partitions.  Two orbits have `g16=0`; the two
retained orbits each have 16 partitions and `g16=4`.  They differ by whether
the cyclic order of the two shared labels agrees (`clean_same`) or reverses
(`clean_reversed`) between the marks.  `clean_total` is their sum.

The C++ enumerator aborts if a retained physical record has `g16 != 4`.
`validate_landing_contract.py` independently reconstructs all 4,140 Bell-8
partitions and the four C4xC4/site-exchange orbits.

## Exact matrix and Schur invariance

At the saved N25 pooled root

```text
p in [0.5926655393282267461514256998465469006335,
      0.5926655393282267461514256998465469006339]
```

the four matrix rows are

```text
(q_axis, E_axis, q_tilted, E_tilted),
```

and the columns are the Bernoulli thermal tangent and the centered clean
source tangent.  The two geometries are `(5,0)` and `(4,3)`.  Every entry is
evaluated by exact integer K profiles and Fraction interval arithmetic.

The original row projection is

```text
M = (q_axis+q_tilted)/2,
Y = (E_axis-E_tilted)/(1152/625).
```

The projected determinant is

```text
det [[M_p,jM],[Y_p,jY]] = M_p*jY-Y_p*jM.
```

Holding M fixed replaces the source column by
`source-(jM/M_p)*thermal`.  This is a column addition, so every `2 x 2`
minor is unchanged exactly.  The residual source response in the Y row is
the projected determinant divided by `M_p`.

| retained block | nonzero raw row minors | P4/root minor | fixed-M source residual | full `T_p/M_p` |
|---|---:|---:|---:|---:|
| clean_same | 6/6 | `+2.7355278002402004e-7` | `+4.643561586827378e-8` | `+7.600979820643663e-7` |
| clean_reversed | 6/6 | `-2.5116785501799413e-5` | `-4.263577227418322e-6` | `+5.0462349846123006e-6` |
| clean_total | 6/6 | `-2.4843232721775393e-5` | `-4.217141611550048e-6` | `+5.806332966676667e-6` |

All displayed signs are certified by outward rational intervals.  The widest
interval in the three projected minors is below `5.7e-40`; the widest interval
in the three mixed Hessians is below `8.8e-40`.

The six unprojected minor midpoints for `clean_total`, in row-pair order, are:

| rows | determinant |
|---|---:|
| axis q, axis E | `+6.208293160316876e-4` |
| axis q, tilted q | `+1.2489552382152975e-4` |
| axis q, tilted E | `+6.654697886038648e-4` |
| axis E, tilted q | `-6.214276684643301e-4` |
| axis E, tilted E | `+1.1224656709745817e-5` |
| tilted q, tilted E | `+6.683692889977057e-4` |

## Enumeration and controls

The modified producer performs six exact `2^24` fixed-origin traversals:
three landing modes on each geometry.  Each traversal takes about 3.5 seconds
on the local Mac; there is no Monte Carlo, OpenMP, server job, or sampling
error.  The two orbit profiles add exactly to `clean_total` in every K slice
for `S2`, `q*S2`, and `E*S2`.

The unmodified `all` mode was also rerun.  Its axis and tilted CSV files match
the saved producer files byte for byte, with SHA256 respectively

```text
d1369a8b6fdb782be78cbe52ae2923e41f993325717a6f113eea02a0b19bdcce
ea464ad416054aefba92a45ba459e64b506257a0c54ea49d17b3a462696c553f
```

The independent scorer reproduces the published full-source value
`J2/A_N=-5.905706006949678e-5` exactly at the displayed precision.

## What is missing for the canonical “ordinary four-arm” claim

The checked-in next-target note names the matrix test but does not define its
finite state space.  Before promoting this provisional counterexample, the
repository must freeze these semantics:

1. **Allowed Bell-8 source landings.** State whether private ports must be
   singleton, whether a shared component may touch more than one port at a
   mark, and whether adjacent marks sharing a physical edge are included.
2. **Meaning of “no extra branch”.** A Bell partition records terminal
   incidence only.  If the exclusion concerns branching away from the eight
   ports, the current exact files cannot decide it.  The producer must retain
   per shared component its off-port degree/branch flag (and, if relevant,
   cycle/homology flag).
3. **Thermal landing event.** The present matrix uses the complete Bernoulli
   K derivative, which is exactly the thermal tangent but does not label the
   flipped site z.  A transition-resolved contract needs, for every selected
   `(configuration,x,y,z)`, `K`, the source Bell key before and after z, the
   two `g16` values, z's cyclic four-port partition, rank before/after, q/E
   before/after, geometry, and the C4 character weight.
4. **Transition row basis.** Freeze whether rows are `(0->1,0->2,1->2)`, the
   four original q/E geometry readouts used here, or another basis.  Without
   that declaration “all minors” is not a unique calculation.
5. **Angular projection.** Freeze whether C4 means an invariant sum over the
   two nonzero local orbits, a signed H4 character, or only the existing P4
   axis/tilted geometry projection.  This audit uses the invariant local sum
   and the original P4 geometry row.

If the intended event is exactly the six-block port-level contract above,
the finite rank-one route is already dead and should stop.  If the intended
event adds a global branch-free condition or a z-local landing restriction,
the extra fields above are the smallest code extension needed; the aggregate
N13 flip counters and current N25 K profiles cannot reconstruct them after the
fact.

## Reproduction

The self-contained artifact lives in
`experiments/p537-landing-matrix-preflight-20260901`.  Its verifier rebuilds
eight exact N25 traversals in a temporary directory, checks every committed
profile byte for byte, recomputes the Bell orbit table and the compact N13
transition result, and leaves no binary in the repository.

```bash
python3 experiments/p537-landing-matrix-preflight-20260901/verify.py
```

Machine-readable exact outward enclosures and all six minors per mode are in
`result.json`; Bell orbit validation is in `landing_contract_validation.json`.
