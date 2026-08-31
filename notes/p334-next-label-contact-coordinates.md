# P334: most safe R0 insertions close loops; observed direct double births use one carrier

The already sampled next-label stream now has exact checkpoint contact
coordinates. Across the four N/orientation rows, **63.44–63.66% of safe R0
insertions create contractible graph cycles**. All197 sampled direct0->2
draws arise from one old occupied component supplying two independent winding
directions through at least three contacts. The distinct two-component2+2
architecture is absent from this finite sample.

## New finite geometry

| Size / orientation | R0 sampled label draws | Safe R0 | Safe R0 with new contractible cycles | Fraction among safe | Direct0->2: one component /2+2 |
| --- | ---: | ---: | ---: | ---: | ---: |
| N325 first | 85216 | 76690 | 48824 | 63.6641% | 63 /0 |
| N325 second | 86128 | 77436 | 49222 | 63.5648% | 61 /0 |
| N425 first | 89552 | 82090 | 52167 | 63.5485% | 39 /0 |
| N425 second | 85872 | 78805 | 49992 | 63.4376% | 34 /0 |

The197 direct-birth draws represent187 distinct
`(N,counter,next_label,orientation)` events. Their contact structures split as

| Contact structure | Number of direct0->2 draws | New graph cycles |
| --- | ---: | ---: |
| e=3,c=1: all three contacts in one old component | 140 | 2 |
| e=4,c=1: all four contacts in one old component | 45 | 3 |
| e=4,c=2: three contacts in the active component, one in another | 12 | 2 |
| Two distinct components with two contacts each and independent lines | 0 | 2 |

Thus the observed direct births are not merely high-contact sites: within the
component, relative lift addresses span both period directions. The45
four-contact/one-component events also create one independent nullhomologous
cycle in addition to the two-dimensional ambient image.

Safe loop closure and merging old components are not mutually exclusive.
For example, N325 first has15055 safe draws with neither,12811 with merging
only,35682 with loop closure only, and13142 with both. The corresponding full
four-cell tables are saved for all orientations. A single categorical
"loop versus merger" label would discard the overlap; the raw e/c fields
retain both coordinates separately.

## What is reconstructed

Source forks are `e32a85939279b8574278024d647b56d2d1485247`, generated with
producer/backend `a3249a59` on original prefix archive `9c495ab1`. Only the
original counter permutation and its first k0 insertions are reconstructed.
The next labels are read from the existing fork rows; no next-label RNG,
suffix sampling/replay, full birth-clock solve or DP is invoked.

For each sampled next label v, count occupied incident edge contacts e and
touched old components c. Then `e-c` is the number of newly independent graph
cycles. For R0, the theorem `e67d9b90` makes their ambient image exactly the
span of the integer windings

```
P^(-1)[alpha_e-alpha_anchor(C)],   alpha_e=step_e-p_C(u_e),
```

where every difference stays within one old component. Absolute potentials
from different components are never subtracted. The kernel rank is
`e-c-contact_rank`; for a safe R0 insertion this equals all e-c cycles.
For direct rank2, architecture1 means one component's own generators already
have rank2; architecture2 means two rank1 components provide independent
directions, which with four NN ports requires2+2 contacts.

For R1/R2, e,c,e-c and the recorded old/after ranks are retained. Their
`contact_rank`, `r0_null_cycle_rank`, and `r0_rank2_arch` are-1. Parallel
essential cycles can be created while an R1 rank stays unchanged; these are
not silently called contractible. Likewise no spanning-tree-dependent R1
contact image is named as a new physical observable.

## Artifact and inference unit

Raw lock: `959a7fa26677c416b874d272f1ba66523fb38f73`.
Code: `b044e6452d3342f496d77a67f13319caa06c92fe`.
The forty compact gzip files occupy about6.7MB and contain640,000 wide rows
(two orientations each) keyed by
`N,batch,counter,quartet,group,next_label`. Each orientation supplies
`oldrank,rank_after,e,c,new_cycles,contact_rank,r0_null_cycle_rank,r0_rank2_arch`.
Repeated labels remain repeated draws; duplicate suffix rows are suppressed.
Caching reduces the actual topology queries to608,964 unique prefix/label
pairs across the two sizes. Reconstruction took0.979s and1.094s locally.

The data and original batch IDs have been handed directly to the shared
response/covariance analysis. Counts here are descriptive draws from the
same original prefix/fork dependency block, not independent binomial trials,
a statement that2+2 is impossible, or a continuum arm exponent. No new
response model or covariance scan is performed in this contact note.
