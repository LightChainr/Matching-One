# Frozen macro-window joint-U collector prototype

`scripts/regular_pair_macro_joint_u_sampler.cpp` implements the raw collector
specified by [`analysis/regular_pair_macro_joint_u_contract.json`](../analysis/regular_pair_macro_joint_u_contract.json).
It is a compile-ready prototype; this delivery did not compile it, run the
pilot, regenerate the kernel, allocate a server, or score a response.

## Fixed invocation and pairing

One invocation handles the axis and tilted geometries at one size so their
occupation bit vector and anchor-index list cannot silently diverge:

```bash
c++ -O3 -std=c++17 scripts/regular_pair_macro_joint_u_sampler.cpp \
  -o /new/path/regular_pair_macro_joint_u_sampler

regular_pair_macro_joint_u_sampler \
  100 2026090101001001 2026090101001002 100 500 16 \
  analysis/regular_pair_spatial_kernel.tsv /new/path/N100.csv

regular_pair_macro_joint_u_sampler \
  400 2026090101004001 2026090101004002 100 500 16 \
  analysis/regular_pair_spatial_kernel.tsv /new/path/N400.csv
```

The executable rejects another size, seed, batch count, configurations per
batch, or anchor count. Occupations use `mt19937_64`, one word per canonical
site index, with threshold `10934234699625173385/2^64`. Anchors use a separate
`mt19937_64` stream and exact uint64 rejection sampling. Partial Fisher-Yates
selects16 distinct indices without replacement. At each N, the same bits and
indices are passed to the two geometries; their E-then-N quotient-BFS index
lists use the same ordering rule. N100 and N400 streams are disjoint.

## Exact window and all-displacement evaluation

For `(a,b)`, the code builds the Gaussian quotient with period basis
`(a,b),(-b,a)`. For every nonzero class it obtains the shortest representative
from the four floor/ceiling choices in that orthogonal basis. Ties use
lexicographic `(dx,dy)`. It includes the class exactly when

```text
16*r2 >= N  and  25*r2 <= 4*N,
```

so both endpoints of `1/4 <= distance/sqrt(N) <= 2/5` are retained. Entries
are sorted by `(r2,canonical_dx,canonical_dy)`; z and -z remain distinct. The
program writes `OUTPUT.axis.window.csv` and `OUTPUT.tilted.window.csv` before
occupations. A formal runner must hash these tables and the pinned g16 TSV,
record the hashes, and validate included/rejected classes. This prototype does
not embed a platform-specific SHA256 implementation; the external gate remains
mandatory.

Every configuration exhausts the complete displacement table at each of16
anchors. There is no pair RNG, pair subsampling or distance grid. If B16 is the
signed sum over all anchor/displacement pairs, then

```text
H = B16/(16*A*N),  A=16,
Cbar = B16/(16*A*|D|).
```

The first is unbiased for the full ordered macro-window sum divided by N^2
against every translation-invariant mark. The second is a raw diagnostic from
the same integer, not a separate sample.

## Topology, Bell8 sewing and support controls

Black components use occupied NN adjacency and white components use vacant
matching adjacency. The original marks are

```text
q = black_components - white_components
    - (K - occupied_NN_edges + full_occupied_square_faces),
E = q^2.
```

Virtual colour joins never alter them. Ports are ordered
`(xN,xE,xS,xW,yN,yE,yS,yW)`. N/E incidences create physical edge IDs and S/W
reuse the reverse ID. The general adjacent/shared-edge rule is retained even
though the frozen window asserts disjoint port sets.

Occupied neighbors use their black-NN root; vacant neighbors use the physical
edge ID in a disjoint namespace. First-appearance labels form the24-bit key.
Missing valid keys are exact zero and signed g16 is retained.

For each endpoint-vacant pair the collector counts exterior components meeting
both port groups. It accumulates total, exactly s=2 and s>=3. It requires zero
g16 whenever s<=1 and total=s2+sge3 both configurationwise and in every
batch/K cell. The strata reuse the same data and are correlated coordinates,
not evidence votes.

## Batch-by-K thermal interface

The CSV has one row per `(N,batch,geometry,K)`, including zero-count cells. It
retains `count,sum_q,sum_E` and, for each total/s2/sge3,

```text
sum_B16, sum_qB16, sum_EB16,
eligible_pair_count, nonzero_pair_count.
```

Controls are `s_le1_pair_count`, `s_le1_nonzero_g16_count`, and
`total_minus_s2_minus_sge3_B16`. Metadata columns include both seeds, anchors,
window count and `source_denominator=16*A*N`. Compiler, peak RSS, binary/kernel/
window hashes and commands belong in the external receipt; stdout supplies
elapsed seconds and dimensions.

The K table avoids a one-p Taylor approximation. Reweight every K cell by

```text
(p/p_ref)^K * ((1-p)/(1-p_ref))^(N-K).
```

After each aligned batch deletion, reaggregate the surviving geometries, solve
the pooled root in the frozen bracket and differentiate the finite K sums.
Baseline q/E need p jets through order2; H/qH/EH need order1. These feed all
four centered direct/root/slope terms. S, qS and ES use the same configuration
B16 and are never independently resampled. The ESS and root-bracket gates
remain mandatory; likelihood reweighting does not create unsampled K support.

## Complexity and boundary

Each geometry/configuration costs `O(N alpha(N)+16*|D|*8)`. The dense lookup
uses32 MiB; quotient storage is O(N^2), about0.64 MiB at N400; batch sums use
signed128-bit integers. Anchors and pairs are within-configuration readouts;
the50,000 paired configurations per size are inference units.

The thermal interface is primary because
[`eed2190c`](https://github.com/LightChainr/Matching-One/blob/eed2190c04b67084ab5aef5827e00377853a0bca/notes/p337-critical-spatial-summability.md)
proves raw canonical gxy absolutely summable at exact square-site p_c. Another
raw-C distance grid is not the mechanism test. The K table probes global
thermal likelihood response; it does not itself prove local pivotal support.

This collector emits no root, U, field contrast or decision. The scorer must
enforce frozen ESS, delete-one bracket, covariance and stop gates. The bilocal
window is not an unfiltered homogeneous epsilon derivative and C4 is not a
pure-spin projector. No pilot, production top-up, server action, test grid or
commit was performed here.
