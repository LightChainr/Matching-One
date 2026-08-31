# Regular-pair spatial producer and completed fixed run

`scripts/p337_regular_pair_spatial_sampler.cpp` is the minimal fresh-iid reader of the prescribed two-site regular-pair Q-jet kernel. Contract `3210aeb338ca7bb52c799d1de9048232f50ab921:analysis/p337_regular_pair_spatial_contract.json` fixes `L=32,64`, `p_ref=0.592746050790` and 200 batches of 1000 configurations per size. Seeds are `2026083123593201` and `2026083123596401`, respectively. After acceptance of kernel `32ff99fa` and theorem `7f60e92d`, the root gave GO and these two fixed blocks completed once. Raw results and the complete execution receipt are in `results/p337-regular-pair-spatial/`; no scoring was done by this producer task.

## Fixed geometry and observable

Square periodic tori `L=32,64`; `r=L/4=8,16`. The 16 anchors are `(i*r,j*r)` for `i,j=0..3`. Every anchor contributes the horizontal pair `((x,y),(x+r,y))` followed by the vertical pair `((x,y),(x,y+r))`. Their 32 values are averaged within each iid configuration. They are correlated readouts, never 32 independent samples.

Each configuration builds one occupied NN DSU. Either occupied insertion endpoint forces zero. Otherwise the eight ports are ordered `x:N,E,S,W,y:N,E,S,W`. An occupied neighbor receives its occupied-NN component root; every vacant neighbor's incident edge receives a distinct singleton ID. For these separations all eight physical incident edges are distinct. Canonical restricted-growth labels start at zero, in port order, with key `sum(label[i] << (3*i))`.

The exact lookup is a sparse tab-separated file with named `packed_key` (also accepts `key`) and `g16` columns; additional columns are ignored. Only the 1874 nonzero Bell8 entries are stored out of 4140 partitions, and an absent key is exactly zero by the frozen lookup convention. `g=g16/16` retains its sign, including negative values; no absolute value or clipping is used. No winding/rank, `q/E`, old source, or old data pipeline is included.

## Batched sufficient output

Each CSV row is one batch of configurations. It stores `samples`, `pairs=32*samples`, the integer `sum_g16`, and five integer contributions `sum_g16_shared0` through `sum_g16_shared4`. Shared count is the number of distinct outside components touching both four-port groups. Contributions sum to the total; divide **each by `16*32*samples`**, retaining the full population denominator. Also retained are the five eligible-pair counts, total eligible pairs and nonzero-pair count. Eligibility means both endpoints vacant. The omitted occupied-endpoint pairs have exactly zero kernel; eligibility is not a conditioning of the estimand.

SE/covariance will be calculated across the original batches. The kernel lookup and production contract hashes belong in the root's run receipt. This producer saves an adjacent `OUTPUT.metadata.json` only after completing the CSV, with the exact RNG rule, implemented probability, seed, sample count and elapsed time.

## Reproducible invocation

```sh
clang++ -std=c++17 -O3 scripts/p337_regular_pair_spatial_sampler.cpp -o /private/tmp/p337_regular_pair_spatial_sampler
/private/tmp/p337_regular_pair_spatial_sampler \
  --L 32 --p 0.592746050790 --seed 2026083123593201 \
  --batches 200 --samples-per-batch 1000 \
  --lookup FROZEN_LOOKUP.tsv --output FROZEN_OUTPUT.csv
```

The L64 command changes the size and uses seed `2026083123596401` and its separate output. No defaults silently select the scientific values. The fixed budget is 200 batches of 1000 configurations per size. No L16 sample or benchmark has been run.

The generator is `std::mt19937_64`, one full uint64 word per site in x-fast row-major order; occupancy is `word < 10934234699625173385ULL`. The exact implemented Bernoulli probability is `10934234699625173385/18446744073709551616`, differing from the frozen decimal reference by less than `2^-64`. The CLI accepts only that reference value. N/E/S/W are `(0,+1),(+1,0),(0,-1),(-1,0)`. Batches continue the same size-specific stream without reseeding. Existing outputs are refused. Preparing or syntax-checking this code does not authorize a production run.
