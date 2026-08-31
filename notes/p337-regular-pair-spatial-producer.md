# Regular-pair spatial producer: prepared, not run

`scripts/p337_regular_pair_spatial_sampler.cpp` is the minimal fresh-iid reader of the prescribed two-site regular-pair Q-jet kernel. It has not sampled configurations. The final exact lookup, `p_ref`, two distinct size seeds, budget and root GO are still required before execution.

## Fixed geometry and observable

Square periodic tori `L=16,32`; `r=L/4`. The 16 anchors are `(i*r,j*r)` for `i,j=0..3`. Every anchor contributes the horizontal pair `((x,y),(x+r,y))` followed by the vertical pair `((x,y),(x,y+r))`. Their 32 values are averaged within each iid configuration. They are correlated readouts, never 32 independent samples.

Each configuration builds one occupied NN DSU. Either occupied insertion endpoint forces zero. Otherwise the eight ports are ordered `x:N,E,S,W,y:N,E,S,W`. An occupied neighbor receives its occupied-NN component root; every vacant neighbor's incident edge receives a distinct singleton ID. For these separations all eight physical incident edges are distinct. Canonical restricted-growth labels start at zero, in port order, with key `sum(label[i] << (3*i))`.

The exact lookup is a tab-separated file with named `key` and `g16` columns (additional columns are ignored). A missing encountered key stops the run; it never becomes zero. `g=g16/16`. No winding/rank, `q/E`, old source, or old data pipeline is included.

## Batched sufficient output

Each CSV row is one batch of configurations. It stores `samples`, `pairs=32*samples`, the integer `sum_g16`, and five integer contributions `sum_g16_shared0` through `sum_g16_shared4`. Shared count is the number of distinct outside components touching both four-port groups. Contributions sum to the total; divide **each by `16*32*samples`**, retaining the full population denominator. Also retained are the five eligible-pair counts, total eligible pairs and nonzero-pair count. Eligibility means both endpoints vacant. The omitted occupied-endpoint pairs have exactly zero kernel; eligibility is not a conditioning of the estimand.

SE/covariance will be calculated across the original batches. The kernel lookup and production contract hashes belong in the root's run receipt. This producer saves an adjacent `OUTPUT.metadata.json` only after completing the CSV, with the exact RNG rule, implemented probability, seed, sample count and elapsed time.

## Reproducible invocation, pending frozen values

```sh
clang++ -std=c++17 -O3 scripts/p337_regular_pair_spatial_sampler.cpp -o /private/tmp/p337_regular_pair_spatial_sampler
/private/tmp/p337_regular_pair_spatial_sampler \
  --L 16 --p FROZEN_P --seed FROZEN_L16_SEED \
  --batches FROZEN_BATCHES --samples-per-batch FROZEN_BATCH_SIZE \
  --lookup FROZEN_LOOKUP.tsv --output FROZEN_OUTPUT.csv
```

The L32 command changes the size and uses its separately frozen seed/output. No defaults silently select the scientific values. Candidate budget is 200 batches of 1000 configurations per size, pending the root's final budget decision.

The generator is `std::mt19937_64`, one word per site in x-fast row-major order; occupancy is `(word>>11) < floor(p*2^53)`. Thus the exact implemented Bernoulli probability is the recorded integer threshold divided by `2^53` (within `2^-53` of requested p). N/E/S/W are `(0,+1),(+1,0),(0,-1),(-1,0)`. Batches continue the same size-specific stream without reseeding. Existing outputs are refused. Preparing or syntax-checking this code does not authorize a production run.
