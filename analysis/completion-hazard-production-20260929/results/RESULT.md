# Completion-hazard production: fixed-cohort readout

Square L=512; N=262144; fixed insertion cutoffs a=154646, b=155385, c=156120.

Scored 140,000 filtrations in 14 batches. This is one shared random block, not independent evidence for each column.

## Primary estimates

| Estimand | Estimate | Aligned batch SE | Status |
|---|---:|---:|---|
| Pooled survival: early minus late | 0.006857066 | 0.004444348 | scoreable |
| Integrated-hazard prediction: -I early + I late | 0.01550219 | 0.003505008 | scoreable |
| Survival difference minus integrated prediction | -0.008645121 | 0.006389162 | scoreable |
| Instantaneous nu difference, window-scaled (not integrated) | -0.01504115 | 0.003301171 | scoreable |
| Exact-D overlap-weighted survival difference | -0.004149993 | 0.004511942 | scoreable |
| Exact-(D,nu_b) overlap-weighted survival difference | 0.001216127 | 0.004763563 | scoreable |
| Early map residual: mean Z + mean I - 1 | -0.009410927 | 0.006314693 | scoreable |
| Late map residual: mean Z + mean I - 1 | -0.0007658063 | 0.004128975 | scoreable |

The nu contrast has the hazard sign (early minus late) and is only instantaneous. It is not an integrated survival prediction.

## Initial cohort means

| Cohort | Count | mean Z | mean I | mean nu_b |
|---|---:|---:|---:|---:|
| early | 24218 | 0.4426873 | 0.5479018 | 99.40206 |
| late | 30505 | 0.4358302 | 0.563404 | 101.5868 |

## Exact common-support conditional estimands

| Stratification | Common cells | Supported E / L | Dropped E | Dropped L | Dropped total |
|---|---:|---:|---:|---:|---:|
| D | 6 | 24218 / 30490 | 0 | 0.0004917227 | 0.0002741078 |
| (D, integer nu_b) | 1202 | 24105 / 30015 | 0.004665951 | 0.01606294 | 0.01101913 |

Weights are e*l/(e+l) in each exact shared cell. These targets differ from the pooled difference and from each other. No mediated/explained fraction is defined. Missing support is reported, not imputed or binned.

## Direction-resolved descriptive means

| Primitive unoriented D | E n | L n | E Z | L Z | E I | L I | E nu_b | L nu_b |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| (0, 1) | 11168 | 13113 | 0.4595272 | 0.4607641 | 0.5356463 | 0.5455196 | 94.1918 | 93.52048 |
| (1, -2) | 2 | 4 | 0 | 0 | 0.5654377 | 0.6576087 | 197 | 353.75 |
| (1, -1) | 925 | 2072 | 0.2648649 | 0.2890927 | 0.6603045 | 0.6938002 | 155.0865 | 155.2375 |
| (1, 0) | 11184 | 13150 | 0.4544886 | 0.4587072 | 0.5390584 | 0.5393133 | 95.34505 | 92.65224 |
| (1, 1) | 938 | 2148 | 0.2782516 | 0.2886406 | 0.6889624 | 0.6927524 | 154.4733 | 151.9046 |
| (1, 2) | 0 | 10 | not_scoreable | 0.1 | not_scoreable | 0.5775421 | not_scoreable | 302.4 |
| (2, -1) | 1 | 3 | 0 | 0.3333333 | 0 | 0.3531556 | 301 | 186.3333 |
| (2, 1) | 0 | 5 | not_scoreable | 0 | not_scoreable | 1.243709 | not_scoreable | 250.6 |

## Interpretation and uncertainty

- E[Z+I | initial cohort]=1 is an expectation identity, not a pathwise equality. Closure assesses the implemented geometry-to-exit map; it does not establish nu_b closure.
- The (D,nu_b) finite-lag contrast is a distinct necessary prediction of a state-sufficiency model. A zero weighted contrast can hide opposite cell effects and does not prove recursive Markovness.
- All delete-one replicates remove the same whole batch for every metric. All stratum weights and shared-cell support are recomputed. Full covariance and replicate vectors are in result.json; no pairwise subset covariance is used.
- Per-direction rows are descriptive. No claim about causal mediation, the iid-label clock, or a large-size non-Markov limit is made.

## Input provenance

Each listed batch was read once, checked against its SHA256 and expected sample count, then retained only as sufficient statistics for resampling. Checksums do not prove tau independence or independent RNG design; those are generator-contract assumptions.

| Batch | Samples | Seed | SHA256 |
|---|---:|---|---|
| production-square-L512-b00.csv.gz | 10000 | 18419746559907930925 | 32f16cc37728d67fdbe3446a48da1bb4d36c48fc4f051648ac8e36e2595d4026 |
| production-square-L512-b01.csv.gz | 10000 | 13442365983088886094 | fdc708228c1e0535787aa3426ea725b276ecef0b8ba193215b1bf93512518892 |
| production-square-L512-b02.csv.gz | 10000 | 12192384056273429425 | 57b5f9ef9bcc6f21a5a14f5e7b6c24fccb94388f349943738916cc1700ab4897 |
| production-square-L512-b03.csv.gz | 10000 | 16795263932832913938 | 65ec5cee09f76a03c69b702a21bb6f27dd70930b48e3b7f047de33103032eaa4 |
| production-square-L512-b04.csv.gz | 10000 | 5606540329541986450 | 033189fab591caa248c66c9710d63104036dc3d26b88e495d965615a84ecc4fb |
| production-square-L512-b05.csv.gz | 10000 | 176734671850227651 | b38b95839aab1c6f6575a9b94702c41a032a2cc00430c6eb904840472488fbd4 |
| production-square-L512-b06.csv.gz | 10000 | 14922363405798405830 | 330829d82ec4dad029f53dc07fe8bda5252bb35118aa9b507012c73df88e5a35 |
| production-square-L512-b07.csv.gz | 10000 | 3540279280158270918 | a292a7a8be0fb5f79f2449745709c51942fc4261fd50e542569d33d04005ee30 |
| production-square-L512-b08.csv.gz | 10000 | 182843979474712324 | 113325911eb95469467bfe73496d5b57448413125dcea383b55b0ff7875db4ac |
| production-square-L512-b09.csv.gz | 10000 | 4164285986391024372 | b5712154670f8179ddf32487f35a794ea16b5ef4f0cdad2025c3f8e447f321da |
| production-square-L512-b10.csv.gz | 10000 | 15500924042285991404 | 98bb914abc5188216cd5b9b4584ba46539bd4c16c1b4d745b6a8033748243668 |
| production-square-L512-b11.csv.gz | 10000 | 17186215622830953477 | 73aaa04f05ee279532f7c781a97acda5191840cf79fc2cee72f2a6c9caaed7f1 |
| production-square-L512-b12.csv.gz | 10000 | 11279241988091952515 | 75f556befdd4986e26f6594f7d2ec44fa6517b9b19021b9cbaa612116471875c |
| production-square-L512-b13.csv.gz | 10000 | 10329155733035115099 | d5ffad750a152bd62325ae22375496b861da96f8687ded8ad7c0a8220d06a202 |

No old archive was pooled into these estimates.
