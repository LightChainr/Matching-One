# Fixed-root landing-marked pivotal H4 prototype

## Outcome

The bounded prototype succeeds as an engineering bridge. It records a local
four-arm landing mark that is absent from threshold-rank histograms, while
retaining the ordinary total pivotal mass only as a Russo regression control.

The exact axis `L=4`, `R=1`, `p=1/2` oracle enumerated all `2^15` non-root
configurations. Fixed-root mass times 16 equals the independent all-site exact
pivotal mass for both terms:

```text
primal:   1.52392578125, difference 0
matching: 2.58642578125, difference 0
```

Configuration-level C4 rotation, reflection, and 45-degree landing-registry
controls all have zero violations. `R=1` is only a schema and symmetry oracle,
not a physical scaling point.

## N=65 pilot

The paired Gaussian representations `(8,1)` and `(7,4)` used the same 200,000
counter-keyed non-root Bernoulli fields, split into 100 aligned batches. The
counter interval `[12000000000,12000200000)` is independent of the existing
full-curve runs. The frozen local radius is `R=3`; the generation code commit
is `f9bc73a7cbbdc2b52fbda0167bc8ccc0eecb73ae`.

| orientation | mu4 | delete-one SE | a4 | delete-one SE | landing acceptance |
|---|---:|---:|---:|---:|---:|
| `(8,1)` | 2.566200 | 0.035355 | 0.492976 | 0.005710 | 0.621079 |
| `(7,4)` | 1.884025 | 0.036136 | 0.370700 | 0.006590 | 0.608246 |

The paired first-minus-second differences are

```text
mu4: 0.682175 +/- 0.050675
a4:  0.1222766 +/- 0.0089797
```

There are 16,017 and 15,638 landed pivotal events, respectively, so the local
counter has adequate acceptance for a future preregistered size study. These
numbers are not an exponent test. A single `N` and a frozen square annulus
cannot distinguish a continuum four-arm coefficient from finite-radius lattice
response.

The unmarked controls are `mu0=8.381425 +/- 0.046304` and
`8.355750 +/- 0.052570`. They are derivative/Russo checks, not a second result.
`mu4`, `a4`, `mu0`, and both orientations share configurations and are retained
with their full delete-one covariance.

## Reproduction

```bash
mkdir -p build/marked-pivotal-h4
c++ -std=c++17 -O3 src/marked_pivotal_h4_mc.cpp \
  -o build/marked-pivotal-h4/marked_pivotal_h4_mc

python3 scripts/marked_pivotal_h4_reference.py --exact-axis-l4 \
  --output results/local-20260829/P100-marked-pivotal-h4/analysis/exact-axis-l4-r1.json

build/marked-pivotal-h4/marked_pivotal_h4_mc \
  --samples 200000 --batches 100 --threads 1 \
  --seed 2026106201 --replica-offset 12000000000 \
  --git-commit "$(git rev-parse HEAD)" \
  --radius 3 --p 0.592746050790 \
  --output-prefix results/local-20260829/P100-marked-pivotal-h4/raw/n65_r3_200k

python3 scripts/analyze_marked_pivotal_h4.py \
  --batches results/local-20260829/P100-marked-pivotal-h4/raw/n65_r3_200k.batches.csv \
  --metadata results/local-20260829/P100-marked-pivotal-h4/raw/n65_r3_200k.metadata.json \
  --exact results/local-20260829/P100-marked-pivotal-h4/analysis/exact-axis-l4-r1.json \
  --output results/local-20260829/P100-marked-pivotal-h4/analysis/score.json
```

Add `-fopenmp` and choose `--threads` when the local compiler provides OpenMP.
