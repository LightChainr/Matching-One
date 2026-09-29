# Independent cloud block: ready-to-run handoff

Run from the copied task repository root, using only the allocated TV2N0X
task directory. This subtask did not access or run on the server.

The local engine has already passed all 9! L=3 permutation controls once per
lattice. Its unchanged SHA256 is
`72997c726bc734da7695d104fc8fb51ed5672506c967c2a88312bdc092925506`.
The completed local manifest SHA256 is
`d223ad486d153239a40c02fb534e03efe33fd9b5dde741f206b32d147c608f4e`.

Acquisition (14 measured-available workers, 140,000 paths per cell):

```bash
python3 analysis/birth-gap-20260929/run_pilot.py \
  --compiler g++ \
  --sizes 64 128 256 512 --batches 14 --samples 10000 --workers 14 \
  --role cloud-tv2n0x-20260929-v1 \
  --skip-exact-control \
  --exact-control-reference 'local L3 square+triangular all 9!; data/run.json sha256=d223ad486d153239a40c02fb534e03efe33fd9b5dde741f206b32d147c608f4e' \
  --git-base-commit d31fa5fe9584b77595e4a78b9560aa2ff77cf0b3 \
  --budget-seconds 300 \
  --output analysis/birth-gap-20260929/cloud-data
```

This first compiles with the specified compiler, skips only the already-supplied
exact control, benchmarks 200 permutations per lattice at the largest requested
L, and stops **before production** if the size-scaled benchmark predicts more
than 300 seconds at the requested worker count. The estimate is not an enforced
hard timeout; CPU contention can make actual time longer. Per-batch data and the
partial manifest are saved as jobs finish.

If the benchmark stops, choose L=64,128,256 and a **different output directory**
(the existing manifest is never overwritten). Use a new role if a production
block already began; otherwise the original production seeds have not been used.
The production namespace is hashed with lattice/L/batch. Different namespace
means separate pseudorandom streams; changing output path alone does not.

Analyze the cloud block without changing local results:

```bash
python3 analysis/birth-gap-20260929/analyze.py \
  --data analysis/birth-gap-20260929/cloud-data \
  --epsilons .00625 .0125 .025 .05 .1 .2 \
  --output-json analysis/birth-gap-20260929/cloud-summary.json \
  --output-md analysis/birth-gap-20260929/CLOUD-SUMMARY.md
```

This epsilon grid is declared before cloud results. It applies to both hard
near-diagonal CDF and soft Z(delta), in prescribed-power and empirical-width
clocks. Source compilation and analysis use only C++17/Python stdlib; no scipy,
NumPy, OpenMP or installations beyond the parent's compiler provisioning.

The analyzer saves full paired delete-one covariance, re-estimates empirical
birth widths within every deletion, and explicitly marks reused prior exact
controls rather than claiming they ran again on the cloud. The optimized beta
recurrence terminates only its floating-point-underflowed far right tail; at
L<=512 and epsilon<=.2 the starting binomial mass is representable. It raises
an error instead of silently substituting a normal approximation outside that
regime.

Interpretation remains finite-size: a decreasing microscopic atom and a stable
mean gap do not exclude a partial diagonal atom. The cloud block adds an
independent larger-size test of the soft near-diagonal behaviour, not a proof
of square-site exponents or process-limit noncoalescence.
