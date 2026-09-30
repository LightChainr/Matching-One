# Joint births: an executed finite-size pilot, not another readiness exercise

Date: 2026-09-29. Claim level: new exploratory independent Monte Carlo block.

**The useful result is a separation of scales.** Single-insertion double births
become rare, but a fixed small near-critical neighbourhood of the diagonal
still contains substantial paired-birth mass. The new soft statistic behaves
approximately linearly in its resolution on this pilot's grid, rather than
showing a resolved positive intercept. This supports prioritising quantitative
near-diagonal control; it does not prove that the scaling limit has no double
births.

## What actually ran

- NN square-site and triangular-site tori, L=16,32,64,128; N=L² sites.
- Eight independent pseudorandom batches × 2,000 uniform site permutations per
  lattice/size: **128,000 production filtrations** in total.
- Both births are recorded in the **same** monotone filtration:
  `J1=first rank>=1`, `J2=first rank=2`, `D=J2-J1`.
- The rank is the integer winding-image rank of the union of all occupied
  components. It is not two directional wrap Booleans. A diagonal rank-one
  cycle is still rank one.
- Triangular edges add only the consistent `(1,1),(-1,-1)` diagonal to the
  four square neighbours. This is triangular site percolation on a primitive
  rhombic torus, **not** square matching-lattice percolation with both diagonals.
- Local ARM64, Apple clang 21, C++17 and Python 3.13.7 standard library;
  maximum four simultaneous processes. Acquisition including compilation,
  control and benchmark took **6.89 seconds**. Sum of individual production
  process wall times was 17.43 seconds; that is not measured CPU time.
- Before production, all 9! permutations were enumerated once for each L=3
  lattice. Square reproduced `P(D=0)=3/35`, `E D=3/2`, `E D²=43/14`;
  triangular reproduced `2/35`, `3/2`, `81/28`. The largest requested cell was
  then benchmarked with 200 separate samples per lattice. Control/benchmark
  samples are excluded from production estimates.

The displacement convention is adapted from
[`src/threshold_rank_axis_mc.cpp`](../src/threshold_rank_axis_mc.cpp).
The new engine accumulates a global basis from every non-tree cycle, using
the winding-vector determinant. All edges of an insertion are processed before
checking rank, so an actual `0->2` step is not artificially split by edge order.

## Results that change the next question

| Lattice | L | Direct atoms / 16,000 | E[D]/L^(5/4) | Z(.025), width-normalised labels | Z(.05), width-normalised labels |
|---|---:|---:|---:|---:|---:|
| square | 16 | 97 | .4230 ± .0032 | .01815 ± .00123 | .02996 ± .00217 |
| square | 32 | 24 | .4206 ± .0029 | .01389 ± .00055 | .02616 ± .00089 |
| square | 64 | 8 | .4273 ± .0026 | .01402 ± .00052 | .02684 ± .00111 |
| square | 128 | 1 | .4278 ± .0030 | .01365 ± .00054 | .02718 ± .00080 |
| triangular | 16 | 38 | .4032 ± .0018 | .01439 ± .00043 | .02634 ± .00065 |
| triangular | 32 | 13 | .4031 ± .0027 | .01372 ± .00031 | .02665 ± .00049 |
| triangular | 64 | 4 | .4075 ± .0030 | .01360 ± .00079 | .02709 ± .00127 |
| triangular | 128 | 1 | .4056 ± .0031 | .01435 ± .00103 | .02813 ± .00137 |

Errors are eight-batch delete-one SE, not confidence limits. Here
`W=Q.75-Q.25` of the equal-weight empirical J1/J2 mixture,
`G=(N+1)(T2-T1)/W`, and `Z(delta)=E[(1-G/delta)_+]`. W is re-estimated in each
delete-one sample. Its discrete quantile can jump at small L; the unusually
large L=16 square width-normalised uncertainties are retained, not smoothed
away or interpreted as a new physical effect.

At L=128, the exponent-free label-gap CDF at `.025,.05,.1,.2` is

```text
square:     .02728, .05426, .10843, .20801
triangular: .02832, .05534, .10771, .20563
```

The corresponding Z values are

```text
square:     .01365, .02718, .05440, .10634
triangular: .01435, .02813, .05487, .10610
```

These values are compatible with a finite near-origin gap density over the
resolved range: halving the resolution roughly halves Z. This is a descriptive
finite-grid observation, **not** a fitted exponent, a hypothesis test, or an
exclusion of an unresolved smaller-scale diagonal mass. The two lattices have
different conformal torus shapes; their numerical proximity after width
normalisation is not a same-modulus universality test.

One observed direct atom out of 16,000 means frequency `0.0000625`, with
two-sided 95% Wilson interval approximately `[0.0000110,0.0003540]`. Even zero
observations would only give a one-sided 95% upper bound `0.0001872`; neither
would be a theorem of zero probability.

The full table also reports the prescribed `L^(5/4)` insertion scale and
`L^(3/4)` label scale. They are useful fixed reference scales, not an established
square-site critical-exponent theorem. The width-normalised readout avoids
assuming that exponent. No free power, lattice metric factor or cutoff was
fitted to improve agreement.

## Exact use of the paired histogram

Uniform random ordering is independent of the order statistics of iid uniform
labels. Conditional on `D=d>0`, `T2-T1 ~ Beta(d,N+1-d)`; d=0 is an atom.
The analyzer integrates this conditional law exactly up to floating arithmetic,
using integer-parameter beta/binomial identities. No extra label Monte Carlo
noise is added. Thus the same histogram supplies:

1. the direct atom and discrete near-diagonal CDF;
2. continuous-label near-diagonal CDF at the same four resolutions;
3. both discrete and label-integrated Z(delta);
4. birth correlation, gap mean and self-normalised widths.

Z was added after acquisition, following the theory agent's exact statistic
request; its grid is the original `.025,.05,.1,.2`. It is not an independently
preregistered second experiment. All these coordinates share one dependency
group per lattice/size/batch. Complete delete-one covariance and all deleted
estimates are saved, not treated as independent evidence votes.

The link to a genuine noncoalescence target is proved in
[`birth-gap-noncoalescence-20260929.md`](birth-gap-noncoalescence-20260929.md):
under pair tightness, every limit has distinct births iff
`lim_delta->0 limsup_L Z_L(delta)=0`. The pilot does not supply pair tightness
or this double limit. A stable mean gap alone would not exclude a partial
diagonal atom either.

## A useful independent next block, if the parent allocates it

Do not use a decreasing direct atom as the success criterion. Two competing
behaviours to distinguish are a stable intercept of Z at shrinking resolution
versus near-origin continuous mass with Z decreasing with resolution. The
smallest useful extension is independent L=128/256 (512 if its benchmark is
cheap), at a declared finer width-normalised grid such as
`.00625,.0125,.025,.05`. Repeated large samples at L=16 are less relevant to
that distinction. This is a suggestion, not an executed cloud campaign or
a requirement that other research wait for it.

For an independent run, use a **new random-seed namespace**, keep its batches
separate, and carry both CDF and Z under the same paired covariance. Reject
the continuous-only interpretation at the examined scales if a clearly
resolved stable intercept develops; retain the unresolved range if precision
or size drift prevents a distinction. Do not invent a numerical rejection
threshold after seeing that new block, and do not claim a finite intercept
distinguishes all possible logarithmically slow approaches to zero.

## Reproduction and handoff

The exact executed command was

```bash
python3 analysis/birth-gap-20260929/run_pilot.py \
  --sizes 16 32 64 128 --batches 8 --samples 2000 --workers 4
python3 analysis/birth-gap-20260929/analyze.py
```

The acquisition runner refuses to overwrite an existing run; for reproduction
choose `--output analysis/birth-gap-20260929/reproduction-data` and pass that
directory to the analyzer with `--data`, plus distinct output paths. Repeating
the original seeds reproduces the block; it is not new independent data.

Files:

- [`engine.cpp`](../analysis/birth-gap-20260929/engine.cpp): transparent C++17 engine;
- [`run_pilot.py`](../analysis/birth-gap-20260929/run_pilot.py): bounded local
  acquisition (this run used four); subsequently extended for an explicitly
  named independent cloud block with a portable compiler and up to 14 workers;
- [`analyze.py`](../analysis/birth-gap-20260929/analyze.py): stdlib readout;
- [`data/run.json`](../analysis/birth-gap-20260929/data/run.json): per-batch
  seeds, commands, timings, compiler, base commit, source/data SHA256;
- [`SUMMARY.md`](../analysis/birth-gap-20260929/SUMMARY.md) and
  [`summary.json`](../analysis/birth-gap-20260929/summary.json): full CDF/Z grid,
  raw counts, uncertainty and covariance;
- `data/pilot-*.json.gz`: full sparse `(J1,J2,count)` histogram for every
  batch; all production data together are below one megabyte.

For ARM Linux the engine is directly portable, with no OpenMP dependency:

```bash
g++ -O3 -std=c++17 -DNDEBUG analysis/birth-gap-20260929/engine.cpp -o /task/path/birth-gap
/task/path/birth-gap square 256 10000 2026092900256000 /task/path/square-L256-b00.json sample
```

The example seed is a **new** cloud block seed, not a local batch seed. The
extended runner now accepts `--compiler g++ --role NEW_NAMESPACE --workers 14`,
`--skip-exact-control --exact-control-reference ...`, and an explicit
`--git-base-commit` for copied source trees. See
[`CLOUD-HANDOFF.md`](../analysis/birth-gap-20260929/CLOUD-HANDOFF.md) for the full
bounded invocation and finer prospective epsilon grid. The local manifest and
generated results were not rewritten after extending the runner; its saved
runner SHA records the acquisition-time revision, while the engine is unchanged.
No local Monte Carlo was rerun. No cloud access, commits, pushes or navigation
edits were performed by this pilot subtask.
