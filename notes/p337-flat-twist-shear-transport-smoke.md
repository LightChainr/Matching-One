# F3 flat-twist shear transport: exact gate, not a new experiment

Status: fresh local N65 low-sample smoke after the representation matrix was
frozen at `3cb06ab`.  No Huawei resource was used and no larger run is
recommended from this result.

## Archive audit

The repository has one production-format projective-birth archive before this
smoke: P334 N65 at `1714141`.  It pairs the physically distinct Gaussian
orientations `(8,1)` and `(7,4)`.  It does not contain identity and period-basis
shear copies of either shape, so it cannot provide an end-to-end archive test
of the `3cb06ab` transport convention.

The new smoke uses 10,000 fresh counter permutations, 20 aligned batches,
seed `202608337`, counter range `[3370000,3380000)`.  Each of the two Gaussian
shapes is run as the pair

```text
P_identity,
P_shear = P_identity T^-1,  T=[1 1;0 1].
```

Explicitly,

```text
(8,1): [[8,-1],[1,8]] -> [[8,-9],[1,7]],
(7,4): [[7,-4],[4,7]] -> [[7,-11],[4,3]].
```

Right multiplication by `T^-1` leaves the period lattice, canonical HNF and
finite graph unchanged.  It only changes the period basis.  Consequently the
projective birth line must transform pathwise as

```text
ell_shear = primitive(T ell_identity).
```

## Result

Every sparse `(tau1,tau2,kind,ell,count)` cell in all 40 shape/batch blocks
obeys the exact line transport.  The fresh same-counter orientation contrast
is

```text
C_identity = [H,A,D]
           = [-0.0056598827,-0.0004602931,-0.0002831596].
```

The matrix frozen at `3cb06ab` gives

```text
M_T C_identity
 = [-0.0001252523,-0.0043738678,-0.0036304151],
```

while the directly reconstructed shear archive gives

```text
C_shear
 = [-0.0001252523,-0.0043738678,-0.0036304151].
```

The full three-dimensional mean residual is

```text
[+4.88e-19,-8.79e-18,-7.98e-18],
```

and the maximum absolute batch residual is `6.55e-17`.  The JSON retains the
complete 6x6 covariance of `[C_identity,C_shear]` and the 3x3 residual
covariance.  A residual chi-square is deliberately undefined: the residual is
a pathwise identity, so its covariance has rank zero up to floating roundoff.

The fresh identity mean differs from the old 20k exploratory mean, as expected
at this sample size.  It was not used to change the frozen matrix.  The object
tested here is the parameter-free transport relation, not equality of two
Monte Carlo point estimates from disjoint counters.

## Scientific consequence

This closes the implementation/convention chain

```text
period basis -> primitive line -> F3 bins -> [H,A,D] -> A4 standard action.
```

It does **not** add evidence for a charged continuum field.  A modular
right-basis shear is a gauge-like relabeling of the same quotient graph, so
more samples of this pair can only repeat an exact identity.  The nontrivial
next discriminator must alter the source while retaining the covariance
contract—for example, an explicit finite-field twist/defect insertion whose
identity and shear versions are not obtained by deterministic relabeling.
Only that experiment can test whether a physical charged response follows the
same A4 transport.

## Reproduce

```bash
clang++ -O3 -std=c++17 src/threshold_rank_integer_period_mc.cpp \
  -o build/p337-shear-smoke/threshold_rank_integer_period_mc

build/p337-shear-smoke/threshold_rank_integer_period_mc \
  --samples 10000 --batches 20 --seed 202608337 \
  --replica-offset 3370000 --threads 1 \
  --first-matrix 8 -1 1 8 --second-matrix 8 -9 1 7 \
  --first-rep 8 1 --second-rep 8 1 --projective-births \
  --git-commit 3cb06ab \
  --output-prefix results/local-20260830/P337-flat-twist-shear-smoke/shape_a

build/p337-shear-smoke/threshold_rank_integer_period_mc \
  --samples 10000 --batches 20 --seed 202608337 \
  --replica-offset 3370000 --threads 1 \
  --first-matrix 7 -4 4 7 --second-matrix 7 -11 4 3 \
  --first-rep 7 4 --second-rep 7 4 --projective-births \
  --git-commit 3cb06ab \
  --output-prefix results/local-20260830/P337-flat-twist-shear-smoke/shape_b

python3 scripts/score_flat_twist_shear_transport.py \
  --shape-a-births results/local-20260830/P337-flat-twist-shear-smoke/shape_a.births.csv \
  --shape-a-metadata results/local-20260830/P337-flat-twist-shear-smoke/shape_a.metadata.json \
  --shape-b-births results/local-20260830/P337-flat-twist-shear-smoke/shape_b.births.csv \
  --shape-b-metadata results/local-20260830/P337-flat-twist-shear-smoke/shape_b.metadata.json \
  --json results/local-20260830/P337-flat-twist-shear-smoke/score.json \
  --markdown results/local-20260830/P337-flat-twist-shear-smoke/REPORT.md
```
