# Periodic-torus sufficient statistics for the Camia--Feng lattice log pair

Status: executable Phase-A0 protocol and exact algebra oracle for Issue #234.
No universal logarithmic coefficient is claimed or fitted here.

Production entry point: `scripts/run_triangular_energy_logpair_mc.py`.

## Fields and the limit that must be preserved

For critical site percolation on the triangular lattice, Camia--Feng assign an
independent symmetric sign `sigma_C` to each black cluster and define

```text
S_x = sigma_C  if x belongs to black cluster C,
S_x = 0        if x is white.
```

For endpoints on one lattice axis,

```text
E_a(z)       = S_(z-a) S_(z+a) - <S_(z-a) S_(z+a)>,
E_a^delta(z) = S_(z-delta) S_(z+delta)
               - <S_(z-delta) S_(z+delta)>.
```

Their normalized fields are

```text
phi_a       = a^(-5/4) |log a|^(-1) E_a,
eta_a^delta = pi_a^(-2) E_a^delta,

hat(phi)_a^delta
  = (2 delta)^(-25/24) eta_a^delta
    + kappa log(2 delta) phi_a.
```

The paper's convenient normalization is `kappa=C1 CL/C2`, but the log partner
is defined only up to adding the bottom field.  This protocol therefore does
not insert a numerical `kappa`: `C1`, `C2`, `CL`, and the one-arm convention
must first be matched in the same periodic geometry.

The limits are noncommuting and are frozen as

```text
for each fixed physical delta: a=1/L -> 0;
only after that:              delta -> 0.
```

A simultaneous free fit in `a` and `delta` is not a realization of the
Camia--Feng construction.

Primary source: Camia and Feng,
*The percolation energy field and its logarithmic partner*,
arXiv:2508.16047v2, equations (1.2), (1.8), (1.15), and (1.18)--(1.23).

## Cluster signs can be integrated out exactly

Fix a percolation configuration `omega`.  For an endpoint pair `A=(x,y)`,
write `X_A=S_x S_y`.  Independence and symmetry of the cluster signs give

```text
u_A(omega) := E_sigma[X_A | omega]
            = 1{x and y are in the same black cluster}.
```

For two endpoint pairs `A` and `B`, define

```text
h_AB(omega) := E_sigma[X_A X_B | omega].
```

This is another zero/one observable.  It equals one exactly when all four
endpoint slots are black and every black cluster touched by those slots has
even multiplicity.  Equivalently, the allowed endpoint partitions are the
all-four block or a union of two two-blocks.  Counting slots, rather than
distinct sites, makes the identity valid even if a tiny periodic quotient
identifies endpoints.

Thus every cluster-sign draw can be removed.  The estimator is a strict
Rao--Blackwellization of an explicit-sign estimator and needs only the black
cluster labels already produced by a union--find pass.  Production can reuse
the triangular neighbor convention in `scripts/square_bond_kappa3.py`; no
new connectivity or winding engine is required.

## Unbiased centering on one stream

The desired covariance is

```text
<E_A E_B> = E_omega[h_AB] - E_omega[u_A] E_omega[u_B].
```

Subtracting two sample means from the same configurations introduces an
`O(1/n)` bias.  For `n` independent configurations, the exact unbiased
second-order U-statistic is

```text
C_AB = (sum_i h_AB,i)/n
       - [(sum_i u_A,i)(sum_i u_B,i)-sum_i u_A,i u_B,i]
         / [n(n-1)].
```

On a periodic torus, the primary protocol averages over all `V=L^2`
translations before the configuration enters the U-statistic.  To keep the
archive integral, record translation sums `U_A,i` and `H_AB,i`.  The formula
then becomes

```text
C_AB = (sum_i H_AB,i)/(n V)
       - [(sum_i U_A,i)(sum_i U_B,i)-sum_i U_A,i U_B,i]
         / [n(n-1)V^2].
```

One block must retain only these integer sufficient statistics:

```text
n,
sum U_L1, sum U_L2, sum U_D1, sum U_D2,
sum U_L1 U_L2, sum U_L1 U_D2, sum U_D1 U_L2, sum U_D1 U_D2,
sum H_LL, sum H_L1D2, sum H_D1L2, sum H_DD.
```

They produce, in frozen order,

```text
C_LL = <E_a(z1) E_a(z2)>,
C_LD = 1/2 [<E_a(z1) E_a^delta(z2)>
            +<E_a^delta(z1) E_a(z2)>],
C_DD = <E_a^delta(z1) E_a^delta(z2)>.
```

Compute this three-vector independently in at least 50 random-stream blocks.
The ordinary covariance-of-the-mean across block vectors is the required full
3 by 3 covariance; it automatically retains the strong common-stream
correlations among LL, LD, and DD.

## Frozen periodic geometry

Use the triangular basis with periods `(L,0),(0,L)` and rescale the physical
torus so `a=1/L`.  The primary local axis is `e1`; the two centers differ by
`(L/2,L/2)`, and `L` is even.  For fixed physical `delta`, the bilocal radius
is the unique nearest integer to `delta L`.

The protocol uses irrational fixed cutoffs

```text
delta = 1/(8 sqrt(2)), 1/(12 sqrt(2)), 1/(16 sqrt(2)),
L     = 64, 96, 128, 192 at every fixed delta.
```

This prevents `delta/a` from being an integer and prevents nearest-vertex
ties, exactly matching the paper's cutoff condition.  It also makes the
realized endpoint displacement differ from the declared physical cutoff by
only `O(a)`.

## Tiny exact oracle

`scripts/triangular_energy_logpair_stats.py` exhausts all `2^18=262144`
states of a 6 by 3 periodic triangular quotient at `p=1/2`.  It explicitly
enumerates the signs of every endpoint-touched cluster to check the parity
formula 1,048,576 times.  The exact centered `(LL,LD,DD)` vector is

```text
584500095/68719476736,
714261375/68719476736,
584500095/68719476736.
```

These numbers validate the algebra, periodic indexing, and centering
convention only.  The quotient intentionally has non-scaling radii and the
values are not continuum amplitudes.

## Runnable Monte Carlo and local smoke

The Phase-A runner assigns every batch its own SplitMix-derived stream, so
worker scheduling does not affect the archive.  Each configuration performs
one black-cluster union--find pass, evaluates the translation-summed `U/H`
statistics, and never draws a cluster sign.  It writes both a CSV containing
the integer sufficient statistics and a JSON analysis containing

```text
(LL,LD,DD), their full 3 by 3 covariance,
the two-field matrix [[LL,LD],[LD,DD]],
Delta_J = LL*DD-LD^2,
J = LL*DD/LD^2 - 1.
```

`J` is unchanged by separate nonzero rescalings of the local and bilocal
fields.  It is therefore a useful early Jordan-shape diagnostic before
`pi_a` and `kappa` are known.  A logarithmic two-point Gram form with a
vanishing bottom-bottom entry and nonzero mixed entry has negative determinant
and `J -> -1`; this is a discriminator to measure, not a result assumed by
the sampler.

The committed 40-sample `L=32`, `delta=1/(8 sqrt(2))` smoke produced

```text
(LL,LD,DD) = (0.0003814, 0.0009275, 0.0011047),
J = -0.510 +/- 0.464.
```

Its only role is to show that sampling, integer archival, U-stat centering,
covariance, and nonlinear diagnostics run end to end.  Forty samples have no
scientific evidentiary value.

The first Huawei line is frozen in
`experiments/p234_phaseA_huawei_20260829.yaml`.  For example:

```bash
python3 scripts/run_triangular_energy_logpair_mc.py run \
  --L 64 --delta-denominator 8 --samples 100000 --batches 100 \
  --seed 2026234064 --workers 16 \
  --output-prefix results/server-20260829/P234-phaseA/d8-L64-100k
```

Run the analogous `L=96,128,192` rows before moving to the smaller fixed
cutoffs.  This preserves the required fixed-`delta` first limit.

## A projective Jordan pencil across sizes

### Fixed-delta continuum score from the existing archive

The batch sufficient statistics already retain `sum_D1` and `sum_D2`.
Dividing by the sample count and torus volume gives the probability that the
two endpoints of a bilocal insertion belong to one black cluster.  The
Camia--Feng spin two-point limit implies, at one fixed physical `delta`,

```text
pi_a^(-2) p_conn(a,delta_realized)
  = constant * (2 delta_realized)^(-5/24) * (1+o(1)).
```

The unknown factor is independent of `L`.  Consequently the already archived
data recover the relative `pi_a` normalization across sizes without another
Monte Carlo run.  With

```text
alpha_L = L^(5/4)/log L,
```

the primary fixed-delta continuum vector is, up to two harmless constants,

```text
[alpha_L^2 LL,
 alpha_L beta_L LD/p_conn,
 beta_L^2 DD/p_conn^2],

beta_L=(2 delta_declared)^(-25/24)
       (2 delta_realized)^(-5/24).
```

`scripts/score_p234_fixed_delta_continuum.py` propagates the joint batch
covariance of `(LL,LD,DD,p_conn)` and scores the direct parent-pair target with
one leading analytic `1/L` correction: the first coordinate extrapolates to
zero, while the second and third extrapolate to constants.
This is the correct `a -> 0` question at fixed physical geometry.

A Jordan shear under a dilation of continuum coordinates must not be confused
with a `log L` drift as the ultraviolet mesh is removed after the fields have
already received their `a`-dependent normalization.  The projective pencil
below remains useful for a calibrated change of physical scale/cutoff, but is
not the primary fixed-delta production score.

Once the relative field gauge is fixed, there is a stronger projective
fingerprint than scoring `J` one size at a time.  In a canonical logarithmic
basis the fixed-delta two-point matrix
has the projective form

```text
M(t)/LD = [[0,1],[1,d0+k t]],       t=log L.
```

Therefore two sizes obey the exact matrix-pencil identity

```text
M(t2) M(t1)^-1 = [[1,0],[k(t2-t1),1]].
```

The unknown logarithmic mixing normalization occupies only the lower-left entry.
The other three entries give `T00=1`, `T01=0`, and `T11=1`.  But there is an
important exact degeneracy: if `LL=0` at both sizes, these three equations hold
for *any* sequence of `DD/LD`.  The two-size pencil is therefore only a matrix
re-expression of the per-size null-field condition, not independent evidence
for a constant Jordan flow.

New flow information first appears with three sizes.  For every adjacent pair
define

```text
k_i = T10_i / log(L_(i+1)/L_i).
```

A single dilation generator requires all `k_i` to agree.  Four sizes give
three rates and a two-degree-of-freedom constant-flow challenge.  This is the
actual projective Jordan fingerprint frozen by the scorer.

There is also an exact no-go that matters for the present raw blocks.  If the
two fields at size `L_i` are independently rescaled, `M_i -> D_i M_i D_i`,
then `J_i` survives but the cross-size pencil does not.  Dividing by `LD`
removes the overall scale, not the size-dependent relative gauge.  The scorer
therefore requires one declared ratio `A_L/A_D` per size for physical fields
`A_L L_raw` and `A_D D_raw`.  It propagates the full three-by-three covariance
and reports the source Gram condition number, but refuses to invent this
normalization.  Phase A does not yet contain `pi_a`, so its raw output alone
cannot legitimately score the constant-flow relation.  This identifies the
exact extra scalar that the next measurement must supply and prevents the
automatic two-size unipotent identity from being mistaken for a discovery.

## Remaining normalization boundary

The raw three-correlation block is sufficient for the first finite-size
study and for any later choice of `kappa`.  It is not sufficient by itself to
form the paper-normalized `eta_a^delta`: the periodic one-arm probability
`pi_a` needs a separately frozen macroscopic-radius convention.  Nor can the
three raw covariances uniquely recover `C1 CL/C2`, because adding a multiple
of `phi` changes the representative of the logarithmic partner without
changing the module.  No universal coefficient should be reported until
those normalizations are fixed.
