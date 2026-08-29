# Periodic-torus sufficient statistics for the Camia--Feng lattice log pair

Status: executable Phase-A0 protocol and exact algebra oracle for Issue #234.
No universal logarithmic coefficient is claimed or fitted here.

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

## Remaining normalization boundary

The raw three-correlation block is sufficient for the first finite-size
study and for any later choice of `kappa`.  It is not sufficient by itself to
form the paper-normalized `eta_a^delta`: the periodic one-arm probability
`pi_a` needs a separately frozen macroscopic-radius convention.  Nor can the
three raw covariances uniquely recover `C1 CL/C2`, because adding a multiple
of `phi` changes the representative of the logarithmic partner without
changing the module.  No universal coefficient should be reported until
those normalizations are fixed.
