# Exact conditional-replica martingale tomography

Issue #256 asks where a global topology observable becomes predictable, rather
than guessing another local field.  The exact gate uses the N=5 and N=10 truth
tables from PR #245 at \(p=2/5\).

For a centered two-vector \(Y\) and nested block-count sigma-fields,

\[
m_j=E[Y\mid\mathcal F_j],\qquad D_j=m_j-m_{j-1},
\]

the oracle exhaustively constructs every conditional-mean table and increment
table.  It then verifies, as exact rational identities,

\[
C_j=E[m_jm_j^T]
   =E[Y(X)Y(X')^T],
\]

where \(X,X'\) are independent conditional replicas, and

\[
\Gamma_j=E[D_jD_j^T]=C_j-C_{j-1}\succeq0,qquad
\operatorname{Cov}(Y)=\sum_j\Gamma_j.
\]

All cross-level matrices \(E[D_jD_k^T]\) vanish exactly.  A second computation
sums the four-replica refinement identity over two independent fine states
inside every coarse state and reproduces every \(\Gamma_j\).  Thus the PSD and
telescoping results are not consequences of subtracting two stored decimal
covariances.

## Krawtchouk bridge to PR #245

The first sigma-field is always \(\mathcal F_0=\sigma(K)\), total occupation.
The exact polynomial basis is

\[
K_r(K)=\sum_j {K\choose j}{N-K\choose r-j}
(1-p)^j(-p)^{r-j}.
\]

Its norm is \({N\choose r}[p(1-p)]^r\).  For every coordinate and degree, the
coefficient obtained from the conditional projection equals the average of the
PR #245 p-biased Boolean subset moments at that degree.  The reconstructed
\(m_0(K)\) and its covariance equal the direct block-count result exactly.
This identifies the radial/Krawtchouk slice without claiming that it contains
the spatial filtration.

## Exact structural result

The topology vector is

```text
[primal cross, matching-complement cross].
```

On N=5 it is completely radial: \(\mathcal F_0\) predicts 100% of its
variance and later increments vanish.  N=10 changes the picture sharply.  Its
filtration is

```text
total occupation
-> five norm-2 fiber counts {j,j+5}
-> individual sites.
```

Every topology increment is positive definite and rank two.  The off-diagonal
covariance signs and leading axes are:

| increment | sign | leading axis | trace fraction |
|---|---:|---:|---:|
| radial \(F_0\) | negative | \(-55.19^\circ\) | \(39877/93786\approx42.5\%\) |
| norm-2 fiber counts | positive | \(+55.84^\circ\) | \(2755/51744\approx5.3\%\) |
| singleton residual | negative | \(-61.87^\circ\) | \(37269/71456\approx52.2\%\) |

The principal-axis rotations between successive nonzero increments are
\(68.97^\circ\) and \(62.30^\circ\).  This is an exact finite graph where a
scalar shell summary would miss a genuine two-coordinate rotation.

The symmetry-control vector

```text
[orientation difference, matching-odd cross]
```

behaves differently.  Its cross covariance is zero at every level.  The N=10
ranks are `1 -> 1 -> 2`: radial and fiber information enter only the
matching-odd scalar, while orientation appears at the singleton level.  The
final rank two is therefore a direct sum of distinct symmetry sectors, not a
rotating plane.  This is the necessary representation warning for larger-N
rank claims.

## New risky mechanism conjecture

The exact N=10 pattern motivates a stronger replacement for the failed
constant odd-shell law:

> The large-N odd-shell anomaly is a rotating rank-two primal/matching topology
> plane.  Radial, mesoscopic/cover and microscopic increments retain two
> predictive directions, but their covariance axis shears with scale.  The
> apparent nonconstant signed shell is a moving scalar projection of this
> plane.

The conjecture is falsified if larger tori make the topology increments rank
one in one fixed basis, or if their leading two-dimensional subspace and axis
have no cross-size stability.  Rotation should only be claimed within one
symmetry block; mixing orientation-H4 with a scalar matching-odd coordinate
cannot manufacture it.

Reproduce the complete exact tables with:

```bash
python3 scripts/exact_conditional_replica_filtration.py \
  --output results/conditional-replica-filtration/exact_oracle.json
```

