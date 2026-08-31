# P398: finite projection memory explains the protected-ray crossing

Parent: `552c45d7595ebcb0d04555cec03b2a5bfd8da44a`.
This uses its width-eight positive continuous generator and the same fixed
AP/landing rays. No new samples, fitted kernel, width extension, or repeated
93+93 rank campaign are involved.

**Result:** the two rays have nearly the same memory duration, about .27,
but the initially faster ray has 6.81 times the initial feedback strength.
This reverses their relative decay rates and explains the normalized-ray
crossing. One explicit hidden geometry statistic per ray already captures
most of the crossing location.

## Exact projection, not an assumed non-Markov model

Let M=-G, use the physical stationary inner product
`<f,g>_pi=sum pi(state) conjugate(f(state)) g(state)`, and separately fix
`psi_s=(A+s exp(-i*pi/4)L)/sqrt(2)`, with s=-1 or +1.
Project orthogonally onto psi and its complement in the already established
ray sector. After stationary-metric whitening, the mass matrix is

\[
M_s=\begin{pmatrix}\omega_s&b_s\\c_s&D_s\end{pmatrix},
\qquad \dim D_s=92.
\]

For the normalized stationary correlation u_s(t), Schur elimination gives

\[
\boxed{u'_s(t)=-\omega_su_s(t)+\int_0^t k_s(t-r)u_s(r)\,dr,\quad
k_s(t)=b_se^{-D_st}c_s,\quad u_s(0)=1.}
\]

Equivalently,

\[
\widehat u_s(z)=\frac1{z+\omega_s-\widehat k_s(z)},\qquad
\widehat k_s(z)=b_s(zI+D_s)^{-1}c_s.
\]

This is exact finite-dimensional elimination. The orthogonal initial-force
term has zero correlation with the initial projected source; individual
stochastic trajectories are not asserted to obey a noiseless scalar equation.
The full configuration remains Markov. Remembering omitted configuration
coordinates is not evidence for morphism/path-order memory.

## The hidden force is an explicit T2/R geometry observable

The parent gives `GA=-3A+R` and `GL=-3L+T2`, where T2 is the character-i
size-two-cluster statistic and R is extra adjacent boundary-contact
multiplicity between the two endpoint clusters. Consequently

\[
\eta_s=(M-\omega_s)\psi_s
=(3-\omega_s)\psi_s-\frac{R+s e^{-i\pi/4}T_2}{\sqrt2},
\qquad \langle\psi_s,\eta_s\rangle_\pi=0.
\]

Thus the first omitted force is a named current-configuration statistic,
with no adjustable source coefficient. The curvature is

\[
u''_s(0)=\omega_s^2+k_s(0),\qquad
(\log u_s)''(0)=k_s(0).
\]

The process is nonreversible: the kernel pairs the right force eta with the
left force `(M^dagger-omega)psi`, not eta with itself. In particular,
`k(0)=b c` need not equal `||c||^2=Var_pi(eta)/Var_pi(psi)`.

| Quantity | psi-minus | psi-plus |
|---|---:|---:|
| Bare projected decay omega | 3.368820241 | 3.691415268 |
| k(0), logarithmic curvature | .4514576654 | 3.075443097 |
| Normalized right-force variance | .1286585420 | 2.939987703 |
| k(0) / right-force variance | 3.508959905 | 1.046073456 |
| Integral of k(t) | .1212071174 | .7859741072 |
| Signed first moment / integral | .2721429267 | .2710544714 |
| Smallest real part of hidden decay pole | 3.560987073 | 3.368873464 |

The force-covariance replacement would seriously understate the minus-ray
feedback. The left/right force alignments are .83036 and .83167; the
asymmetry is mainly their unequal magnitudes, not orthogonal forces.

## Actual memory kernels and what their poles mean

| t | k-minus(t) | k-plus(t) |
|---:|---:|---:|
| 0 | .45145767 | 3.07544310 |
| .1 | .30927832 | 2.02904524 |
| .26565732 | .16663410 | 1.05403712 |
| .5 | .07026457 | .43973179 |
| 1 | .01140390 | .07573949 |
| 2 | .00031359 | .00254267 |

There is no kernel sign change on the declared 401-point grid from 0 to 8.
This is not a positivity theorem: the spectral residues are signed.
In fact the plus-ray hidden block has a conjugate pair of decay poles
**5.29412605 +/- .30999968 i**, although the unprojected charge-i masses
were numerically real. The total real kernel remains positive on that grid.

These are poles of the self-energy, equivalently zeros of the source
resolvent `u_hat`, not poles of the full physical propagator. They are not
new oscillating physical states, a complex scaling dimension, or a Jordan
collision. Projection can create this structure without changing the
underlying stochastic generator.

## Why the curves cross

At zero distance the plus ray loses faster by
`Delta omega=.3225950262`. But its logarithmic curvature is larger by
`Delta k0=2.6239854317`. Therefore

\[
\log\frac{u_+(t)}{u_-(t)}
=-.3225950262\,t+1.311992716\,t^2+O(t^3).
\]

This no-fit short-distance expansion predicts a nonzero crossing at
**.2458817204**, versus the complete **.2656573200**. The corresponding
instantaneous decay rates first cross at **.1292994909**; their integrated
difference returns to zero at the later correlation crossing. No hidden
eigenvalue collision is needed.

The minimal explicit extension retains one normalized eta per ray. The
stationary-L2 Galerkin mass matrices, fixed by the original generator, are

```text
minus: [[3.368820241, 1.258628955],
        [ .358690036,4.576458844]]
plus:  [[3.691415268, 1.793638589],
        [1.714639234,3.828950570]]
```

This four-observable description uses only psi-minus, eta-minus, psi-plus,
eta-plus. It exactly matches each ray's value, slope and curvature at zero
but is not closed: the next hidden forward residual norms are 1.79510 and
1.47484. Its no-fit crossing is **.2541213924**, about 4.34% below the full
answer. This is a retrospective mechanistic compression, not a fresh
validation or an exact new two-state model.

The main physical lesson is now specific: R/T2 feedback is stronger in the
plus ray, rather than having a much longer duration. Eliminating that
feedback made the two-observable kernel fail the semigroup law.

## Reproduction and boundary

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /Users/lc/python-envs/research-py311/bin/python scripts/p398_width8_projected_memory.py
```

The script additionally caps runtime BLAS threads at one. The full hidden
poles, residues, covariance/force quantities, kernel grid and one Volterra
consistency evaluation are saved. The latter agrees at about 3e-15.
These are deterministic float64 projections of an inherited exact finite
generator, not new Monte Carlo evidence, continuum field counting, or
site-Matching/Jordan identification.
