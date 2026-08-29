# A minimal identifiable RG dynamics: the affine clock pencil

Status: exact two-state theorem plus a post-reveal full-state extrapolation for
Issue #180.  This replaces a free transfer-matrix fit by two fixed algebraic
alternatives.

## 1. The two-state theorem

Let `A=exp(-t_A G)` and `U=exp(-t_U G)` be two clocks of one two-dimensional
generator.  If `A` is nonscalar, every matrix commuting with `A` lies in

```text
span{I,A}.
```

Hence one common two-state generator implies

```text
U = alpha I + beta A.                                  (1)
```

This is stronger than matching two separately fitted spectra.  For any frozen
sources and readouts, define context Hankel blocks

```text
H_I[l,r]=c_l^T b_r,
H_A[l,r]=c_l^T A b_r,
H_U[l,r]=c_l^T U b_r.
```

Equation (1) gives the exact, similarity-free gate

```text
H_U = alpha H_I + beta H_A,                            (2)
rank[vec(H_I),vec(H_A),vec(H_U)] <= 2.                 (3)
```

Two noncollinear context entries determine `alpha,beta`; every other entry is
held out.  No PCA, SVD cutoff, eigenvector alignment, or matrix logarithm is
needed.  Equation (1) also implies the #255 mixed rectangle `[A,U]=0`.  In a
minimal observable/reachable rank-two realization the converse holds: a
commuting second clock must obey (1).

## 2. Dyadic clocks give two fixed extrapolators

For the squared clock `U=A^2`, Cayley--Hamilton gives

```text
A^2 = tr(A) A - det(A) I.                              (4)
```

The two live finite-size mechanisms therefore have fixed coefficients.

### Ordinary analytic q=2 mode

After retaining the fixed-point coordinate, the eigenvalues are `(1,1/2)`:

```text
A^2 = -(1/2) I + (3/2) A,
x_(4N) = -(1/2)x_N + (3/2)x_(2N).                      (5)
```

### Rank-2 Jordan mode

After leading scaling is removed, `A=I+K`, `K^2=0`:

```text
A^2 = -I + 2A,
x_(4N) = -x_N + 2x_(2N).                              (6)
```

The familiar scalar q2 and logarithmic second-difference nulls are therefore
not separate curve fits.  They are the two Cayley--Hamilton polynomials of a
minimal two-state dynamics.  If the realization is genuinely common, the
same coefficients must hold for every readout of the full state.

## 3. Direct N145/N290 vector prediction

The executable scorer reconstructs, inside each delete-one block,

```text
x_N=(I_S,I_Du,T_D,T_Su)
```

from the existing P50 histograms and propagates the independent N145/N290
covariances.  The observed states are

```text
x145=(-0.00848091,-0.01132667,0.44674282,1.44580065),
x290=(-0.00607902,-0.01139607,0.72844620,1.95448331).
```

Without fitting a transfer matrix, (5)--(6) freeze the N580 point forecasts

```text
q2:
  (-0.00487808,-0.01143077,0.86929789,2.20882464)

Jordan:
  (-0.00367714,-0.01146547,1.01014958,2.46316598).
```

These are deliberately risky full-vector predictions.  A model does not pass
by matching only `T_Su`; one common affine coefficient pair must transport all
four typed readouts within their full covariance.  The predictions are
post-reveal constructions from N145/N290 and become evidence only on a future
N580 block.

## 4. Exact marked-sector discriminator

The exact #249 endpoint witness has endpoint eigenvalues `(1,2)`, so its
direct double cover obeys `U=A^2=-2I+3A`.  Those endpoint coefficients predict
value `-5` on the charged coordinate where `A_charged=-1`.  The direct-cover
charged endpoint is `U_charged=0`, giving affine defect `+5`.

Thus the same experiment distinguishes two possibilities without increasing
the scalar bulk rank:

```text
unmarked contexts satisfy one common (alpha,beta)
marked/charged context rejects those coefficients
```

The second line is cover-morphism memory, exactly the extra sector anticipated
by #249, rather than evidence for a third scalar correction exponent.

## Minimal two-clock acquisition

Measure `H_I,H_A,H_U` on the #255 two-source by two-readout mixed block.  Use
two entries to solve (1), hold out the other two, and score (2) jointly with
the path-order commutator.  Then add one allowed charged/marked entry.  A
rank-two bulk clock predicts the same coefficients on every unmarked entry;
the risky categorical prediction is that the marked entry alone has a stable
defect.

## Boundary

Cayley--Hamilton, the affine pencil, the commutator implication, and the exact
charged counterexample are algebra.  The four-channel lattice state being one
observable/reachable two-dimensional realization is the hypothesis.  The
N580 numbers are post-reveal predictions, not a claim about data not yet run.
