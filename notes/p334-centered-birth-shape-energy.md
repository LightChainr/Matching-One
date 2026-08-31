# A positive two-birth shape measure behind the centered thermal moment

The new W-squared direction should be placed in the shape of the **complete**
topology transition, rather than compared only with an uncentered J1 whose
linear birth-center term is much noisier. This note supplies an exact,
non-fitted change of observable and its binomial-smoothing correction.

## 1. A probability measure, not an arbitrary moment subtraction

On a complete path let K1 and K2 be the first and second rank births,
`F_i(p)=Pr[Binomial(N,p)>=K_i]` and `A=F1+F2-1`. Both births are between
1 and N. The derivative

```
q(p) = A'(p)/2
     = [BetaDensity(p;K1,N+1-K1)+BetaDensity(p;K2,N+1-K2)]/2
```

is a normalized positive probability density. The derivative identity follows
by differentiating the binomial tail: its telescoping derivative is
`N! / [(K-1)!(N-K)!] * p^(K-1)*(1-p)^(N-K)`. Therefore

```
E_q[p]   = C/(N+1)
E_q[p^2] = [C^2+C+W^2/4]/[(N+1)(N+2)],
```

where `C=(K1+K2)/2`, `W=K2-K1`. Integration by parts gives
`J0=1-2 E_q[p]` and `J1=1/2-E_q[p^2]`. These remain true after averaging
over the original random paths.

At the fixed, inherited reference p0=p_ref, define the canonical shape energy

```
Q_ref = E_q[(p-p0)^2]
      = 1/2-J1-p0*(1-J0)+p0^2.
```

Consequently its paired orientation contrast is exactly
`Delta Q_ref = -Delta(J1-p0*J0)`. The coefficient p0 is fixed by the existing
research reference, not chosen to minimize the observed error. A reduction
in error relative to uncentered J1 concerns a changed, explicitly named
observable; it does not improve the original J1 estimate by fiat.

## 2. Remove the known Beta smoothing exactly

Let Y choose `K1/(N+1)` or `K2/(N+1)` with equal probability, in addition to
the original path randomness. Its fixed-reference squared displacement is

```
R_ref = E[(Y-p0)^2]
      = E[(C-(N+1)*p0)^2+W^2/4]/(N+1)^2.
```

This separates two nonnegative per-path components: displacement of the
birth center from the reference and squared lifetime. Set
`a_N=((N+1)*p0+1/2)/(N+2)`. Without any large-N approximation,

```
Delta R_ref = -(N+2)/(N+1) * Delta(J1-a_N*J0).
```

The finite-N coefficient is slightly different from p0 because J1 includes
the order-statistic smoothing term C. It is still fixed before this readout
and has not been optimized on the orientation contrast.

The canonical and intrinsic variances of these two positive measures obey

```
Var_q(p) = (N+1)/(N+2)*Var(Y) + m*(1-m)/(N+2),
m = E[Y] = E[C]/(N+1),
Var(Y) = [Var(C)+E[W^2]/4]/(N+1)^2.
```

To derive the first line, condition on Y. The corresponding Beta distribution
has mean Y and variance `Y*(1-Y)/(N+2)`; apply total variance. The second
line also shows why lifetime variance alone is insufficient: both the mean
lifetime squared and its variance contribute to E[W^2].

The pooled variance Var(C) is over original paths, not variance of batch
means. Its uncertainty must be propagated by deleting an original batch and
recomputing the pooled centering. It is not a new independent source. Likewise
the difference between R_ref and Var(Y) is the squared displacement of the
ensemble mean from p0; it must not be called broadening automatically.

## 3. Concrete readout and scientific use

`scripts/p334_centered_birth_shape.py` reads the already completed K1/K2 CSVs
at `9c495ab13e65f2bc93dc0849ee3b73f88724c4b1`, without replaying paths,
solving clocks or drawing samples. It supplies original 20-batch vectors for
both orientations with four columns each:

- intrinsic fixed-reference center energy;
- intrinsic lifetime-square energy;
- canonical fixed-reference shape energy;
- J0, to retain the shared birth-center direction.

The existing covariance coordinator will join these columns with the pooled
mean/variance decomposition. This asks whether the N425 lifetime coordinate
survives in complete centered shape, or instead compensates a center-shape
change. Neither outcome requires selecting a higher moment order or fitting
another correction law. H4 normalization and geometry conventions remain
those of the original archive; this is a new analysis of that source, not
independent evidence or a continuum-field identification.
