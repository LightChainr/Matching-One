# A local colour gas carries the bulk law; a winding projection fixes phase weights

**New structural consequence.** The named closed source has, at
`m=exp(t)` integer, an exact representation as a local gas with m² active
colours and one vacant state, multiplied by the bounded topological factor
`m^(-r)`. This factor does not change the thermodynamic pressure at fixed
t, but it cancels the m²-fold full-phase degeneracy exactly. Deleting it
shifts the finite-volume matching root at strong coupling by
`2t/N+O(exp(-2t))` in logit coordinates.

This is a derivation from the fixed action, not a newly fitted source or
a claim that this constrained gas has the standard Potts critical line.
No new samples, enumeration or numerical grid enter the result.

## 1. Exact finite local weights

Keep the honest square torus and the same occupied NN graph. Write N for
its number of sites, K for occupation, B for occupied NN edges, C_B for
occupied components, and r for the ambient image rank. The established
source identity is

```text
S*=2 beta1-r-3K+2N+1
  =2B+2C_B-5K-r+2N+1.
```

Use y=p/(1-p), m=exp(t), Qc=m², and the site activity a=y/m^5. For integer
m>=2 assign each vertex a state in `{0,1,...,Qc}`. State0 means vacant;
the other states are active colours. The local edge weight is

```text
W(c_u,c_v)=1      if at least one endpoint is vacant,
          =Qc    if both are active with equal colour,
          =0     if both are active with different colours.
```

Let Z_col be the sum of `a^K product_e W(c_u,c_v)` over all colours.
For a specified occupied subset A, every occupied connected component
has one common colour. Summing its colourings gives exactly

```text
Z_col(y,m)=sum_A (y/m^5)^K m^(2B+2C_B).
```

The original occupation-odds partition function is therefore

```text
Z_star(y,m)=sum_A y^K m^(S*(A))
          =m^(2N+1) sum_colours a^K product_e W(c_u,c_v) m^(-r(A)).   (1)
```

For the normalized Bernoulli convention, multiply the full partition
function by `(1-p)^N`. It is configuration independent and does not alter
any normalized comparison here. The colouring construction is finite and
nonnegative; colours on adjacent occupied vertices are forced to agree.
It is not the unconstrained finite-temperature Potts edge interaction.

The global factor is the same zero-winding projection identified by the
[two-current representation](closed-source-two-current-representation.md):
one of two currents loses r free ambient coordinates. No new coupling
parameter, colour-fitting procedure or altered source convention is used.

## 2. The winding projection is invisible to bulk pressure, not to q/E

Since r is always0,1 or2, equation(1) gives the exact finite bounds

```text
m^(-2) Z_col <= m^(-(2N+1)) Z_star <= Z_col,
0 <= (log Z_col-log[m^(-(2N+1))Z_star])/N <= 2 log(m)/N.             (2)
```

Thus at any fixed finite t the pressure densities have the same limit,
if that limit exists along the chosen sequence, after the explicit
configuration-independent constant is removed. This does not require
inferring a phase transition or taking an unproved thermodynamic limit.
Any subsequential limits also agree by the bound.

The normalized rank-sector law, however, obeys exactly

```text
P_star(r=j)=m^(-j) P_col(r=j) / sum_(k=0)^2 m^(-k) P_col(r=k).      (3)
```

The pressure bound therefore does not justify dropping the projection
in q=r-1, E=q² or their root/slope-normalized U. All information entering
those observables is precisely in the three weights changed by(3).
This distinction provides a concrete bulk-versus-topological mechanism,
without a new empirical descriptor.

## 3. It exactly removes the ordered colour multiplicity

In the two-state coordinate h=y/m, the local colour activity is
`a=h/m^4=h/Qc²`. Each full single-colour state has weight

```text
a^N Qc^(2N)=h^N.
```

There are Qc such states. Before projection their total weight is
`Qc h^N`, compared with weight1 for the empty state. Full occupation has
r=2, so the factor m^(-r)=1/Qc cancels that colour multiplicity exactly.
The projected full weight is h^N, also the result of the previously
derived defect action `h^K m^(-g)`.

This explains why the projected strong-coupling matching root balances
empty and full states at h->1. It is not an accidental equality of two
energy minima with an overlooked entropy factor. The winding correction
fixes their relative finite-volume weight.

## 4. A precise root shift if the projection is deleted

Define the comparison law by deleting only the factor m^(-r) in(1).
Equivalently its occupation source is
`S_drop=2 beta1-3K+2N+1=S*+r`. It is a specified counterfactual law,
not a proposed replacement source. Both finite laws have unique simple
homogeneous matching roots for t>=0: the argument using supermodularity
and a pivotal configuration applies to 2 beta1 as well as to J.

At a fixed p and t>0, the unprojected law is the projected law tilted by
exp(tr). Its expected r strictly increases, since the finite full-support
law has more than one rank value. This holds separately in each geometry.
Consequently their pooled matching roots satisfy the exact order

```text
p_drop(t) < p_star(t),       t>0.                                  (4)
```

The [four-edge cut bound](closed-source-two-state-turnover.md) makes the
strong-coupling difference explicit. For the unprojected law put

```text
d=h m^(2/N),
eta(A)=g(A)-r(A)+2K/N.
```

Its weight, after removing a common factor, is exactly
`d^K m^(-eta(A))`. Empty and full states have eta=0. Any proper nonempty
occupied set has

```text
eta=Bmix-2C_B+2K/N >= 2C_B+2K/N >= 2+2/N.
```

Hence the unique pooled matching root obeys

```text
d0=1+O(m^(-(2+2/N))),
logit p_drop(t)=(1-2/N)t+O(exp[-(2+2/N)t]).                         (5)
```

For a formal analytic justification take xi=m^(-1/N). All exponents
N*eta are nonnegative integers, and the limiting pooled matching mean
is `(d^N-1)/(1+d^N)` with nonzero derivative N/2 at d=1. The analytic
implicit-function theorem gives(5), without guessing a saddle point.

The fixed original source already has
`logit p_star(t)=t+O(exp(-2t))`. Subtraction yields

```text
logit p_star(t)-logit p_drop(t)=2t/N+O(exp(-2t)).                   (6)
```

For a fixed geometry dilated in both directions by k, N becomes k²N;
the explicitly determined coefficient2/N therefore decreases by k².
For N25 the two limiting logit slopes are1 and23/25; at N100 they are1
and49/50. These are algebraic size predictions for a named comparison,
not newly measured thresholds or a fitted scaling exponent.

## Scientific consequence and limit order

The rank projection changes the pressure by at most2t/N at fixed t,
yet exactly removes the full-phase colour degeneracy and imposes the
finite-root drift(6). A bulk-equivalent local model need not have the
same finite topological matching observer. This gives a precise role
for the rank correction in the already closed source.

All strong-coupling assertions take fixed N before t->infinity. The
remainders are not uniform in N; equations(2) and(6) must not be combined
to infer a thermodynamic critical line, first-order transition, continuum
field or a double-scaling crossover. The next size analysis must retain
this limit order and the winding-sector costs, not promote local colour
degeneracy alone to a universality-class claim.
