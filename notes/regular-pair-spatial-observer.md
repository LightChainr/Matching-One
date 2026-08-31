# The Q-activated spatial observer is a colour contraction, not a mark covariance

The [regular finite local interaction](regular-pair-interaction-result.md)
has a nonadditive two-mark Q derivative on one fixed four-path occupation.
The next question is its spatial continuation after averaging original
occupations. The [fixed production contract](../analysis/p337_regular_pair_spatial_contract.json)
selects exactly L32/r8 and L64/r16, with two new random blocks. It does
not recycle the old N25 mark archive or fit a continuum exponent.

## 1. The exact observable before any simulation

Give two nonadjacent vacant-site summands independent coefficients
lambda_x and lambda_y multiplying the same canonical Kreg. For a fixed
original occupation A, normalize its colour contraction by its unmarked
weight and write

```
F_A(Q)=1+lambda_x beta_x(Q,A)+lambda_y beta_y(Q,A)
         +lambda_x lambda_y beta_xy(Q,A).
```

The old regularity theorem gives `beta_x(1,A)=beta_y(1,A)=beta_xy(1,A)=0`.
Consequently, for the finite partition function after summing occupations,

```
C_Q(x,y) = partial_lambda_x partial_lambda_y log Z
         = <beta_xy>_Q - <beta_x>_Q <beta_y>_Q,
C(x,y)   = partial_logQ C_Q|Q1
         = <g_xy(A)>_(ordinary iid percolation),
g_xy(A)  = partial_Q beta_xy(Q,A)|Q1.                       (1)
```

The background Q derivative multiplies beta_xy(1)=0. Both disconnected
product derivatives also vanish. Thus (1) includes the normalization
and connected subtraction exactly; it is not an uncentered replacement
for a covariance. In particular it is neither `Cov(a_x,a_y)` nor
`Cov(t_x,t_y)`. An occupation with either marked vertex occupied gives
zero and remains part of the probability law.

The lambda parameters are ordinary per-vertex parameters. No factor of
N or N squared is included in C. Uniform site-average epsilon derivatives
would sum over pairs and introduce their own explicit normalization.

## 2. The finite exact spatial kernel

For two nonadjacent sites, partition their eight incident edge-nodes by
the exterior hypergraph components. A vacant neighbour leaves an isolated
edge-node; an occupied neighbour contributes its original NN component.
Let pi be this partition and b its number of blocks. Then

```
g_pi = partial_Q [Q^(-b) sum_colours Kreg_x Kreg_y D_pi]|Q1.
```

One direct exact formula uses coarsenings sigma of pi. A coarsening with
k distinct colour labels has multiplicity `(Q)_k`. For k>=2,

```
partial_Q (Q)_k at Q1 = (-1)^(k-2) (k-2)!.
```

The one-colour term vanishes identically, while for k>=2 the multiplicity
itself vanishes at Q1. The derivatives of the regular kernel entries
therefore disappear. If `k4_x=4 Kreg_x(sigma;Q1)` and similarly for y,

```
16 g_pi = sum_(sigma coarsens pi, |sigma|>=2)
          (-1)^(|sigma|-2)(|sigma|-2)! * k4_x * k4_y.        (2)
```

This finite integer kernel supplies the simulation directly. It does not
estimate a response by subtracting nearly equal Q>1 simulations or build
a new feature regression. Its four-local-colour entries are `k4=-8`
for four distinct colours; `-6/-3` for one opposite/adjacent equal pair;
`-2/-1` for opposite/adjacent two-pair patterns; zero for3+1/all4.

If pi shares no exterior component between the sites, the two contractions
factor, each of order Q-1. If there is only one shared component, fixing
its colour makes each partial contraction constant by colour symmetry;
gluing divides the product by Q, which is regular at Q1. In both cases
g_pi=0. A nonzero signal therefore needs at least two components spanning
the two marked neighbourhoods. Its exact signed weights, rather than a
count of independent single-site marks, determine the observer.

## 3. What the fixed new data can decide

The primary question is whether the occupation-averaged noncontact
coefficient C64 at distance16 is nonzero. This is a possible averaged
cancellation question even though an individual far-connected occupation
has a nonzero kernel. A99% Monte Carlo interval excluding zero rejects
the finite contact-only null. An interval containing zero ends the
specified calculation unresolved; no additional sampling follows it.

Each geometry gets200000 new iid configurations, grouped into200 batches.
Within a configuration the32 translation/direction pairs are averaged,
not treated as separate independent observations. The exact shared-component
decomposition preserves its full joint batch covariance. C32 and the
fixed ratio C64/C32 are secondary spatial summaries, with an unbounded
ratio reported honestly when its denominator is unresolved.

Neither a scalar two-point function nor C4 lattice symmetry identifies
a continuum spin-four field. The fixed-cut colour recoupling already
contains several blocks. No x=17/4 or x=21/4 fit, selected angle, additional
source, or posterior support class is added to this experiment. The
kernel and the physically realized sampling rule are committed before
either production block is revealed.
