# Same high-pass energy without six-level subtraction: an exact bridge

The 20k pilot (`1751e1c`) established an acquisition problem, not absent
high-degree topology. Here is a same-estimand escape route, together with the
specific sampling obstacle that must be overcome before running it.

## 1. Exact positive spectral representation

Let `K(t)=<F,T_t F>=sum_j w_j t^j`, where `w_j=||F_j||^2>=0`, and
`t_l=2^-l`, l=0,...,5. Put `a=h_5=9765/32768`.
The coefficients of the divided difference satisfy

`c_l = a / product_{m != l}(t_l-t_m)`.

Thus the *same* energy measured in the pilot is

`A=sum_l c_l K(t_l)=a [t_0,...,t_5]K`.

The Hermite--Genocchi identity, with `W~Dirichlet(1,1,1,1,1,1)` and
`T=sum_l W_l t_l`, gives

`A = a/5! E_W[K^(5)(T)]`.

The factor 1/5! is necessary: a normalized uniform Dirichlet distribution is
5! times simplex Lebesgue measure. The machine certificate verifies all
original coefficients and multipliers through degree 20 exactly. Equivalently,
`h_j=a h_complete_(j-5)(t_0,...,t_5)` for j>=5.
The classical integral normalization is also given in equation (52) of
[de Boor, Divided Differences](https://ftp.cs.wisc.edu/debooron/deboor2.pdf).

## 2. Fifth mixed differences

Use normalized Rademacher derivatives at p=.5:

`D_S F = 2^-5 sum_{u in {0,1}^S} (-1)^(5-|u|) F(u,X_outside)`.

For any five distinct bonds S, this derivative annihilates degree<=4
**pointwise**, before averaging. Fourier orthogonality gives

`K^(5)(t)=5! sum_|S|=5 <D_S F,T_t D_S F>`.

Consequently

`A=a E_T sum_|S|=5 E[conj(D_S F(X)) D_S F(Y_t)]`.

This removes the six coefficients up to 1984. The supplied kernel computes
each derivative using exactly 32 vertex evaluations. A noisy pair costs 64
F evaluations, or 192 N112 topology classifications. No new simulation is
performed in this bridge commit.

## 3. Positive population is not a positive single-pair sample

With a common ancestor Z and two conditionally independent `sqrt(t)` noise
arms X,Y,

`<D_S F,T_t D_S F> = E_Z |T_sqrt(t) D_S F(Z)|^2 >=0`.

The **conditional expectation squared** is positive. An individual unbiased
pair product need not be: for a degree-six Walsh mode, take five derivative
coordinates; the derivative is the remaining spin, and opposite outside
spins give product -1. The certificate includes this exact counterexample.

Squaring one noisy derivative instead estimates its unfiltered energy and
changes the estimand. Squaring a finite noisy conditional mean adds inner
variance; its unbiased correction can again be negative. A literally
nonnegative same-estimand implementation needs exact conditional integration
or another representation, not a renamed pair product.

## 4. Why uniform five-bond production is the wrong next run

There are `C(224,5)=4,493,032,544` subsets. If S is sampled uniformly, the
importance multiplier is `a*C=1,338,942,345.952...`.

For a single normalized degree-five Walsh mode supported on S0, all
derivatives except S=S0 are exactly zero. The estimator is therefore

- `a*C` with probability `1/C`;
- zero otherwise.

Its mean is the correct `a`, but its variance is exactly
`a^2(C-1)=399,010,376.139...`. This is worse than the old pilot's approximate
degree-five per-replica variance `1.09e6`. It is an explicit sparse-support
obstruction, not a universal lower bound for topology or adaptive sampling.
If the support is known and sampled with probability one, the same control
has constant response a and zero variance.

## 5. The narrow next gate

Reuse a topology-aware active/pivotal candidate set, but record exact proposal
probabilities. For a declared proposal `q(S|X,Y,T)`, the paired integrand is

`a conj(D_S F(X)) D_S F(Y) / q(S|X,Y,T)`.

Full support is required unless omitted derivatives are proved zero. An
epsilon-uniform mixture is a valid support safeguard, not a guarantee of good
variance. Selecting only currently pivotal bonds is **not** automatically
safe: a fifth mixed derivative can be nonzero when no individual coordinate
is pivotal at the starting configuration.

The next finite experiment should compare variance using those actual
inclusion weights. The minimal kernel already accepts a declared probability;
it does not invent or fit a proposal. No new raw block, physical field name,
state count, or N112 continuum prediction is introduced here.
