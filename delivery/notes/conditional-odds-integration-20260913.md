# Anchoring rare-sector odds from a geometric onset and conditional occupancy

Date: 2026-09-13. Completed identity and fixed-size executable reference.
This is an application of thermodynamic/path integration, not a new general method.

## 1. Eliminate absolute rare probabilities from the derivative, not from the model

At homogeneous site fugacity t=p/(1-p), let

    Z_j(t)=sum_{r=j} t^K,       mu_j(t)=E_t[K|r=j].

Then

    d/d log t log(Z_2/Z_0)=mu_2-mu_0.                      (1)

Suppose the lowest powers are a_j*t^{k_j}. Their coefficients are geometric
configuration counts, not fitted amplitudes. Write R_j(t)=Z_j/(a_j*t^{k_j}),
so R_j(0)=1. Integrating (1) with its EXACT lower-end anchor gives

    log(Z_2(t)/Z_0(t))
      = log(a_2/a_0)+(k_2-k_0)log t
        + integral_0^t [(mu_2(u)-mu_0(u))-(k_2-k_0)]/u du. (2)

The integrand is regular at zero: it equals R_2'/R_2-R_0'/R_0. The Bernoulli
normalizer (1+t)^N cancels in odds only. The matching root is where (2) vanishes.

For the 4x4 NN torus, k_0=0,a_0=1,k_2=7,a_2=16. These reproduce the established
row/column-cross onset, and the script verifies them from the exact sector counts.
The integrated conditional-mean equation returns

    p_*=0.59067211233102829689590201143951286962111713272216...

independently of the direct subtraction Z_2-Z_0 root calculation. Four parameter
checks agree with direct log odds to the declared working precision. This is a
70-digit numerical integration control, not a rational quadrature certificate.

## 2. Rigorous finite lower-end anchor error without knowing all coefficients

At positive t0, if only N,k_2,a_2 are known, total configuration-count bounds imply

    0 <= R_0(t0)-1 <= (1+t0)^N-1,
    0 <= R_2(t0)-1 <= sum_{j=1}^{N-k_2} binom(N,k_2+j)t0^j/a_2.

Because log(1+x)<=x, the absolute correction to the simple geometric anchor is
bounded by the maximum of these two explicit quantities. For N=16,k_2=7,a_2=16,
at t0=10^-12 this bound is 8.04375000000715e-10. The bound is exact rational and
is stored separately from the actual correction computed by the oracle.

This permits a nonzero lower endpoint in a future conditional-sampling integration
without inventing its normalization constant. A small t0 is not a guarantee that
conditional sampling or quadrature at other t is efficient.

## 3. Sampling allocation, with the missing costs stated

For a predeclared quadrature g approximately equal to anchor+sum_l w_l(mu_2l-mu_0l),
independent conditional sample means give variance

    sum_l w_l^2 [V_2l/n_2l + V_0l/n_0l],

where V_jl=Var(K|rank j,t_l). If one independent conditional observation has cost
c_jl, minimizing variance under sum c_jl*n_jl=B gives

    n_jl proportional to |w_l| sqrt(V_jl/c_jl),
    variance_min = [sum_jl |w_l| sqrt(V_jl*c_jl)]^2 / B.

This is standard optimal allocation. It does not include anchor error, quadrature
bias, dependence between conditional samples, adaptive-selection bias, or the cost
of generating the conditional normalizers. Those must be controlled before a run
is commissioned. A rejection sampler can bury the original rare-event cost in c_jl.
The exact rank bridge of #733 already KNOWS its normalizer; using it here is an
oracle check, not a way of learning an unknown normalizer for free.

Spatial Hessians require conditional covariances of H=sum h_i n_i, not only K.
The invisible-spatial-marks construction proves that the entire (rank,K) archive
can leave these spatial responses undetermined. Do not request another unmarked K
histogram as a substitute for the required spatial observable.
