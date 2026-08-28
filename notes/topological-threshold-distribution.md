# The matching function as a topological-threshold distribution

Status: exact finite-size probability interpretation plus scaling conjectures.

This viewpoint is potentially more useful than treating `M_L(p)` as merely a root-finding function.

## 1. Exact CDF from monotonicity

Mertens and Ziff prove that for every finite torus the matching function can be written as

\[
M_L(p)=R_L^x(p)-\widehat R_L^x(1-p),
\qquad x\in\{c,b,e,h\},
\]

and is monotonically increasing. For the cross-wrapping channel,

\[
M_L(0)=-1,\qquad M_L(1)=+1.
\]

Therefore

\[
\boxed{F_L(p)=\frac{1+M_L(p)}2}
\]

is exactly a cumulative distribution function on `[0,1]`, with density

\[
\boxed{\rho_L(p)=\frac12 M_L'(p).}
\]

The ordinary matching root

\[
M_L(p_L^*)=0
\]

is exactly the **median** of this finite-size topological-threshold distribution.

This interpretation requires no asymptotic scaling assumption.

## 2. Configuration-level construction from cross wrapping

Assign every site an independent `U_v ~ Uniform(0,1)`. At parameter `p`, make the site black when `U_v <= p` and white otherwise. Connect black sites on the primary lattice and white sites on the matching lattice.

Define

\[
D_c(p)=I\{\text{black cross-wraps}\}
-I\{\text{white matching configuration cross-wraps}\}.
\]

The topology of matching configurations implies

\[
D_c(p)\in\{-1,0,+1\}.
\]

As `p` increases:

- the black cross-wrapping event is increasing;
- the white cross-wrapping event is decreasing;
- black and white cross-wrapping events cannot coexist.

Hence for every fixed uniform field, `D_c(p)` is nondecreasing and has the form

```text
-1  -->  0  -->  +1
```

with either jump allowed to coincide.

Define two random thresholds:

- `T_-`: the `p` at which white cross-wrapping disappears;
- `T_+`: the `p` at which black cross-wrapping appears.

They obey `T_- <= T_+`, and configuration by configuration

\[
\frac{D_c(p)+1}{2}
=\frac12 I\{T_-\le p\}
+\frac12 I\{T_+\le p\}.
\]

Taking expectations gives

\[
\boxed{
F_L(p)=\frac12 P(T_-\le p)+\frac12P(T_+\le p).
}
\]

Thus `rho_L` is not merely an abstract density: it is the equal mixture of the disappearance and appearance threshold distributions.

## 3. Newman-Ziff rank representation

Under the same iid-uniform construction, sort the site labels. The ordering is a uniformly random permutation and is independent of the order-statistic values.

Let

- `K_-` be the occupation rank at which white cross-wrapping disappears;
- `K_+` be the occupation rank at which black cross-wrapping appears.

Conditioned on `K=k`, the corresponding actual threshold is the `k`-th order statistic of `N` iid uniforms,

\[
T\mid K=k\sim \operatorname{Beta}(k,N+1-k)
\]

for the convention in which the transition occurs when the `k`-th site is occupied. Boundary/off-by-one conventions must be fixed by exact small-system tests.

Therefore the entire canonical density is a **finite beta mixture** determined only by the two histograms of threshold ranks:

\[
\rho_L(p)
=\frac12\sum_k P(K_-=k)\,\mathrm{BetaPDF}_{k,N+1-k}(p)
+\frac12\sum_k P(K_+=k)\,\mathrm{BetaPDF}_{k,N+1-k}(p).
\]

This means a Newman-Ziff simulation need not store a dense probability grid to preserve the full matching crossover. The threshold-rank histograms are a sufficient compressed representation for the cross-wrapping channel.

## 4. Immediate algorithmic consequence

For each random site permutation, run coupled forward/reverse connectivity and record only:

```text
K_minus, K_plus
```

plus any additional observables/covariates desired.

From these ranks one can later reconstruct:

- `M_L(p)` at arbitrary `p`;
- the matching root / median;
- `M'_L`, `M'''_L`, ... analytically from beta mixtures;
- moments and quantiles of the threshold distribution;
- the standardized universal profile;
- uncertainty by resampling permutations or histogram counts.

For the cross channel this may be substantially more compact than storing every microcanonical wrapping indicator for every `k`.

Keep the full microcanonical data during the reference stage until the threshold-rank equivalence has been exhaustively tested.

## 5. Scaling-limit random variable

Near criticality,

\[
M_L(p)\to\mathcal M(z),
\qquad z=b(p-p_c)L^{1/\nu},
\qquad \nu=4/3.
\]

Define the rescaled random threshold

\[
Z_L=b(T_L-p_c)L^{3/4},
\]

where `T_L` is chosen with probability `1/2` from `T_-` and `T_+`.

Then the matching-scaling hypothesis is equivalent to convergence in distribution

\[
\boxed{Z_L\Rightarrow Z}
\]

with universal CDF

\[
P(Z\le z)=\frac{1+\mathcal M(z)}2
\]

for fixed torus shape.

Because `mathcal M` is odd, the limiting density

\[
\rho(z)=\frac12\mathcal M'(z)
\]

is even. Hence the limiting distribution is symmetric around zero.

This gives a direct probabilistic meaning to universality of the matching scaling function.

## 6. Median, mode, mean, and finite-size corrections

Several pseudo-critical estimators become ordinary location statistics.

### Median

\[
M_L(p)=0
\]

defines the median.

### Mode

Because

\[
\rho_L'(p)=\frac12M_L''(p),
\]

a root of

\[
M_L''(p)=0
\]

at the central maximum is the mode of the threshold density.

Mertens-Ziff's separate `M''=0` estimator can therefore be interpreted as a mode estimator rather than an unrelated derivative trick.

### Mean

The mean threshold is

\[
\mu_L=\frac12E[T_-+T_+].
\]

In the rank representation,

\[
E[T\mid K=k]=\frac{k}{N+1},
\]

so, modulo the fixed rank convention,

\[
\boxed{
\mu_L=\frac{E[K_-]+E[K_+]}{2(N+1)}.
}
\]

This is extremely cheap to estimate.

Since the limiting distribution is symmetric, mean, median and mode must all converge to `p_c`, but their finite-size biases can belong to different correction sectors.

This creates a new estimator family without introducing arbitrary fitting functions.

## 7. New finite-size diagnostics from distribution asymmetry

Define

\[
\delta_{\rm mm}=\text{mean}-\text{median},
\qquad
\delta_{\rm md}=\text{median}-\text{mode}.
\]

Both vanish in the symmetric scaling limit.

Their size exponents isolate finite-size **skew/asymmetry corrections** rather than the symmetric width of the critical window. This may expose the operator responsible for the different convergence exponents of existing pseudo-critical estimators.

Also compute standardized skewness of the finite threshold distribution. Its asymptotic decay should be predictable from the leading odd correction to the even universal density.

## 8. Global universal invariants may be statistically better than high derivatives

Local derivative ratios such as

\[
\kappa_3=\frac{M'''(p_c)}{M'(p_c)^3}
\]

are interesting but increasingly high derivatives are noisy.

The threshold-distribution view supplies global scale-free shape invariants:

- standardized kurtosis `mu_4/mu_2^2`;
- standardized sixth moment;
- interquantile ratios, e.g. `(q_0.9-q_0.1)/(q_0.75-q_0.25)`;
- peak-height times standard deviation;
- entropy after unit-variance standardization.

For fixed torus modulus, these should be universal across microscopic realizations after the usual scaling limit.

A candidate analytic scaling function should match both local invariants (`kappa_3`, `kappa_5`) and these global invariants.

## 9. Tail prediction

Away from criticality the correlation length scales as

\[
\xi\sim |p-p_c|^{-\nu}.
\]

At fixed scaling variable `z`,

\[
L/\xi\propto |z|^{\nu}=|z|^{4/3}.
\]

Wrapping/crossing failure probabilities in a massive phase are expected to be exponentially small in `L/xi`. This suggests the testable asymptotic form

\[
\log\rho(z)\sim -C|z|^{4/3}+\text{subleading logs/powers}
\]

rather than a Gaussian `exp(-z^2)` or logistic `exp(-|z|)` tail.

This is a scaling argument, not an exact theorem for the torus matching density. It should be tested numerically before using it to constrain analytic families.

If confirmed, the `4/3` stretched-exponential tail would be a strong fingerprint linking the full threshold distribution directly to the correlation-length exponent.

## 10. Practical server outputs

Add, if cheap, the following to every permutation-based simulation:

```text
K_minus
K_plus
```

or aggregated histograms thereof.

For large runs preserve at least:

- histogram counts of `K_minus` and `K_plus` separately;
- their joint covariance / joint 2D histogram if affordable;
- sample count and RNG counter range;
- exact period matrix and topology convention.

The **joint** distribution contains additional information: the gap

\[
G_L=T_+-T_-
\]

measures the sample-level interval between destruction of complementary white cross-wrapping and creation of black cross-wrapping. Its scaling law may be universal and is not determined by `M_L` alone.

That gap is a new observable worth preserving from the outset.
