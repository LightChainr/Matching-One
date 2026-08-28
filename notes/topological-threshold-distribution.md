# The matching function as a topological-threshold distribution

Status: exact finite-size probability interpretation plus scaling conjectures. Literature overlap is marked explicitly.

## 1. Exact CDF from monotonicity

Mertens and Ziff show for every finite torus that

\[
M_L(p)=R_L^x(p)-\widehat R_L^x(1-p),
\qquad x\in\{c,b,e,h\},
\]

is monotonically increasing and tends from `-1` to `+1` across `[0,1]`. Therefore

\[
\boxed{F_L(p)=\frac{1+M_L(p)}2}
\]

is exactly a cumulative distribution function, with density

\[
\boxed{\rho_L(p)=\frac12M'_L(p).}
\]

The ordinary matching root

\[
M_L(p_L^*)=0
\]

is exactly the **median** of this finite-size topological-threshold distribution.

The root/median language is a reinterpretation of their matching function, not a new finite-size identity.

## 2. Configuration-level construction from cross wrapping

Assign every site an independent `U_v ~ Uniform(0,1)`. At parameter `p`, make a site black when `U_v <= p` and white otherwise. Connect black sites on the primary lattice and white sites on the matching lattice.

Define

\[
D_c(p)=I\{\text{black cross-wraps}\}
-I\{\text{white matching configuration cross-wraps}\}.
\]

Matching topology gives

\[
D_c(p)\in\{-1,0,+1\}.
\]

As `p` increases, black cross-wrapping is increasing and white cross-wrapping is decreasing; the two cross-wrapping events cannot coexist. Hence every coupled realization has the monotone form

```text
-1  -->  0  -->  +1
```

with either transition allowed to coincide.

Define

- `T_-`: the parameter at which white cross-wrapping disappears;
- `T_+`: the parameter at which black cross-wrapping appears.

Then `T_- <= T_+` and configuration by configuration

\[
\frac{D_c(p)+1}{2}
=\frac12I\{T_-\le p\}
+\frac12I\{T_+\le p\}.
\]

Averaging gives the exact representation

\[
\boxed{
F_L(p)=\frac12P(T_-\le p)+\frac12P(T_+\le p).
}
\]

Thus `rho_L` is the equal mixture of two concrete topological transition-threshold distributions.

## 3. Newman-Ziff rank representation

Sort the iid uniform labels. Their ordering is a uniformly random permutation independent of the order-statistic values.

Let `K_-` and `K_+` be the occupation ranks of the two topology transitions. Conditioned on `K=k`, the actual threshold is an order statistic,

\[
T\mid K=k\sim\operatorname{Beta}(k,N+1-k)
\]

for the convention in which the transition occurs when the `k`-th site is occupied. Off-by-one conventions must be frozen with exact small-system tests.

Therefore

\[
\rho_L(p)=\frac12\sum_kP(K_-=k)\,\mathrm{BetaPDF}_{k,N+1-k}(p)
+\frac12\sum_kP(K_+=k)\,\mathrm{BetaPDF}_{k,N+1-k}(p).
\]

The two threshold-rank histograms are consequently a compact sufficient representation of the entire **cross-channel** canonical matching curve.

## 4. Algorithmic consequence

For each random site permutation, a coupled forward/reverse connectivity implementation can record

```text
K_minus, K_plus
```

plus desired control variates.

Afterward the beta mixture can reconstruct

- `M_L(p)` for arbitrary `p`;
- the root/median;
- derivatives analytically;
- quantiles and moments;
- a standardized scaling profile.

Keep full microcanonical data during reference validation until the threshold-rank equivalence is exhaustively checked.

## 5. Scaling-limit random variable

Near criticality,

\[
M_L(p)\to\mathcal M(z),
\qquad z=b(p-p_c)L^{3/4}.
\]

Choose `T_L` with equal probability from `T_-` and `T_+` and define

\[
Z_L=b(T_L-p_c)L^{3/4}.
\]

Then scaling of the matching function is equivalent to distributional convergence

\[
Z_L\Rightarrow Z,
\qquad
P(Z\le z)=\frac{1+\mathcal M(z)}2.
\]

Since `mathcal M` is odd,

\[
\rho(z)=\frac12\mathcal M'(z)
\]

is even. The limiting topological-threshold distribution is therefore symmetric about zero for a fixed torus shape.

## 6. Existing and reinterpreted location estimators

### Median

`M_L=0` is the median. Its unusually rapid convergence is empirically near `L^-4` for square-site matching.

### Mode

Since

\[
\rho_L'(p)=\frac12M_L''(p),
\]

the central solution of `M_L''=0` is the mode of the threshold density. Mertens-Ziff explicitly studied this derivative estimator and found a much slower apparent convergence near `L^-1.67` on their sizes.

### Mean — already studied in 2016

The threshold-distribution mean is

\[
\mu_L=\int_0^1p\rho_L(p)\,dp
=\frac12\left(1-\int_0^1M_L(p)\,dp\right).
\]

This is **exactly Eq. (41) of Mertens and Ziff (2016)**. They numerically found this estimator to converge only around `L^-1.65`, so it is not a new fast threshold estimator.

The rank interpretation nevertheless gives the cheap representation

\[
\mu_L=\frac{E[K_-]+E[K_+]}{2(N+1)}
\]

under the fixed rank convention.

Furthermore, because finite `M_L(p)` has integer power-basis coefficients, the exact finite-size mean is rational. This arithmetic observation may be useful for exact-enumeration diagnostics, but the known slow convergence makes it a low-priority route for improving `p_c`.

**Scientific lesson:** the spectacular cancellation in the matching median/root is a special central finite-size structure. Do not assume that global statistics of `rho_L` inherit `L^-4` or `L^-7` convergence.

## 7. Distribution-asymmetry diagnostics

Although the mean itself is slow, differences among location statistics may still isolate correction sectors:

\[
\delta_{\rm mm}=\text{mean}-\text{median},
\qquad
\delta_{\rm md}=\text{median}-\text{mode}.
\]

Both vanish in the symmetric scaling limit. Their exponents probe finite-size skew/asymmetry rather than the width of the critical window.

Also preserve the standardized skewness of the finite threshold density. It may identify why median, mode and mean have very different finite-size biases.

## 8. Global universal shape invariants: useful, but expect slower corrections

The threshold distribution gives global scale-free quantities such as

- standardized kurtosis;
- standardized sixth moment;
- interquantile ratios;
- peak-height times standard deviation.

They remain useful as fingerprints of the universal limiting profile and may be statistically much more stable than fifth or seventh derivatives.

However, the slow convergence of the mean is a warning that global moments can couple strongly to generic correction fields. Use cross-model profile collapse and held-out size tests; do not expect the special root cancellation automatically.

## 9. Tail structure has an integrable-field-theory anchor

The correlation length obeys

\[
\xi\sim |p-p_c|^{-4/3}.
\]

Delfino and Viti (2011) used the integrable `Q -> 1` Potts field theory to derive, for off-critical rectangular crossing in the massive regime,

\[
P_{\rm cross}(L,R)
\sim A\frac{L}{\xi}e^{-R/\xi},
\qquad
A=\frac12(3-\sqrt3),
\]

in the appropriate large-distance limit. Earlier Newman-Ziff simulations also confirmed the `4/3` stretched-exponential character of spanning-probability tails.

For square-like geometry, `L/xi` scales as `|z|^(4/3)`. Therefore a testable qualitative prediction for the matching threshold density is

\[
\log\rho(z)
=-C|z|^{4/3}+O(\log|z|)
\]

with geometry/topology-dependent prefactors. The exact periodic cross-wrapping prefactor need not equal the open-rectangle amplitude above; only the massive scaling mechanism is being imported.

If

\[
1-\mathcal M(z)\sim C_0s e^{-c s},
\qquad s=z^{4/3},
\]

then differentiating suggests a leading density tail of the schematic form

\[
\rho(z)\propto z^{5/3}e^{-c z^{4/3}}
\]

up to topology-dependent powers and subleading terms. Fit the exponential exponent first; do not overinterpret the prefactor without a torus-specific derivation.

## 10. Joint threshold gap: genuinely additional information

The one-dimensional mixture `F_L` loses the pairing between `T_-` and `T_+`. Preserve their joint distribution whenever cheap.

Define

\[
G_L=T_+-T_-\ge0.
\]

The critical window suggests

\[
G_L=O(L^{-3/4}),
\]

so

\[
G_z=bL^{3/4}G_L
\]

may have a nondegenerate universal limiting distribution for fixed shape.

In Newman-Ziff rank units the corresponding gap should be of order

\[
K_+-K_-=O(NL^{-3/4})=O(L^{5/4}).
\]

This is a new observable not determined by `M_L(p)` alone. Its relation, if any, to pivotal-site / multi-arm structure should be investigated rather than assumed.

## 11. Practical server outputs

If cheap, preserve for every permutation or in aggregate:

- histogram of `K_minus`;
- histogram of `K_plus`;
- **joint** histogram or at least moments/covariance of `(K_minus,K_plus)`;
- histogram/moments of `K_plus-K_minus`;
- sample count and RNG counter range;
- exact period matrix and topology convention.

The cross-threshold histograms can later regenerate the whole canonical matching curve, while the joint gap contains information that the ordinary matching function discards.

## References / prior overlap

- Mertens and Ziff, Phys. Rev. E 94, 062152 (2016): exact matching function, monotonicity, scaling function, median/root, mode-like `M''=0` estimator, and the integral/mean estimator Eq. (41).
- Ziff and Newman, Phys. Rev. E 66, 016129 (2002): crossing-probability mean/median/mode estimators and their distinct finite-size corrections.
- Newman and Ziff, Phys. Rev. E 64, 016706 (2001): microcanonical algorithm and stretched-exponential spanning tails.
- Delfino and Viti, J. Phys. A / arXiv:1110.6355 (2011): integrable-field-theory off-critical crossing tails.
