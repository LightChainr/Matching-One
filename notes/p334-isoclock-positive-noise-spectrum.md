# P334 × P437: identical blockade means, strictly different positive noise response

For the specified five-site double star G and C4-plus-isolate H, at fixed
p=1/2 the variance over a shared blockade mask satisfies the exact identity

\[
\boxed{\operatorname{Var}_B m_G(B)-\operatorname{Var}_B m_H(B)
       =\frac{\rho(1-\rho)^2}{64}>0\quad(0<\rho<1).}
\]

Here `m(B)` is the **exact conditional safe probability** under that mask.
The two graphs have identical complete unmarked clocks, identical mean safe
responses for every p, and identical mean responses under every uniform
independent-blockade strength. Yet their same-mask conditional-response
variances separate them at every nontrivial correlation. This connects the
P334 uniform-blockade semigroup to positive Fourier degree energies, not to a
new production spectrum.

## Shared-mask identity

Let f be a finite safe indicator on n sites, B_i independent Bernoulli(a),
and U_i,V_i independent Bernoulli(u), independent of B. Set X_i=B_i U_i and
Y_i=B_i V_i. Conditional on B, X and Y are independent with the same law, so

\[
\operatorname{Var}_B E[f(BU)\mid B]
=\operatorname{Cov}(f(X),f(Y)).
\]

Each replica has marginal occupation p=au, while the standardized coordinate
correlation is

\[
\rho=\frac{\operatorname{Cov}(X_i,Y_i)}{p(1-p)}
=\frac{(1-a)u}{1-au}.
\]

Different coordinates remain independent. In the orthonormal p-biased
product Fourier basis, cross-products from unequal supports vanish and an
equal support S contributes rho^|S|. Hence

\[
\boxed{\operatorname{Var}_B m(B)
=\sum_{k=1}^n\rho^k E_k(p),\qquad
E_k(p)=\sum_{|S|=k}\widehat f_p(S)^2\ge0.}
\]

At fixed p, this family is realized by
`u=p+rho(1-p)` and `a=p/[p+rho(1-p)]`. At rho=0 the mask is deterministic;
at rho=1 the two replicas equal B. Thus its endpoints are0 and Var_p(f).
Positivity concerns this exact population variance and the energies, not
every individual noisy centered-replica product.

This refines the mean closure in
[the uniform-blockade semigroup note at d53db2f3](https://github.com/LightChainr/Matching-One/blob/d53db2f3/notes/p334-uniform-blockade-clock-semigroup.md):
the mean depends only on the old cardinality clock, but retaining the same
spatial mask in two replicas accesses information beyond it. The parameter p
here is occupation among these five selectable sites, not a full-N canonical
occupation parameter after a nonzero-prefix checkpoint.

## Exact32-point spectra

The source graphs and truth tables are unchanged from250c5899:
G edges01,02,03,14; H edges01,12,23,30 with isolated site4. Both have

\[
I(z)=1+5z+6z^2+2z^3,
\quad E_{1/2}f=\frac7{16},\quad
\operatorname{Var}_{1/2}f=\frac{63}{256}.
\]

Applying an integer Walsh transform to each saved32-entry indicator table,
with basis `chi_S(x)=(-1)^sum_(i∈S)x_i`, gives:

| Noise degree k | E_k(G) | E_k(H) | G−H |
| --- | ---: | ---: | ---: |
| 1 | 40/256 | 36/256 | 4/256 |
| 2 | 14/256 | 22/256 | −8/256 |
| 3 | 8/256 | 4/256 | 4/256 |
| 4 | 1/256 | 1/256 | 0 |
| 5 | 0 | 0 | 0 |
| Nonconstant total | 63/256 | 63/256 | 0 |

The constant energy E0 is49/256 for both. Consequently their full variance
difference is `(4rho−8rho²+4rho³)/256=rho(1-rho)²/64`; the strict inequality
above follows without a rho scan. At rho=1/2, the fixed-p realization is
a=2/3,u=3/4, and

\[
\operatorname{Var}_B m_G=\frac{393}{4096},\qquad
\operatorname{Var}_B m_H=\frac{385}{4096},\qquad
\Delta=\frac1{512}.
\]

## What geometric information entered

For a decreasing Boolean safe indicator, let I_v(p) be the probability that
site v is pivotal, with the other sites drawn independently at p. At p=1/2,
the positive singleton Walsh coefficient is I_v/2. The exact pivotal vectors
in the source vertex order are

\[
I_G=\frac18(5,3,1,1,2),\qquad
I_H=\frac18(3,3,3,3,0).
\]

Their sums coincide, so the common mean-response slope is−3/2. Their squared
norms differ, giving `E1=(1/4)Σ_v I_v²=40/256` versus36/256. The clock fixes
the uniform mean slope but not its distribution across sites. Higher-degree
energies carry the corresponding joint-support structure.

These fixed-p pivotal probabilities are not the previously computed final
birth-site probabilities π_v. The latter integrate pivotal probabilities
over the uniform insertion label, `π_v=∫_0^1 I_v(p)dp`. Squaring at fixed p
and squaring that integral are different observables. Both nevertheless
retain spatial information erased by the unmarked mean clock.

## Exact cross-route connection, with separate sources

The common structure with
[P437's positive difference bridge at79988f8d](https://github.com/LightChainr/Matching-One/blob/79988f8d/notes/p437-positive-difference-bridge.md)
is a positive noise-degree expansion `Σ rho^k E_k`, and the need to preserve
same-source replica/conditional-mean semantics rather than square a noisy
estimate indiscriminately. P437's
[coherent decomposition at888af29d](https://github.com/LightChainr/Matching-One/blob/888af29d/scripts/score_p437_coherent_decomposition.py)
separates a fixed-support coefficient from positive higher-support energy.

No P437 production data were recomputed or combined here. That N112 C3
observer, its fixed five-bond support, and its production dependency group
are distinct from this constructed five-site pair-trigger indicator. In
particular, **E5=0 in both present examples**: this is not evidence for the
P437 fifth-or-higher-degree signal. The bridge transfers an exact identity
and an identifiability mechanism, not measured energies or field labels.

## Artifact and lifecycle

`scripts/p334_isoclock_positive_noise_spectrum.py` reads the already-saved
truth tables and applies a32-point integer Walsh transform per graph. It
saves all Walsh numerators, exact E0..E5, pivotal vectors and the specified
variance in `results/p334-isoclock-positive-noise-spectrum/exact_positive_spectra.json`.
No new truth-table family, N425 sample, MonteCarlo run, network solve or server
connection is introduced. This is a theoretical finite counterexample based
on the P334 factor, not a fresh empirical evidence block.
