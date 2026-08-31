# A source-normal response hidden behind negative self-curvature

The raw own-source center curvature is negative in all four size/geometry
readouts. Its decomposition reveals two opposing terms: a larger negative
term from the first-score component of source curvature, and a **positive
response orthogonal to both original source scores**. The latter survives
at about3.5–5.8 original-batch SE and is therefore a new finite-data
mechanism coordinate beyond the first-order source span.

## Measured result on the original00 new64 population

All entries below are `×10^-8`, mean±one original20-batch SE. Every
covariance is paired; rows retain the original full20000-prefix denominator.

| N / own receiving geometry | raw H_oo C | first-score tangent part | source-normal part |
|---|---:|---:|---:|
| 325 / first | -6.222±1.208 | -10.338±0.830 | +4.116±0.949 |
| 325 / second | -5.985±1.148 | -9.218±0.982 | +3.233±0.669 |
| 425 / first | -8.493±1.014 | -11.793±0.506 | +3.300±0.954 |
| 425 / second | -6.520±0.836 | -10.497±0.639 | +3.977±0.682 |

For the own-source `A(p_ref)` curvature the signs reverse in the same
decomposition. Its source-normal components are respectively
`(-4.873±0.935,-3.814±0.811,-3.713±1.153,-4.456±0.718)×10^-7`.
Thus the center and canonical A readouts give a coherent response shape:
the normal perturbation delays the center and lowers A at the fixed
reference point. They are correlated functions of the same birth clocks,
not independent confirmations.

Mixed fs components remain without a consistent cross-size pattern, both
before and after the source projection. The dominant new information is
the opposing own-source components, not a detected mixed-source interaction.
No claim of exact additivity is made.

## Exact source-space decomposition

Everything in this section is conditional on a complete recorded prefix Z.
Let `s=(s_f,s_s)` be the two exact first density scores and `T_ij` the
exact second density score from the fixed two-parameter family. Under the
uniform next-label law define

```
G_kl = E[s_k s_l],
M_ij,k = E[T_ij s_k],
alpha_ij = M_ij G^{-1},
phi_ij = T_ij-alpha_ij s.
```

The full census already establishes invertibility of G on every original00
prefix. This is a two-dimensional exact **source Gram**, not the singular
high-dimensional Monte Carlo covariance. No response estimate is used to
choose alpha. The source-normal score satisfies

```
E[phi_ij s_k]=0    for k=f,s,
E[phi_ij | class a]=0    for every joint-safe class a.
```

Consequently

```
H_ij F = alpha_ij H_1 F + E[phi_ij F].
```

The first term reflects source-space skewness acting through the existing
first responses. The second is the response to a distinct orthogonal
label perturbation. Both contributions were computed at each prefix
before population averaging; substituting population-averaged Grams or
response matrices would answer a different question.

For exact class counts n_a, summed marks S_i and d vacant labels,

```
M_ij,k = d^-4 sum_a sum_{u in a}
          (n_a L_i(u)-S_i)(n_a L_j(u)-S_j)(n_a L_k(u)-S_k).
```

All numerator arithmetic is integer. Equivalently
`M_ij,k=sum_a pi_a^4 kappa_a(L_i,L_j,L_k)`, whereas
`G_kl=sum_a pi_a^3 Cov_a(L_k,L_l)`.
The second-score covariance-subtraction term drops out of M because
each first score has zero within-class mean. These powers of pi differ
and must not be replaced by a single unweighted mark Gram.

## An operational Euler-invisible source

The normal score need not remain only a decomposition. It defines a new
finite, class-mass-preserving one-parameter source:

```
q_eta(u|Z) = pi_a exp[eta phi_ij(u)]
             / sum_{v in a} exp[eta phi_ij(v)]
```

inside each safe class, with outside probabilities unchanged. Its first
score at eta=0 is exactly phi_ij because every class mean is zero.
Both geometries' instantaneous rank/Euler joint law is preserved for
every finite eta. The measured positive own-center normal response is
therefore the derivative of an explicitly constructible perturbation
orthogonal to the first two loop scores. No eta trajectory was generated
in this analysis; the derivative uses the existing uniform samples.

Source orthogonality does not imply orthogonality or linear independence
of the future observer responses, and does not count continuum fields.

## What this excludes, and what it leaves open

If for every prefix the conditional future label mean had the form

```
m_F(Z,u)=c_a(Z)+b_f(Z)s_f(Z,u)+b_s(Z)s_s(Z,u),
```

with arbitrary class intercepts but a common pair of slope coefficients
across classes, every normal response above would vanish. The nonzero
measurements give evidence outside that description. They can arise from
different susceptibilities in different contact classes, or from higher
within-class mark structure. The present decomposition does **not** decide
between those two mechanisms. This is also a different target from the
previous approximately20% residual of a four-feature, across-prefix
loading projection; the two percentages must not be identified.

It also does not by itself refute a successful prediction of the first
Jacobian by `J(Z)=B G(Z)`. A component orthogonal to both first scores is
invisible to J. What is excluded is promoting that first-response relation
to the stronger full conditional-label-mean closure written above. The
PR509 first-response prediction task and this normal-response result can
therefore coexist and constrain different pieces of the mechanism.

Raw Hessians acquire first-derivative terms under nonlinear source-coordinate
changes. This normal component has a stronger transformation property.
For an invertible smooth reparameterization t=t(z), at the same base law,

```
T'_ab = J_ia J_jb T_ij + (partial_a partial_b t_i) s_i.
```

The first-score span is unchanged. Projecting out that span removes the
inhomogeneous last term, leaving
`phi'_ab=J_ia J_jb phi_ij`. Thus vanishing of the entire normal response
tensor is invariant under such a reparameterization of this same source
family. A named component and its sign still depend on the physical axes.
This finite source-space statement is not spacetime curvature or path memory.

## Reuse receipt

- Reader `8fa6fa01`; numerical result `2c3a5ca2`,
  [`score.json`](../results/p334-source-normal-curvature/score.json).
- Exact source Gram/class counts: `1cfa4ae8`; complete mark census and
  saved new64 first-response means: `8ad30617`; new64 second-response
  prefix means: `c48fa360`.
- The two prefix NPZ files retain exact cubic numerators, source Grams,
  coefficients, first/second responses and tangent/normal contributions.
- Each size keeps all72 scalar coordinates with aligned original20
  deleted-batch LOO/factors; no new prefix, suffix, fork-gzip pass, model
  fit, source-weight experiment, determinant or conditional-shape rerun.

Lifecycle: exploratory finite-size source-normal birth response / paired
N325,N425 original00 / original source population and already collected
new64 conditional streams / same shared dependency group / branch-delivered.
