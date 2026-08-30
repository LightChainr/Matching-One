# P250: spatial positive spectrum is not a hidden-field count

Date: 2026-08-30. Related issue: #406; upstream #250, #249, #370.

## Result and claim boundary

For the **complete, uniformly anchor-averaged raw P250 autocorrelation** on the
N101 parent of the N505 experiment, every nonzero parent momentum has strictly
positive spectral weight, for either hand and either measured charge, at every
Bernoulli parameter `0<p<1`. Consequently the full spatial moment matrix has
rank at least 100. This follows from a Fourier/cyclotomic proof and two exact
configurations of the actual N505 geometry, not a fitted spectrum.

This excludes a globally exact eight-state realization of that complete
endpoint series. It does not reject a useful low-rank approximation on a
finite measured window, imply 100 statistically resolvable modes, count CFT
fields, or exclude Jordan structure of a physical transfer/dilation operator.
The witness probabilities may be extremely small. No production covariance
or Monte Carlo sample was used to compute this result.

The canonical spatial shift is semisimple. A separate RG, thermal, or
connectivity action need not be. This distinction is the proposed change of
research variable, not a change to any frozen rank score.

## 1. Audited observable, not a generic toy kernel

Input semantics are pinned to
`33c557b9aebed1bc9c07019b9cd5cee6c04be947`:

- `scripts/z5_projective_leg_multiseparation_mc.py`, `ProjectiveLegIndex.scalar`:
  black NN rank-one membership is +1; white matching rank-one membership is
  -1; other component ranks give zero.
- `scripts/z5_projective_leg_cross_scale_mc.py`, `charged_rows`: five-fiber DFT
  of that real field, on children `19+12i` and `21-8i`, both of order 505.
- `scripts/z5_charged_threepoint_mc.py`, `dft_charges`: division by 5 and phases
  `zeta5^(-r*f)`. Thus opposite charges are complex conjugates.
- `scripts/z5_projective_leg_bivariate_mc.py`, `gauge_charged_rows` and
  `pair_value`: multiply by fixed fifth-root site phases, then evaluate
  `F_r(x)*F_(-r)(x+d)`.

The parent group is

```
H = Z^2 / <(10,1),(-1,10)> = Z/101,
j(x,y) = x-10y mod 101.
```

Uniform anchor averaging is part of the target below. The finite pseudorandom
sample itself is not claimed to be an exact product measure. In particular,
the implementation's modulo-101 anchor sampler has the usual tiny reduction
bias if interpreted as a uniform 64-bit draw; the exact theorem concerns the
intended uniform target or an explicitly complete anchor average.

## 2. Finite-group Bochner identity: proof

For a fixed hand and charge, let `f_omega(x)` be its fixed-gauge field. Define

```
C(d) = E_omega[ (1/101) sum_x f_omega(x) conj(f_omega(x+d)) ],
fhat_omega(k) = sum_x f_omega(x) exp(-2*pi*i*k*x/101).
```

Fourier inversion and character orthogonality give

```
C(d) = sum_k w_k exp(-2*pi*i*k*j(d)/101),
w_k = E[abs(fhat_omega(k))^2]/101^2 >= 0.
```

No translation covariance of the section-dependent field is needed: the
explicit anchor average supplies stationarity. The difference-indexed matrix
`T[u,v]=C(v-u)` is a positive Gram matrix. Its eigenvalues are `101*w_k` up to
frequency permutation. Periodicity and `C(-d)=conj(C(d))` follow as well.
The complete sum-indexed group Hankel matrix `H[u,v]=C(u+v)` differs from T by
row inversion, so it has the same rank. A *truncated* sum-Hankel matrix need
not be positive semidefinite and is not interchangeable with a principal
Toeplitz block.

The canonical spectral realization has

```
U = diag(exp(-2*pi*i*k/101)), V=U^(-10),
U^101=I, U^10 V=I, U^(-1) V^10=I.
```

It is commuting and unitary. More generally, an exact cyclic spatial action
whose minimal polynomial divides `t^101-1` is diagonalizable over C, since
that polynomial has no repeated roots. This conclusion does not refer to a
physical row-transfer matrix or scale generator.

The positivity applies to a fixed-hand, fixed-charge raw autocorrelation.
It must not be imposed on arbitrary signed hand differences, thermal
responses, or cross-field kernels. A Wilson-line/parallel-transport insertion
would also require its own explicitly changed observable contract.

## 3. Cyclotomic full-support lemma: proof

Let `K=Q(zeta5)`. The compositum of K and `Q(zeta101)` is `Q(zeta505)` and

```
[Q(zeta505):K] = phi(505)/phi(5) = 400/4 = 100.
```

Hence the degree-100 polynomial `Phi101(t)=1+t+...+t^100` remains irreducible
over K. For a nonconstant vector `f in K^101`, put
`P_f(t)=sum_{j=0}^{100} f(j)t^j`. If `P_f` vanished at any nontrivial 101st
root, irreducibility would force `P_f=c*Phi101`, including the possibility
`c=0`. That would make f constant, a contradiction. Thus **all 100 nonzero
Fourier frequencies of any such nonconstant vector are nonzero**.

This is not a numerical FFT tolerance argument. It is the elementary degree
argument over the *specific coefficient field* of the acquisition. It would
not be valid for arbitrary complex field values.

## 4. Two exact witnesses in the production geometry

In the plus child occupy the staircase of 19 east steps followed by 12 north
steps, excluding its repeated endpoint. In the minus child use 21 east steps
and 8 south steps. Both are closed by their declared period columns.

The generator uses adjugate-pair quotient keys and a lifted BFS. The separate
verifier uses axial cyclic labels and weighted union-find with integer cycle
displacements. It imports none of the generator's geometry/rank code.

| hand | occupied | black NN components | white matching components | nonzero parent field / zero |
|---|---:|---|---|---|
| plus, 19+12i | 31 | one, rank 1 | 474 vertices, one, rank 1 | 29 / 72 |
| minus, 21-8i | 29 | one, rank 1 | 476 vertices, one, rank 1 | 29 / 72 |

The last column holds for all four nontrivial deck charges. Use the legitimate
parent section `(j,0)`, `j=0,...,100`, with fiber points `(j+10f,f)`, `f=0,...,4`.
The certificate retains the occupied lifts, component sizes/ranks, exact zero
and nonzero values, and hashes of all cyclotomic coefficient vectors.
Coefficients represent `5*F_r` in the basis `1,zeta5,zeta5^2,zeta5^3`; the common
factor 5 does not affect Fourier support.

A change of section shifts the fiber cyclically and hence multiplies a
charged value by a fifth root. The repository's fixed C4 gauge does the same.
Neither operation changes its zero set. A vector with both zero and nonzero
entries remains nonconstant, so the support conclusion holds in the actual
C4 gauge without selecting a favorable gauge from the data.

Each witness has probability `p^m(1-p)^(505-m)>0` under independent Bernoulli
occupation. Nonnegative spectral weights cannot cancel its contribution.
Therefore `w_k>0` for every `k!=0`, and the complete spatial rank is at least
100. Whether the zero frequency is positive is not needed or claimed.

A useful corollary is that every proper principal Toeplitz block is strictly
positive definite: a null vector extended by zero to H would have vanishing
Fourier coefficients at all nonzero frequencies, so it would be constant on
H; its missing coordinates force that constant to vanish. Again, this is an
exact statement, not a lower bound on a statistically measurable eigenvalue.

## 5. A constructive alternative to endless rank escalation

**Archived data:** form a difference-indexed Toeplitz block on the ten
monomials of total degree at most three. All needed differences have Manhattan
radius at most six and are already in the recorded displacement domain.
Then fit/bound a nonnegative vector on the **known 101 spatial frequencies**:

```
E[y] = A w, w>=0.
```

Use the original raw moments and the complete covariance, not the already
truncated SVD projectors. The minus-hand R2 physical-coordinate map must be
applied before evaluating `j(d)`; do not copy a plus-hand abstract-coordinate
period constraint into the wrong chart. Partial data may admit many spectra.
Report positive-cone feasibility, bounds on spectral masses, and predictive
intervals rather than a unique reconstructed density without identification.
For a singular covariance, impose its exact null-space equalities separately.
A numerical cone distance needs statistical calibration; neither an optimizer
failure nor a naive chi-square degrees-of-freedom subtraction is a proof.

The exact endpoint audit gives 41/41 distinct parent vertices at radius four,
61/61 at radius five, but **85 labels and only 77 vertices at radius six**.
For example `(5,0)` and `(-5,-1)` differ by the parent period `(10,1)`.
Cross-campaign aliases remain independent noisy estimators of one mean, not
additional spatial points or magically paired replicas. Preserve their real
covariance and independent seed blocks when constraining equal means.

**Future acquisition efficiency, not new production here:** after constructing
the two homology indices, evaluate all 101 charged parent-field values. Average
all origins, or equivalently accumulate the nonnegative periodogram
`abs(fhat(k))^2/101^2`. One such array determines all 101 displacement means.
For the vector Y of random-anchor pair estimates and the configuration omega,

```
Cov(E[Y|omega]) <= Cov(Y)
```

in positive-semidefinite order, by the law of total covariance. This is a
Rao-Blackwell variance guarantee at a fixed number of configurations, not a
wall-time speedup claim. Benchmark its extra field/DFT work before claiming
better variance per second. Keep per-batch periodograms and joint covariance
across hand/charge. This turns repeated shell acquisition into one full-group
measurement and makes empirical PSD much easier to preserve.

## 6. A compressed commutator is not microscopic curvature

For commuting U,V and any projection P, Q=I-P,

```
[PUP,PVP] = PVQUP - PUQVP.
```

To prove it, insert `P=I-Q` between U and V on each side; the two uncompressed
products cancel. The exact three-dimensional control uses
`U=diag(1,-1,1)`, `V=diag(1,1,-1)`, and `P=I-11^T/3`.
Its compressed commutator is nonzero although U and V commute.

An ordered connected two-morphism experiment is still useful. It needs an
actual intermediate operation or observation; swapping two endpoint labels
cannot provide it. Before interpreting an effective commutator, distinguish
microscopic order, intervention-induced order, and eliminated-state leakage.
No such ordered P250 data are claimed here.

## 7. Research snapshot and arXiv-informed continuations

These are scoped suggestions, not a replacement control plane. Immutable
frontier anchors read in this contribution are:

- P250 `93eaab1640429d22353777b3756cf31049929024`: rank-eight R2 kernel-plane
  bridge rejected; it is not evidence of noncommuting translations.
- P321 `533ea6f8b659c4a0c3ae297c492971461de1d450`: deep equal-area pool retains
  the conditional E4 relation, but four 99% Fieller ratio sets are unbounded
  and the sole bounded D/C set includes zero. There is no identified empirical
  five-modulus subleading ratio curve.
- P321 `d0fca795c89338e16b732334c87d7885bc7fcfd0`: a graded two-closure sum is
  exact; its regular Q derivative is diagonal, not a forced Jordan extension.
- P334 `a251bd000cf1573e4a1cf9cd502466f581cab679`: the old one-mark reservoir
  fails at N9. Its image is the nonuniform sum
  `12*0+8*67+4*70+24*164=4752=11M`, against `16M` demand. Two-mark release
  reaches the full `48M` target image in these bounded rows. An arbitrary-HNF
  one-mark saturation conjecture must not be revived.
- P337 `cfb3eadcc5fb4b5a7017403696262eb9ba207c77`: finite-energy gluing works
  for globally typed, separated theta-six and figure-eight-eight events;
  ordinary arm colors alone do not determine global homology.

Three further proposals follow, separately from the proved spatial result:

1. **P321, test relations without dividing noisy amplitudes.** Once a physical
   two-closure thermal calculation supplies a candidate `f_j`, test the
   homogeneous constraint `D_N,j - rho_j*f_j*C_N,j=0` in the full joint C/D
   covariance. A free common coupling belongs in a profiled model image, not
   in post-reveal per-geometry ratios. Unbounded Fieller sets do not logically
   forbid every joint model test; they forbid claiming measured ratios. This
   supplies no missing F_t by itself and makes no power promise.
2. **P337, make global gluing probability the next unknown.** Study
   `P(typed theta closure | separated six arms)` and the corresponding typed
   eight-arm closure probability at fixed modulus/near-critical coordinate.
   The conjecture is a positive scaling limit on compact nondegenerate moduli,
   not an automatic consequence of finite energy. A first triangular control
   can use established near-critical/pivotal tools. The square NN/matching
   dictionary still requires its own proof or explicit universality assumption.
   Archive exterior component IDs and relative deck addresses when acquisition
   is warranted; the current threshold-only archive cannot invent them.
3. **P334, analyze target capacity before flow.** Treat the exact nonuniform N9
   mark-fiber signature as a symbolic capacity obstruction, and ask which
   incidence invariant predicts the next deficient fiber. Two-mark saturation
   in the checked rows is not a universal theorem. Extra bookkeeping labels
   do not create new unmarked target capacity.

## 8. Primary-literature scope ledger

arXiv pages/version metadata checked on 2026-08-30. These are selected relevant
sources, not a claim of an exhaustive survey or of reproducing every paper.

| source | useful input | limit of what is imported |
|---|---|---|
| [Yang-Xie-Stoica 1505.02510v3](https://arxiv.org/abs/1505.02510v3) | multidimensional positive Toeplitz/Vandermonde methods | their low-rank uniqueness condition is not assumed for a high-rank partial matrix; our finite-group identity is proved above |
| [Widder-Zimmer-Schilling 2503.20457v6](https://arxiv.org/abs/2503.20457v6), April 2026 revision | rigorous Mori projection/GLE and memory-kernel formulation for the specified setting | not a percolation field identification or a theorem for arbitrary projections |
| [Liu-Jacobsen-Saleur 2403.19830](https://arxiv.org/abs/2403.19830) | physical Jordan blocks can emerge from finite-size diagonalizable Potts/loop transfers | a different action from unitary finite-group spatial shifts |
| [Camia-Feng 2508.16047v2](https://arxiv.org/abs/2508.16047v2), June 2026 revision | critical triangular-percolation energy/log-partner lattice fields and scaling correlations | does not identify P250 or transfer rigor automatically to square site |
| [Ang et al. 2604.05503](https://arxiv.org/abs/2604.05503) | proposed exact critical-loop sphere three-point formulas with several validation routes | not a stand-alone citation for torus one-points or lattice singular vectors |
| [Roux-Ribault-Jacobsen 2604.24491](https://arxiv.org/abs/2604.24491) | torus one-points via sphere four-points at a different central charge; modular bootstrap | does not provide P321's missing lattice closure calibration |
| [Ikhlef-Morin-Duchesne 2602.15742](https://arxiv.org/abs/2602.15742) | critical ADE RSOS/TL modules, local connectivity operators and lattice difference relations | an ADE construction strategy, not an already established square-site thermal-Q4 operator |
| [Jacobsen-Ribault-Saleur 2208.14298v3](https://arxiv.org/abs/2208.14298v3) | twisted O(n)/Potts state spaces and diagram-algebra branching | representation information, not equality of a formal trace with physical homology weights |
| [Cardy 2201.00478v2](https://arxiv.org/abs/2201.00478v2) | TTbar-deformed modular forms and one-point structure | not a derivation of unknown homology thermal amplitudes |
| [He-Sun 2004.07486v2](https://arxiv.org/abs/2004.07486v2) | first-order TTbar torus correlation/Ward framework | does not remove the sector-resolved product-rule obstruction |
| [Grans-Samuelsson et al. 2302.08168v2](https://arxiv.org/abs/2302.08168v2) | combinatorial-map-labeled loop correlators and tested bootstrap-basis conjecture | the critical completeness statement is a conjecture, not an automatic dictionary |
| [Garban-Pete-Schramm 1305.5526v4](https://arxiv.org/abs/1305.5526v4) | triangular near-critical/dynamical scaling limits from pivotal measures | does not prove square-site typed torus gluing |

## 9. Reproduce and validate

```
python3 scripts/p250_spatial_spectrum_certificate.py --output /tmp/p250-spectrum.json
python3 scripts/verify_p250_spatial_spectrum.py /tmp/p250-spectrum.json
python3 -m unittest discover -s tests -p test_p250_spatial_spectrum_certificate.py -v
```

The focused suite has 19 passing tests, including deliberately corrupted
occupation, topology, charged coefficients, alias counts, and conclusion
boundaries. Python byte-compilation also passes. The full repository suite
and production scores were not run in this contribution. The machine verifier
checks the finite witnesses independently; the written Fourier/cyclotomic
argument supplies the theorem rather than pretending a finite scan proves it.

No remote execution environment was contacted. All additions are independent
new files; existing results, freezes, navigation and branch histories remain
unchanged. The positive-cone archived-data scorer and a periodogram acquisition
benchmark are next analyses, not completed results of this contribution.
