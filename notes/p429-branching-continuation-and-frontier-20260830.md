# Matching One: trace prediction is not autonomous state closure

Date: 2026-08-30. Research issue: [#429](https://github.com/LightChainr/Matching-One/issues/429).
Status: **exact finite-volume counterexample, independently checked; scaling programme is a conjecture**.

## Decision

Do not extend the cooperative-trigger ladder merely from c3 to c4. Test whether a proposed current-state summary preserves **branching continuations**, not only an unbranched lifetime distribution. We now have two rank-one configurations whose complete lifetime laws agree but whose next predictive-state distributions differ.

The distinction matters only when an autonomous coarse Markov process, intermediate state-conditioned operation, or intervention is claimed. The full survival vector remains an exact and useful predictor of the original rank-only future trace. A more complex state is unnecessary for that narrower objective.

## Current repository boundary

The initial main snapshot was `c64bfabde8c8bb13290d4f5d8e5f44c1779d4d30` (13:13:48 UTC), including the N10 killed-survival certificate. The final branch base is `f854f823a90be3395178fdfce04f1bae21435d46` (13:46:07 UTC). The intervening two commits add adjacent-annihilator asymptotics and exact Gaussian channel fingerprints; they do not change this rank oracle or continuation observable.

Read alongside the active research discussions:

| Existing frontier | Consequence for the next analysis |
|---|---|
| [#401](https://github.com/LightChainr/Matching-One/issues/401), [#403](https://github.com/LightChainr/Matching-One/issues/403), PR #415 `d09f925`: two-/three-site completion, overlap terms, killed survival | Do not re-propose c2, c3 or full survival. The new question is closure of the *recomputed* signature under updates. |
| [#398](https://github.com/LightChainr/Matching-One/issues/398#issuecomment-5468948610): width-four gr1 filled by seven rooted coordinates, rank 13/13, determinant 3072 | The next module calculation is the affine/endpoint/radical/Gram/source closure of that fixed registry, not another list of marks or deeper Q orders. |
| [#406](https://github.com/LightChainr/Matching-One/issues/406), [#418](https://github.com/LightChainr/Matching-One/issues/418): spatial Fourier positivity and exact section/phase-mask structure | Diagnose coordinate and gauge effects before reading additive spatial rank as a field count. This is not the same generator as occupation-rank evolution. |
| [#405](https://github.com/LightChainr/Matching-One/issues/405), [#408](https://github.com/LightChainr/Matching-One/issues/408): theta/figure-eight census and typed-arm acceptance | The remaining question is the scale dependence of the *typed acceptance probability*. Local surgery and a bare six-arm exponent do not determine it. |
| [#419](https://github.com/LightChainr/Matching-One/issues/419): observer degree/Walsh bandwidth | Cheap low-rank response may reflect the observer and generator. Do not interpret every small state count as a physical field space. |

These developments belong to their existing authors/threads. They are not new claims of this note. The contribution here is the exact trace-versus-branching counterexample and its acquisition consequence.

## Exact N16 witness

Use P=diag(4,4), with physical coordinates modulo four:

```
A={(0,0),(1,0),(2,0),(3,0),(1,1),(3,1),(0,3),(1,3)},
B={(0,0),(1,0),(2,0),(0,1),(1,1),(2,1),(3,1),(0,3)}.
```

Both have k=8, rank one and primitive ambient line (1,0). With d=8 vacancies, let

```
b_m(S)=#{U subset vacant(S): |U|=m, r(S union U)=1},
s_m(S)=b_m(S)/binom(d,m).
```

The complete count vector is the same:

```
b=(1,7,18,20,8,0,0,0,0),
s=(1,7/8,9/14,5/14,4/35,0,0,0,0).
```

Independent enumeration of all 8! future permutations for each witness gives the same counts for exiting at steps 1,...,5:

```
(5040,9360,11520,9792,4608).
```

Thus this is not a missing higher-order coefficient. Every unmarked future rank trace has the same distribution.

After one common uniformly random vacant insertion, the number x' of one-site exits at the successor has the following choice counts:

| initial | absorbed | x'=1 | x'=2 | x'=3 |
|---|---:|---:|---:|---:|
| A | 1 | 3 | 2 | 2 |
| B | 1 | 1 | 6 | 0 |

Clone that successor and add one independently chosen vacant site in each clone. Both remain rank one with probabilities

```
B_11(A)=(3*36+2*25+2*16)/392=95/196,
B_11(B)=(36+6*25)/392=93/196,
difference=1/98.
```

There are 190 versus 186 successful choices out of 392; branch choices may select the same site across the independent copies. The unbranched two-step survival is 9/14 for both. If cloning occurs immediately at A/B, all joint rank-only probabilities are products of equal survival probabilities and cannot distinguish them. The shared update is essential.

## Failure under the actual uniform-permutation law

Strong non-lumpability alone need not imply failure under one particular initial law. We therefore checked the actual law separately.

Condition at k=8 on line (1,0) and the full b-vector above. There are 192 current subsets and 192*8!=7,741,440 selected prefixes. Let E mean that the ninth configuration remains rank one and has x'=3. The exact conditional probabilities are:

| first birth K1 | selected prefixes | P(E) |
|---|---:|---:|
| 4 | 110592 | 1/6 |
| 5 | 442368 | 1/6 |
| 6 | 1198080 | 2/13 |
| 7 | 2442240 | 8/53 |
| 8 | 3548160 | 2/11 |

All five cohorts have the same next-step absorption probability 1/8. Nevertheless the transition to a future signature with x'=3 differs by 1/66 between K1=8 and K1=4. K1 is measurable from the past rank observations; E is measurable from the next survival signature. Hence the **recomputed microscopic signature process** is not Markov under uniform growth.

This does not refute a filtered belief-state representation based only on rank observations. In that representation the observer shifts and renormalizes a survival curve rather than learning the actual hidden successor's signature.

## The bounded refinement calculation

The Python analysis uses the identity

```
m*b_m(S)=sum_{safe one-site children T of S} b_(m-1)(T).
```

Each safe m-set is counted once for each of its m possible first sites. Exact integer divisibility is checked. Strong Markov classes are computed separately, backwards in occupied-count order: configurations have the same class precisely when they have the same k/line labels and the same multiplicities of successor classes. Rank-two successors are one cemetery state. Induction over layers proves this is the coarsest strong Markov partition preserving the declared labels.

| quotient | rank-one states | survival classes | strong Markov classes | split survival classes |
|---|---:|---:|---:|---:|
| 2+i | 10 | 2 | 2 | 0 |
| 3 | 162 | 10 | 10 | 0 |
| 3+i | 310 | 16 | 16 | 0 |
| 3+2i | 2340 | 62 | 62 | 0 |
| 4 | 19932 | 210 | 214 | 4 |
| 4+i | 38896 | 346 | 390 | 42 |

Counts exclude the cemetery. These six quotients are the entire reported search, not an exhaustive all-HNF minimality claim. The class counts are not continuum field counts.

## From a linear predictor to an observable algebra

For a finite stochastic kernel P, a partition C is strongly lumpable iff the algebra of block-constant functions is P-invariant. Proof: the algebra has basis 1_C over blocks C; P1_C(S) is exactly the probability of transitioning into C. Invariance is exactly constancy of every such probability on each source block.

The minimal autonomous state therefore corresponds to the smallest unital pointwise algebra containing the declared observation functions and closed under P. This is different from a linear space containing P^m 1_alive. Taking products corresponds to independent forks of the same microscopic state; applying P corresponds to an ordinary shared update.

For successor survival functions f_i,

```
Gamma_ij=P(f_i*f_j)-(P f_i)*(P f_j)
```

is a conditional covariance matrix, so it is positive semidefinite. For the displayed witness, the relevant Gamma entry is 1/14 at A and 3/49 at B. Their difference is again 1/98.

No abstract novelty is claimed for the finite-state algebra/lumpability criterion. The new contribution is its exact percolation witness and a way to measure its physical consequence. This pointwise real function algebra is reduced, and is not #398's potentially nonreduced Hankel or Q-adic object. Identical words such as rank, algebra or nilpotent do not make the generators or representations identical.

## An error-controlled compression target

An approximate block map should be evaluated on successor probability laws, not only mean predictions. A simple sufficient finite-horizon guarantee is available without importing a continuous-time theorem. Suppose a declared representative coarse kernel Q_k satisfies

```
TV(pushforward(P_k(S,.)), Q_k(Z(S),.)) <= epsilon_k
```

uniformly over the tested microscopic states in every block. Successive maximal couplings imply a mismatch probability at most

```
1-product_k(1-epsilon_k) <= sum_k epsilon_k
```

for the finite coarse trace, provided observation labels agree at coupled states. For a finite branching experiment, a union bound applies over its transition edges; the shared prefix is counted once. This is a uniform guarantee, not something established by a small average regression residual.

A practical next experiment estimates a few shared-prefix/fork kernels at frozen scaled horizons, stratified by a declared geometry summary. Current full-survival equality is an exact null for root forks, but not for delayed forks. The new N16 control must pass before a large pilot is interpreted.

At large N, near-critical windows, geometry injection radius, current-state snapshots and covariance domains must be typed. The candidate N^(5/8) site-count horizon assumes the usual two-dimensional thermal scaling; an independently calibrated intrinsic clock avoids making the exponent assumption part of the basic test. Clones, shared prefixes and initial checkpoints form a nested sampling hierarchy and are not independent configurations.

## Literature scan: what can be imported, and what cannot

The scan covered current arXiv entries across probability, process abstraction, algebraic realization and LCFT. Relevant theorem sections were read for behavioural witnesses, approximate bisimulation and hierarchical percolation; other entries supply source-checked scope and directions, not unverified imported theorems.

1. **Turkenburg--Beohar--van Breugel--Kupke--Rot**, *Constructing Witnesses for Lower Bounds on Behavioural Distances*, [arXiv:2504.08639v2](https://arxiv.org/abs/2504.08639), revised 2025-10-13. Constructive quantitative modal witnesses suggest that a failed coarse state should return a small distinguishing experiment, not just a residual. The present 1/98 is directly enumerated, not an invocation of their metric normalization.

2. **Spork--Baier--Katoen--Klueppelholz--Piribauer**, *Approximate Probabilistic Bisimulation for Continuous-Time Markov Chains*, [arXiv:2505.15587v2](https://arxiv.org/abs/2505.15587). Transition-law error and clock-rate error are separated; bounded-reachability errors are controlled. The current model is layered discrete growth, so their CTMC bounds need a declared time/model map. The discrete coupling bound above is derived separately.

3. **Chen--Clerc--Panangaden**, *Two behavioural pseudometrics for continuous-time Markov processes*, [arXiv:2511.21621](https://arxiv.org/abs/2511.21621), submitted 2025-11-26. Transition-sensitive and path-sensitive distances are genuinely different targets. This is useful language for the trace/branching distinction, not a torus result.

4. **Garban--Pete--Schramm**, *The scaling limits of near-critical and dynamical percolation*, [arXiv:1305.5526](https://arxiv.org/abs/1305.5526). Full geometric near-critical states and pivotal updates provide the continuum context. Their Markov construction does not imply a finite rank/line or reliability quotient is Markov.

5. **Widder--Zimmer--Schilling**, *On the generalized Langevin equation and the Mori projection operator technique*, [arXiv:2503.20457v6](https://arxiv.org/abs/2503.20457), revised 2026-04-24. Projection-induced memory motivates asking which current geometry was discarded. The exact finite result here does not identify a GLE kernel or asymptotic memory law.

6. **Alves--Baldasso--Moreira--Teixeira**, *Percolation on hierarchical lattices*, [arXiv:2606.11503v1](https://arxiv.org/abs/2606.11503), submitted 2026-06-09. Edge-replacement graphs under explicit seed hypotheses admit rigorous critical and near-critical analysis. A small hierarchical control could test whether the proposed abstraction machinery detects known recursive structure. It would calibrate methodology, not approximate square-site exponents by assertion.

7. **Balle--Panangaden--Precup**, *A Canonical Form for Weighted Automata and Applications to Approximate Minimization*, [arXiv:1501.06841](https://arxiv.org/abs/1501.06841). Minimal linear realizations address a declared string/test language. Extending that language to forks is not accomplished by renaming the old Hankel matrix.

8. **Martin--Senecal--Spencer**, *Cell modules for the Temperley-Lieb algebra in mixed characteristic*, [arXiv:2601.17445](https://arxiv.org/abs/2601.17445), submitted 2026-01-24. Submodule structure and Jantzen-like layers motivate #398. The repository's measured J2=0 should take precedence over an expectation of deeper layers.

9. **Bernardi--Jelisiejew--Reig Fite**, *Hankel and Multiplication Tensor Completions for Cactus Rank*, [arXiv:2606.30600](https://arxiv.org/abs/2606.30600), submitted 2026-06-29. Flat extension and multiplication-tensor completion are equivalent under the stated algebraic setup. This does not make every critical translation kernel a finite Artinian realization.

10. **Giaquinto--Mastnak**, *The Center of the Temperley-Lieb Algebra*, [arXiv:2607.12247](https://arxiv.org/abs/2607.12247), submitted 2026-07-14. Cellular filtrations and deformation give exact centre information. They concern a specified TL algebra; rooted connectivity closure and physical logarithmic coupling still require their own maps.

11. **Camia--Feng**, *The percolation energy field and its logarithmic partner*, [arXiv:2508.16047v2](https://arxiv.org/abs/2508.16047), revised 2026-06-01. Concrete lattice logarithmic observables are a stronger field-identification control than a free logarithmic fit. They do not equate a finite probability-state algebra with a Virasoro module.

12. **Paul Roux--Ribault--Jacobsen**, *Torus one-point functions in critical loop models*, [arXiv:2604.24491](https://arxiv.org/abs/2604.24491), submitted 2026-04-27. Modulus-sensitive one-point structure gives an orthogonal physical fingerprint. The first author's name is Paul, not Augustin as transcribed in an older issue.

## Reproduction and provenance

```
python scripts/p429_branching_continuation.py --output /tmp/p429-exact.json
python -m unittest discover -s tests -p 'test_p429*.py' -v

g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic \
    scripts/verify_p429_n16.cpp -o /tmp/verify-p429
/tmp/verify-p429 /tmp/p429-n16-ranks.txt
```

The imported Python topology oracle is **byte-identical** to Git blob `62e06795fdfa91a956aedd62b7344e84aa5efc5c` at commit `fee33287`. It is not another topology implementation or a changed production observable. The separate C++ verifier uses row-major coordinates, potential union-find, direct subset counts and selected-prefix permutation enumeration, with no prefix DP. All 65,536 N16 rank/line states agree after exact relabeling. The selected-prefix enumeration covers 7,741,440 prefixes; the future-trace tests additionally enumerate 2*8! suffixes.

Focused validation: 13 tests passed, including the C++ check; no skipped test in the recorded run. One development test initially exposed a tuple-versus-JSON-list serialization mismatch; the schema was normalized without changing any scientific count. No repository-wide CI or production simulation was run.

## Claim boundary

The new exact statements are finite trace equivalence, failure of update lumpability, the 1/98 branching witness, the 1/66 actual-law history witness, and the bounded class census. The proposal is to study approximate branching-aware closure at meaningful scaled horizons. No continuum field, asymptotic memory exponent, universal state dimension or minimal HNF obstruction is identified.
