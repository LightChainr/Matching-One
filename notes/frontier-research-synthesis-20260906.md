# Frontier research synthesis after #579/#581/#582

**Date:** 2026-09-06  
**Scope:** research synthesis and experiment/theory design only. No new lattice evidence, no issue-lifecycle claim, no production authorization by implication.  
**New focused follow-ups:** #584, #585, #586.  
**Existing issues amended rather than duplicated:** #249, #263, #580, #583.

## Executive view

The latest repository results and the 2024–2026 percolation/loop/CFT literature point in the same direction.

The project should stop asking first whether the observed finite-size state is exactly rank 1, rank 2, one Jordan pair, one modular form, or one extra scalar correction. The better working object is now

```text
low-dimensional transferable base
      +
small context/topology/cover/combinatorial fiber,
```

with several distinct notions of rank depending on which experiment category is being represented.

Three repository results make this more than an aesthetic reframing:

1. **#579 / N=580:** denominator-free projective inference exposes a negative aspect-ratio curvature that no frozen ray predicts. The old `A8/A4 << 1` assumption was not a direct angular-amplitude measurement. The next angular experiment must therefore measure the nuisance harmonic rather than assume it.
2. **#582:** the full threshold law has a genuine non-affine shape flow. One held-out-transferable shape direction removes about 99% of the covariance-weighted flow, but a small, highly resolved, structured remainder remains and does not close at rank 2.
3. **#581:** the exact square-bond critical-Q tangent splits into a duality-even Betti component and an ambient-homology component, and these two source types are already spectrally different at Boolean degree one on tiny exact controls.

The external frontier simultaneously says:

- torus loop correlations are organized by **modular-covariant solution spaces plus combinatorial maps**, not by one modular function per field;
- exact charged three-point/OPE data are now available for generic loop fields;
- annulus propagation can contain infinitely many real/complex/logarithmic spectral contributions;
- a known percolation logarithmic pair exists as an explicit lattice object;
- continuum Jordan blocks can emerge while every finite transfer matrix remains diagonalizable;
- finite noninvertible symmetry/defect sectors have a natural tube-algebra/open–closed language;
- noisy linear Hankel rank and positive/Markov realization dimension are distinct inference problems.

The recommended program is therefore to **type the state before naming the field**.

---

## 1. What the N=580 modulus anomaly now means

The important result of the N=580 reanalysis is not that `bare_aspect_ratio` is a winner.

A model that predicts proportions predicts a ray. Once the full three-rung raw vector is scored against each ray with the measured covariance, seven of eight frozen competitors are strongly excluded. More importantly, the linear second divided difference

```text
f[1,2,4] = (m(4) - 3 m(2) + 2 m(1))/6
```

is negative at roughly three standard deviations, while the frozen linear families predict zero and the frozen modular families predict positive curvature.

Thus the current physical statement is simply:

> the three-point aspect response is concave, and the frozen model list is not.

The natural objection is spin-8 leakage. Solving the per-rung leakage signs exactly gives required `|A8/A4|` values of order 3–18 for the surviving nonconstant families, not the much larger leading-order numbers first reported. The qualitative conclusion survives: the nuisance harmonic would have to be comparable to or larger than the nominal spin-4 amplitude for most candidates.

But the load-bearing assumption was never a direct measurement: the old H4-beats-H8 result was a **homology-character model comparison**, not an angular Fourier-amplitude measurement.

### Consequence for #583

N=650 with three orientations per family is worth doing because it directly estimates

```text
C, A4, A8
```

with one correlated response block.

However, three angles and three coefficients interpolate exactly. They **cannot validate** the truncated angular law

```text
C + A4 cos(4 theta) + A8 cos(8 theta)
```

without an additional constraint.

Therefore split the scientific questions:

```text
N=650:
  measure A8 instead of assuming it;

held-out angular closure:
  add one overdetermining angle/geometry at the cheapest arithmetic-compatible size.
```

This correction is recorded on #583. A fourth angular coordinate is more informative than another long aspect-ratio ladder if the extraction of `A4` itself is not yet validated.

If `A8` is genuinely large and the overdetermined angular closure passes, rescore the old modulus data with the measured harmonic content. If `A8` is small and closure passes, the N=580 concavity becomes a cleaner continuum/modulus problem. If closure fails, stop interpreting the fitted `A4/A8` coefficients as physical harmonics until the angular representation is redesigned.

---

## 2. #582 changes the state question

The whole-law Wasserstein/quantile analysis is now one of the most informative objects in the repository because it starts from the complete one-dimensional threshold distribution rather than a selected derivative or moment.

The result is not W0/W1/W2/Wbroad in the original decision table.

After projecting the exact affine tangent `span{1,Q_N}`:

- center + width alone are overwhelmingly rejected;
- the non-affine part is only a few percent of the raw displacement but is resolved at very high signal-to-noise;
- one shape direction learned without the held-out transition removes roughly 99% of the covariance-weighted flow;
- random directions do not come close to that performance;
- the remainder remains highly significant and structured;
- pairwise remainder angles form reproducible clusters that are not explained by the already-corrected spin-4 residue, production campaign, scale multiplier, or monotone size.

The scientific object is therefore better written as

```text
v_shape(context)
   = a(context) g_base
   + r_fiber(context),
```

where `g_base` is strongly transferable and `r_fiber` is small but not noise.

The mistake would now be to add an unconstrained third PCA vector and call the result rank 3.

### New #584

#584 asks whether the remainder is indexed by labels that already existed before the residual was examined:

```text
Smith type,
deck group,
Gaussian cover word,
parent lineage,
primitive homology data,
period/orientation class,
other declared context metadata.
```

The test is permutation-controlled alignment/classification of the residual directions, not a free exponent fit.

If one label predicts the residual, freeze a **base + typed fiber** representation and take it to #580. If no label beats permutation/random controls, add at most one or two crossed transitions designed to swap the leading candidate label while holding the rest fixed.

### Why this is compatible with the predictive-state no-go theorems

The cut-network results #435/#491/#549/#550 prove that increasingly rich branching continuations can require growing predictive classes even when simpler survival summaries agree.

That is not in conflict with a low-dimensional marginal threshold-law flow.

A coherent picture is

```text
simple endpoint/marginal base
+
context information exposed only by richer continuations.
```

Thus a successful #584 does not imply a finite universal state dimension. It identifies a useful factorization of one declared experiment family.

---

## 3. #581 gives a physical type system for finite tangent information

On exact-critical square-bond tori, the critical-manifold Q score obeys

```text
T - T* = X,
T = V/2 + B/2,
B - B* = 2X,
```

and hence

```text
T = constant + B_even/2 + X/2.
```

The new exact gate is especially valuable because it proves that `B_even` and `X` are distinguishable information channels rather than two arbitrary coordinates:

```text
B_even: all Boolean degree-1 coefficients vanish;
X:      degree-1 coefficients are uniform and nonzero.
```

On the tiny L=3 control the ambient channel already puts a large fraction of its spectral weight at degree 1, while the duality-even Betti part begins at higher degree and has a different tail.

This does **not** identify `B_even` with the local energy field or `X` with a continuum topological defect. It does something more modest and more useful: it supplies a finite-lattice type system that can be tested before continuum naming.

### Immediate empirical route

Continue #581 on a square-bond control with the martingale/noise decompositions already proposed in #256/#227:

```text
where does B_even become predictable?
where does X become predictable?
which declared observable couples to which scale profile?
```

Then use only the intrinsic square-site `X_site=r-1` on the square-site side unless a separate bulk coordinate has an explicit definition/lift.

The highest-information possibility is that the #582 base and remainder acquire different scale fingerprints. That would turn the recurrent “two-state” language into a concrete finite statement such as

```text
base  -> bulk/mesoscopic signature,
fiber -> ambient-topology/quotient signature,
```

or the reverse.

---

## 4. New organizing hypothesis: a graded/fibered representation

The earlier proposal of a unified finite representation should now be sharpened.

Instead of one tensor product chosen in advance, treat the experiment category as graded by physically distinct labels:

```text
Bernoulli / Krawtchouk / threshold-law coordinate
x spatial / annular information scale
x angular irrep
x deck / Smith / homology character
x combinatorial-map / connectivity sector
x intervention / generator label
x insertion order / Q-tangent grade.
```

A useful schematic is

```text
low-dimensional continuum-like base
       ⋉
discrete/topological/combinatorial fiber.
```

The semidirect-product symbol is only mnemonic: no group action is claimed. The point is that a small base can coexist with nontrivial context memory.

This hypothesis makes several old results mutually compatible:

- scalar closure repeatedly fails;
- a low-dimensional direction can nevertheless transfer extremely well;
- observer/generator changes can rotate the residual;
- endpoint Hankel rank can be small while spatial or branching predictive complexity is high;
- positive/reversible microscopic dynamics can mimic low-order/Jordan-like scalar traces;
- topology/deck/map information may first appear only in marked/twisted/multi-point channels.

The next measurements should determine **which grading is necessary**, not how many arbitrary latent vectors can be fit.

---

## 5. Map-resolved torus tomography replaces a list of modular guesses

### External frontier

Roux--Ribault--Jacobsen, **arXiv:2604.24491**, construct torus one-point functions in critical loop models from sphere four-point data. The relevant lessons for this project are structural:

- modular covariance does not select one universal torus shape;
- even simple primary families admit multiple modular-covariance solutions;
- torus one-point functions are infinite non-chiral conformal-block combinations;
- logarithmic blocks are included;
- the correlation-function definition includes a **combinatorial map/connectivity pattern** in addition to the surface and field labels.

Grans-Samuelsson et al., **arXiv:2302.08168**, already developed combinatorial/ribbon maps as an organizing basis for loop-model correlation functions.

This strongly supports #249's L3 possibility: apparent extra state may be map memory rather than another local bulk primary.

### New #585

#585 replaces the next round of isolated `E4`, `r^a`, single-log-block guesses by a covariance-weighted fit to **symmetry-allowed map-resolved solution spaces**.

The sequence is:

```text
lattice observable semantics
 -> allowed map/representation sectors
 -> modular-covariant solution subspace
 -> projective/subspace score on existing moduli
 -> held-out modulus/readout
 -> sphere four-point / OPE cross-check.
```

The amplitudes are nuisance coordinates inside a physically declared subspace. The subspace itself must be fixed by field/map/symmetry semantics, not enlarged until it fits.

### Why #156 and #244 are useful here

- #156's C3 primitive-homology characters can eliminate torus sectors by representation content instead of scalar shape alone.
- #244 forbids nontrivial deck-character linear response in an unmarked invariant scalar, so a solution requiring such a charge cannot be assigned to the wrong lattice observable.
- #576 keeps standard Pinson wrapping as a semantic control but not an explanation of the N=290 anomaly.

If no defensible map dictionary can be supplied for the current scalar observable, the right next experiment is a marked/connectivity-resolved observable, not another scalar torus.

---

## 6. Exact three-point data makes OPE spectroscopy realistic, but the old charged tables should not simply be rerun

Ang et al., **arXiv:2604.05503**, propose exact normalized three-point constants for generic critical-loop primaries `V_(r,s)`, where for `r>0` the field carries `2r` legs and `s` is a momentum/phase for the winding legs. The formulas are checked by bootstrap, transfer-matrix calculations and CLE/LQG methods.

Jacobsen--Nivesvivat--Ribault--Roux, **arXiv:2510.04701**, also develop nonzero-spin three-point functions and show that lattice extraction can be complicated by degeneracies in the finite Jones--Temperley--Lieb modules.

The repository's #250 program was directionally correct, but the existing broad charged/spatial tables already revealed richer finite rank and should not be reset to the original simple R2 story.

The improved route is:

1. let #585 or another representation analysis identify a specific map/charge sector;
2. construct a lattice insertion faithful to the generic-loop `V_(r,s)` prescription;
3. measure a normalization-invariant three-point ratio;
4. compare

```text
spin/dimension
+ representation tensor
+ combinatorial map
+ normalized OPE coefficient.
```

This is far more identifying than another one-point amplitude ratio.

If the field depends on Q, its Q tangent must include both measure and explicit field/OPE/projector derivatives; that belongs to #586/#263, not to a raw score covariance alone.

---

## 7. Graded Q->1 tangent CFT should restart through positive controls, not through Matching One

The old #263 architecture remains the right formal decomposition:

```text
partial_Q correlator
 = measure derivative
 + projector/invariant-tensor derivative
 + field/normalization derivative
 + OPE derivative
 + dimension/central-charge/block derivative.
```

What was missing was an executable calibration sequence.

### New external controls

Gefei Cai, **arXiv:2603.28161**, gives exact boundary four-point functions for `CLE_kappa`, `4<kappa<8`, including critical percolation at `kappa=6`. The whole cross-ratio function lives in a differentiable one-parameter conformal family. One can differentiate the exact function and the BPZ/fusion differential equation itself, yielding an inhomogeneous tangent ODE.

Camia--Feng, **arXiv:2508.16047**, provides explicit triangular-lattice energy/four-arm fields whose two- and three-point scaling limits realize a logarithmic pair.

Together with #581 and the exact 2026 three-point constants, this creates a clean positive-control ladder.

### New #586

#586 is intentionally narrower than reopening all of #263:

```text
A. exact CLE boundary Q tangent;
B. lattice square-bond/triangular estimator split into B_even and X;
C. known bulk logarithmic-pair graded signature;
D. one complete projector/OPE differentiated crossing check;
E. only then spin-4 Matching One candidates.
```

The output is a graded tangent table, not “Jordan yes/no”.

This route is especially valuable because #333 has already shown that generic-Q lifts are not unique. A successful square-bond tangent is a positive control, not a unique continuation of the square-site observable.

---

## 8. Annulus/radial spectroscopy remains important, but high rank is no longer a failure condition

Sun--Xu--Zhuang, **arXiv:2410.04767**, derive exact annulus crossing laws in which the two-open-path channel has a transcendental leading backbone exponent and a countable family of complex subleading exponents; the thin-annulus expansion contains logarithmic corrections at every order.

This changes how #253 should be interpreted.

A radial state that needs several modes, rotates, or fails a two-state Jordan fit may still reflect real continuum spectrum rather than poor finite-size physics.

Therefore distinguish at least

```text
r_endpoint  rank seen by endpoint/context Hankel data,
r_radial    rank needed for composable annulus propagation,
r_positive  positive hidden-state dimension,
r_CFT       representation/operator-sector complexity.
```

Do not collapse these into one `rank`.

The annulus route should only be revived through an actual composable boundary/cut-network state satisfying a test of

```text
U(t1+t2) ?= U(t1) U(t2),
```

not another static rectangle or shell-constant fit.

A useful constrained two-mode diagnostic from #253 remains the recurrence discriminant:

```text
Delta > 0  two real modes,
Delta = 0  repeated/Jordan boundary,
Delta < 0  conjugate complex pair.
```

But an eventual failure of all low-order recurrences should now be interpreted against the exact infinite annular spectrum rather than as automatic evidence against continuum structure.

---

## 9. Realization methodology: weighted linear state, positive state and intervention-stable state are different objects

### Weighted realization

He et al., **arXiv:2505.19639**, revisit realization from range-space and null-space least-squares viewpoints and motivate weighted least-squares estimators. Their assumptions do not directly match the repository's correlated jackknife blocks, but the lesson is straightforward: when the covariance of context-Hankel entries is already available, it should enter the realization fit, not merely a post-hoc singular-value plot.

#249 has been amended accordingly.

### Positive realization

Taghavian--Sjolund, **arXiv:2502.21102**, study minimal positive Markov realizations. Positive realization is a stricter problem than an unconstrained linear input/output realization.

This is directly relevant to P398/#580, because the microscopic control is positive/reversible while its scalar traces can mimic nonnormal/Jordan phenomenology.

Report separately when meaningful:

```text
r_linear,
r_positive,
r_transport.
```

A positive hidden-state model may require more states than a signed linear realization; conversely, a low linear Hankel rank is not a physical state count.

### Intervention-stable state

#580 asks the strongest finite test in this hierarchy: freeze the state/readout functions in one environment and change only the physically intervened generator.

After the exact P398 calibration, the first square-site candidate should no longer be an arbitrary low-rank basis. It should consume

```text
#582 common base
+
#584 typed fiber if one survives,
```

and compare it against fresh per-intervention refits.

A state that transports controlled interventions is a substantially stronger finite object than a same-environment compression, while still not proving latent/continuum identity.

---

## 10. Defect/tube-algebra route remains a serious alternative to bulk-field naming

Choi--Rayhaun--Zheng, **arXiv:2409.02159**, develop generalized tube algebras and symmetry-resolved partition functions for finite noninvertible symmetries, including twisted boundary states and generalized open–closed duality.

This is a natural external framework for the repository's older #244/#233/#321 cluster:

```text
#244 exact charged-mode null for invariant scalars,
#233 possible Q->1 derivative defect,
#321 torus homology balance vs cylinder/TL sector crossing.
```

The key change from the old language is that a matching-odd object need not be an ordinary `Z2` parity eigenvector. It may be a matrix element or tangent in a noninvertible defect/tube sector.

The first serious test should still be a positive control on square-bond Potts/TL where defect pull-through/fusion is exact. Only then should one ask whether a local square-site/matching interface exists.

A proof of nonexistence under a fixed locality/bond-dimension class would be scientifically valuable and would push the interpretation toward intrinsic homology/nonlocal map structure.

---

## 11. Previous recommendations retained, revised, or deprioritized

### A. Unified finite representation — **retained and sharpened**

Old form:

```text
Bernoulli/Krawtchouk x spatial/deck character x angular irrep.
```

New form:

```text
base distributional/thermal coordinates
x spatial/noise scale
x angular irrep
x deck/Smith/homology
x combinatorial map/charge
x intervention
x insertion/Q-tangent grade.
```

The point is not to fit the full tensor product; it is to identify the smallest grading required by held-out prediction.

### B. P491 cut network x terminal algebra — **retained as theorem route**

The exact cut-network theorem and the mature terminal algebra should still be collided. A natural goal remains a canonical quotient/congruence of the cut network under terminal operations, in the spirit of Myhill--Nerode/bisimulation for the declared continuation language.

This is complementary to #584: #584 asks whether a small fiber explains current marginal flows; the cut-network route asks what information is required for arbitrary compositional futures.

### C. Rao--Blackwell original-U using cut-network conditioning — **retained**

Conditioning an expensive signed original-U estimator on an exact cut-network/terminal state remains a high-value variance-reduction route. It should be treated as estimator design rather than evidence for the sufficiency of that state.

### D. Global multi-observer latent model + proof-carrying elimination — **retained but retyped**

The old proposal

```text
Y_(observer,generator,context) = C_observer^T A_generator B_context
```

is still useful, but the latent columns should now be constrained by #582/#581/#584 rather than learned as unconstrained PCA modes.

#370-style exact/model-class elimination is especially valuable after a candidate graded representation has explicit algebraic constraints. Use it to prove that a whole model class cannot fit, not to certify one numerically convenient basis.

### E. Distributional RG / quantile transport — **promoted; #582 is the realization**

The earlier idea of studying `F_(mN)^(-1) o F_N` has effectively matured into the current Wasserstein/quantile tangent analysis. The next step is state typing (#584), not another family of scalar distribution moments.

### F. Optimal observable synthesis — **retained, later**

Once two surviving mechanism subspaces are explicit, construct a readout that maximizes covariance-weighted model separation subject to exact symmetry/semantic constraints. Do not optimize before the candidate spaces are typed; otherwise the readout can learn the post-reveal residual.

### G. P250 spatial spectrum x Gaussian/deck transfer — **retained as fiber diagnostic**

The high spatial spectrum should be compared with the smaller endpoint/full-law state, not treated as a contradiction. Ask whether the extra spatial modes load the same deck/map fiber that #584/#249 identify.

### H. Complex zeros — **retain only as finite fingerprints**

Use zeros/Galois/defect polynomials to classify finite analytic structure or distinguish model families. Do not use finite complex-zero patterns as a shortcut to continuum exponents unless a theorem supplies the bridge.

### I. Branching reliability signatures — **retained**

The exact growing predictive-class results make branching interventions particularly valuable. They are better tests of state sufficiency than more unbranched endpoint samples.

### J. Finite transfer operator rather than scalar Jordan fits — **retained**

But now require either a minimal realization (#249), composable annulus state (#253), or intervention-stable representation (#580). Do not fit a 2x2 matrix merely because two scalar channels are available.

### K. Machine-assisted theorem discovery on exact catalogs — **retained as support**

Use exact terminal/cut/topology catalogs to conjecture congruences, monotonicities or obstruction identities, followed by symbolic/exact proof. Avoid turning catalog complexity into a surrogate physical claim.

### L. Null-field/logarithmic residue (#252) — **retain as separate high-risk higher-point project**

Use the triangular Camia--Feng control first and test fixed cross-ratio residue shapes. Do not use rank-3/log-squared terms as rescue competitors for the same first-order Matching-One sequence.

---

## 12. Recommended execution sequence by information gain

This is a scientific sequence, not a permission hierarchy.

### Immediate / existing data or exact work

1. **Close the missing N=580 covariance entry deterministically** so the last covariance-sensitive frozen ray has a definitive projective score.
2. **#584:** identify or falsify a discrete/context index for the #582 remainder using existing data and permutation controls.
3. **#581:** move from the exact gate to one square-bond multiscale/noise control, keeping `B_even` and `X` separate.
4. **#249:** compare current realization with covariance-weighted fitting and report distinct rank notions.
5. **#586 Phase A:** derive the exact CLE boundary tangent and inhomogeneous ODE before any Matching-One use.

### First new stochastic block

**#583 N=650** remains the best immediate new production if compute is spent, but its claim should be narrow:

```text
measure A8 jointly with A4;
do not call a 3-angle interpolation a validation of the angular form;
pre-register one overdetermining angular holdout.
```

### Theory in parallel

- **#585:** construct map-resolved torus solution spaces and score existing modulus data.
- **#586:** positive-control graded tangent crossing.
- **#250:** only after a specific map/charge field is selected, build the exact normalized OPE observable.
- **defect/tube route:** square-bond positive control before square-site interface claims.

### After the state is frozen

Run **#580** exact P398 intervention first. If it succeeds, feed the #582/#584 state to the first square-site intervention without relearning the basis.

### Reserve until a compositional map exists

- annulus/radial spectroscopy beyond current static rows;
- full open–closed/tube realization;
- higher-point rank-3/null-residue program;
- large new spatial/charged campaigns.

---

## 13. What not to do next

The following now have low information relative to the alternatives above:

```text
another ordinary larger-N one-dimensional sequence;
another free exponent fit;
another scalar A+B log N fit;
a new modular function added after the N=580 reveal without a map/representation origin;
calling three fitted angles a test of a three-parameter harmonic model;
identifying Hankel rank with positive-state dimension or CFT field count;
relearning a latent basis independently in every intervention;
using an unmarked invariant first-order observable to infer a forbidden deck charge;
rerunning the old broad charged table instead of constructing a field-faithful OPE insertion;
interpreting high annulus rank as failure before comparing with the exact complex spectrum.
```

---

## 14. Literature read in this pass

Primary/recent anchors most directly relevant to the recommendations:

- Roux, Ribault, Jacobsen, **Torus one-point functions in critical loop models**, arXiv:2604.24491.
- Ang, Cai, Jacobsen, Nivesvivat, Roux, Sun, Wu et al., **Exact solution of three-point functions in critical loop models**, arXiv:2604.05503.
- Jacobsen, Nivesvivat, Ribault, Roux, **Three-point functions in critical loop models**, arXiv:2510.04701.
- Grans-Samuelsson et al., combinatorial-map formulation of critical loop correlations, arXiv:2302.08168.
- Cai, exact boundary CLE four-point/connectivity functions, arXiv:2603.28161.
- Sun, Xu, Zhuang, exact annulus crossing laws / complex subleading spectrum, arXiv:2410.04767.
- Camia, Feng, explicit percolation energy/logarithmic-partner lattice fields, arXiv:2508.16047.
- Liu, Jacobsen, Saleur, emerging Jordan blocks in lattice models, arXiv:2403.19830.
- Choi, Rayhaun, Zheng, generalized tube algebras and symmetry-resolved partition functions, arXiv:2409.02159.
- He et al., weighted least-squares viewpoints for realization, arXiv:2505.19639.
- Taghavian, Sjolund, minimal positive Markov realization, arXiv:2502.21102.
- Radhakrishnan, Tassion, strict inequalities between arm exponents, arXiv:2410.23250.

The literature supports the methodology and supplies positive controls. It does **not** identify the Matching-One observable with any one of these continuum constructions without the lattice semantic maps described above.

---

## Bottom line

The strongest current research hypothesis is no longer “there is one hidden Jordan partner.” It is also not “the state is hopelessly high-dimensional.”

The data support a more precise and falsifiable middle position:

```text
there is a strong transferable low-dimensional base,
but a small structured fiber survives;
that fiber may be topological, quotient, map, charge or intervention context,
and the repository now has exact tools to distinguish those possibilities.
```

The next phase should determine the **type and compositional law of that fiber**. Continuum field naming should come after that finite structure survives angular holdout, map-resolved modulus prediction, positive-control tangent/OPE checks, and intervention transport.