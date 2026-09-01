# External-review synthesis: six distinct September 1 packages

**Status.** This is a read-only synthesis of seven supplied ZIP archives, of
which two are byte-identical.  It records evidence, corrects snapshot drift,
and orders possible follow-up work.  Commands, priority proposals, and
"do/not do" language inside the archives are external opinions, not repository
instructions, production authorization, issue locks, or closure decisions.

## Archive identity and evidence boundary

| supplied archive | SHA256 | treatment |
|---|---|---|
| `matching_one_p275_complete_audit_20260901.zip` | `6d073fe5046d3f2a5b39341e1ea1b13c3cf3fa4fdcdcd801c05ed63be4f3c11d` | distinct review A |
| `matching_one_updated_review_complete.zip` | `ad2963f6b26d4d1611761483c985f86f7a6bd1a474d93cf262bb4aad5049febf` | distinct review B |
| `matching_one_p275_review_20260901.zip` | `61b64fa8f52c1971d28304b20993a0671d7c054eb0c0c06e25ac4118101f0a9f` | distinct review C |
| `matching_one_observable_review_20260901.zip` | `35b109c872565ef4b6539db71248401d778818289b62b23ab3f2ce8bc08b3f6d` | distinct review D |
| `Matching_One_更新研究与验证完整包_20260901.zip` | `12c92dbb7fdc5f1fb6d52736958b86b87c649670876a09010379027878d7cfe5` | distinct review E |
| `matching_one_integrated_audit_20260901.zip` | `ad3b017217d8c09629f39bea6a5bfd135afb97799bce7dd6578fcfc8fd0b9465` | distinct review F |
| `matching_one_updated_review_complete 2.zip` | `ad2963f6b26d4d1611761483c985f86f7a6bd1a474d93cf262bb4aad5049febf` | byte-identical duplicate of review B; not a seventh opinion or evidence block |

The packages are unusually reproducible for external reviews: they contain
source, typed inputs or field extracts, manifests, deterministic results, and
run receipts.  Their package-level manifests were checked during intake.  A
manifest verifies package bytes, not upstream raw Monte Carlo bytes or the
mathematical interpretation.  In particular:

- A and C use the same two N25 exact CSV blobs (`2f881f25...` and
  `c1c3a602...`); agreement is a useful implementation cross-check, not an
  independent data block.
- D reuses the one N112 2M/100-batch common-field archive at `2402a333` and a
  transcribed `9 x 9` covariance.  Its model tests are correlated post-reveal
  analyses of that block.
- E reuses the P156 archives at `9de4b016`; its C/Q reanalysis does not create
  a new prospective generation.
- F transcribes the published N65 four-cell result and covariance.  Its four
  original TSV shards, P45 baseline bytes, and full/selected joint delete-one
  factors were not present.

## Snapshot chronology

All six reviews correctly recognized Issue #275 comment `5490477497`
(`2026-09-01T07:30:37Z`) as the P0 observable/normalizer/identifiability
turn.  They do not all describe the same later repository state.

| review | last explicitly pinned navigation/science state | relation to the current September 1 chain |
|---|---|---|
| A | PR #267 `f6a86759` | after the observable-identifiability routing; before paired-C3 contract/result `9eeb6700`/`0b9e89c9`, forward score `7ab2bb1b`, and H8/H0 interpretation `93e47066` |
| B | PR #267 `e4448552`; #537 comment cutoff `5490511026` | before `f6a86759`, `9eeb6700`, `0b9e89c9`, `7ab2bb1b`, and `93e47066` |
| C | Issue #275 snapshot updated `07:31:21Z`, PR #509 `ac5f5fe7`; no later paired-C3 result in its source ledger | before the paired-C3 and forward-score results above |
| D | PR #267 `e4448552`, N112 `2402a333` | before the current routing and all later paired-C3/forward-score commits |
| E | PR #267 `e4448552`, P156 `9de4b016` | before the current routing and all later paired-C3/forward-score commits |
| F | N65 chain through `f9ba1ff6` and comment cutoff `07:33:31Z` | includes the contact-stage/full-J65 update, but precedes `f6a86759` and the paired primitive-C3/forward-identifiability results |

The current interpretation therefore adds three facts unavailable to these
snapshots:

1. the frozen equal-area primitive-C3 gate at `0b9e89c9` rejects H4
   (`73.6412/1`, `p=9.369e-18`) and retains H8
   (`1.1122/1`, `p=.291603`);
2. `93e47066` shows that a signed-scalar H0 line also survives
   (`1.32069/1`, `p=.250468`), so the result is an H8/even branch against H4,
   not a unique H8 local field;
3. `7ab2bb1b` excludes fixed q2 `kappa=.5` (`p=.0189631`) while a Jordan
   affine image remains compatible (`p=.169092`) but unidentified.  The
   subsequent existing-data audit also excludes the naive primitive-H8 radial
   transplant `kappa=2^(-11/8)` (`p=.00964044`).

## Adoption, correction, and deferral matrix

| review | evidence-bearing contribution | adopt into the scientific map | correction required by current state | retain as lower-attention option |
|---|---|---|---|---|
| A: complete P275 audit | exact N25 `12 -> 2` sector-jet forward map; direct reroot control; general known-complex-gain H4/H8 theorem; Q-normalizer term; N65 arithmetic reconciliation | the full moving-root formula, exact source units, general complex-gain degeneracy condition, and the distinction between formal rank and physical mechanisms | its proposed second calibrated rotation is now completed for the primitive observer; `15 degrees` is optimal only for its stated common-complex-gain class, not for arbitrary signed-real gains | N25 map as a calibration oracle and Q-path rule; do not count its N65 marginal-SE reconstruction as batch covariance |
| B: updated frontier review | positive N25 transverse-jet nonidentifiability certificate; completion allocation changes under common K gauge while full U is invariant; q/E normalizer recovery through two thermal derivatives | the statement that critical traces plus a common slope do not determine U; require a horizontal/gauge representative before interpreting a stage share | its C3 acquisition proposal predates `0b9e89c9`; translate the archive's prohibitions into attention ordering rather than locks | its bounded scorer remains useful if an actual candidate supplies the typed jet/gain contract |
| C: P275 review | independent L3/L4 site-colour oracle; relative normalizer and cross-multiplied q/E null; normalized differential operator; exact N25 forward interface; conditional C3 design | the specified-lift normalizer identities and the requirement that a Ward claim include rank-1/total-trace and transverse jets | its "missing second rotation" statement is superseded for primitive C3; its N65 unit conversion is a support audit, not new P275 evidence | preserve the exact finite model as a semantic control, not a Kreg/generic-Q/continuum result |
| D: observable review | first complete model-elimination score on the existing N112 nine-vector; exact small-bond controls; same-area integer-rotation feasibility audit | the N112 exclusions below and the conclusion that the next model must change the winding-resolved rank-1 distribution, not only its mass or common denominator | these are post-reveal parametrization exclusions, not rejections of all normalizer, all H4, or all H8 mechanisms; N2800 is a feasibility result, not a production requirement | alternate modulus/readout designs are possible after a theory vector exists; the N2800 construction is only one restricted embedding class |
| E: round-8 observable audit | two-coordinate normalizer/numerator control; fixed-source mixed-thermal-jet counterexample; covariance-correct P156 two-generation reconstruction | the mixed-jet requirement and the fact that P156 already supplies one calibrated norm-2 transfer experiment | P156 and the later N65 primitive-C3 gate use different observer coordinates, radial constraints, modulus/size changes, gain nuisances, and dependence structures; they are not contradictory votes | the P156 C/Q extension is exploratory reuse; use it to constrain a typed transport law, not as another independent result |
| F: integrated N65 audit | prospective N65 sign/determinant retention; exact exposure-unit diagnosis; full J65 arithmetic; selected/full coverage boundary; commuting-table counterexample | correct the exposure convention, preserve the frozen signed matrix, and rebuild full/selected/remainder jointly from existing factors if the raw archive is available | the high `p=.95893` split-power fit is post-reveal and nonunique; `theta=-1` is an orthant identity; selected share is an allocation share, not a causal fraction | the horizontal/gauge sidecar is a valuable existing-data refinement, but not a prerequisite that blocks #275 theory or other observer work |

## N112 model elimination: what was actually removed

Review D scores the existing N112 square-bond vector
`(E, Re H, Im H)` across three rho children with its full `9 x 9` covariance.
Each model was proposed after the data and retains its stated nuisance
amplitudes; their p-values are therefore descriptive elimination scores rather
than prospective discovery probabilities.

| declared model class | nuisance dimension | constraint dimension | score | Gaussian-batch Hotelling reference |
|---|---:|---:|---:|---:|
| common denominator only, one real scale per geometry | 3 | 6 | `T2=215.10524` | `p=1.52e-21` |
| rank-1 mass only, one fugacity per geometry | 3 | 6 | `T2=334.36860` | `p=5.16e-28` |
| E free and H restricted to a real rescaling of its reference direction | 6 | 3 | `T2=19.08778` | `p=6.46e-4` |

Two further positive winding-fugacity parametrizations also fail after
profiling the three rank-1 masses: `cos(4 theta)` gives
`T2=143.71930/5`, and `cos(8 theta)` gives `146.57259/5`.  Their whitened
tangent cosine is `.9998168`, so they are both poor explanations and a poor
pair for discriminating mechanisms.

The supported inference is narrow and useful: the observed finite N112 change
contains a rank-1 **internal directional-shape** component.  It is not enough
to change only the common normalizer, only total rank-1 mass, or one common
cosine fugacity.  It does not by itself identify a local field, transfer this
square-bond result to square-site U, or exclude a theory whose Ward/modulus
relations predict a different winding-resolved correction.

## N65 contact-stage exposure unit correction

Review F finds a direct code-unit mismatch:

```text
N25 positive exposure: P
N65 positive exposure: P/N
```

The N65 scorer already restores the fixed-empty-site and Bernoulli importance
weight and then divides exposure by `N`; comparing it directly with N25 adds
an unintended factor 65.  Under a common unweighted-P convention, the
first-birth/double-contact decomposition becomes:

| coordinate | N25 | corrected N65 | two-point finite decay |
|---|---:|---:|---:|
| positive exposure P | `.03269409005` | `.01772483009` | `.64073313` |
| conditional signed intensity K/P | `-8.98416e-5` | `-5.20606e-6` | `2.98083687` |
| signed mass K | `-2.93729e-6` | `-9.22766e-8` | `3.62157001` |

Using P/N on both sizes instead gives `1.64073313 + 1.98083687`; the signed
mass exponent is unchanged.  The repair therefore changes the proposed
allocation of the decay, not the signed four-cell matrix, determinant, frozen
sign gate, or full `J65=-.0016225098893862522`.

The fixed split `3,29/8,3,3` has `Q=.636436/4`, nominal `p=.958930`, but was
chosen after the N65 reveal.  Values from extra power `1/4` through `1` are
also compatible in the same conditional profile.  It is a mechanism
fingerprint to retest, not an exponent certification.  The selected cells
have central signed share `.025505` of full J; absent full/selected
cross-covariance and a gauge-invariant carrier definition, neither the
selected share nor its complement is a physical causal fraction.

## Common normalizer and mixed-thermal-jet contract

All six reviews converge on the same object firewall:

```text
local/source or propagating representation
  -> unnormalized restricted traces Z0,Z1,Z2
  -> geometry-specific normalizer
  -> normalized q and E
  -> geometry projection and pooled moving root
  -> original U or source response.
```

For the explicitly specified integer-Q site-cluster colour lift only,

\[
\frac{Z_\sigma}{Z_{id}}
=\frac{E_{id}-q_{id}}{E_\sigma-q_\sigma},
\]

and the cross-multiplied common-denominator null is

\[
Q(E_\sigma+q_\sigma)(E_{id}-q_{id})
-\operatorname{Fix}(\sigma)(E_\sigma-q_\sigma)(E_{id}+q_{id})=0.
\]

These identities do not automatically transfer to Kreg, a closed-source
lift, generic-Q critical FK, or another cluster fugacity.

At a pooled root, the common forward coordinates are

\[
\delta\widehat Y|_M=jY-RjM,
\qquad
\frac{\delta U}{A_N}
=\frac{jY_t-RjM_t-R_tjM}{M_t}.
\]

Consequently a complete same-source candidate input is not a critical trace
value alone.  It is the three-sector source/thermal jet, equivalently
`{delta Z_r, partial_t delta Z_r}` together with the baseline `Z_r` thermal
jets and the physical normalizer.  The repeated positive counterexamples in
B/C/D/E show that identical critical traces, identical normalized q/E, and
even an identical pooled slope can coexist with different U responses.  This
is the consensus reason to request the mixed thermal jet; it is not a general
denial of Ward or representation constraints that actually determine that
jet.

## P156 versus the new N65 primitive-C3 gate

The two results constrain different statistical models and should coexist in
the map rather than be pooled as H4/H8 votes.

| feature | P156 norm-2 chain (`9de4b016`) | new N65 equal-area gate (`0b9e89c9`, interpretation `93e47066`) |
|---|---|---|
| finite observer | square-bond primitive-homology `C_nontrivial_real`, with Q reflection and S scalar controls | complex primitive real-C3 coordinate after the declared intrinsic baseline subtraction |
| geometry operation | multiplication by `1+i`: rotate by `pi/4` and double area at each generation | two equal-area N65 Gaussian quotients, physical angle `atan2(5,12)` at the square modulus |
| data dependence | six independent size/lineage blocks; middle generations reused in adjacent residuals, with their covariance retained | one paired 2M/100-batch common-random-field block across the two orientations |
| frozen nuisance class | fixed radial transfer, primary H4 ratio `-1/2`; positive-phase ratios `+1/2,+1/4,+1/8` fixed | one profiled signed-real gain per candidate line; common character convention retained |
| result | H4 `7.330306/4`, nominal `p=.119429`; fixed local-H8 `+1/8` gives `461.227/4` | H4 rejected; H8 survives; post-reveal H0 also survives |
| correct boundary | supports an alternating negative rank-4 norm-2 transfer under its fixed scale law | selects an H8/even phase branch for this finite primitive observer, not a unique local H8 primary |

P156 rejects its **fixed positive radial laws**.  It does not exclude an
observer-dependent H8/even phase with a different modulus form factor and a
profiled signed gain.  Conversely, the N65 phase result does not erase the
P156 alternating chain or transfer its conclusion to rho-child E_top,
global K1/K2, or original square-site U.  An observer transport map is needed
before any joint mechanism score.

## Unified future operation order

The reviews and the later commits support the following attention order.
Items can proceed in parallel when their dependencies are independent; the
ordering is not a permission system.

1. **Freeze semantics, not tasks.** Record source, observer, geometry,
   normalizer, nuisance amplitudes, and dependency group for N112, P156,
   N65 contact-stage, paired primitive-C3, rho-child E_top, and global K1/K2.
   Preserve the current branch-only/open-PR status of every result.
2. **Use the N112 elimination as a model requirement.** A candidate for that
   observer should predict the winding-resolved rank-1 probability correction
   `delta pi_l`, hence both `delta E` and complex `delta H`; another common
   denominator or total-mass coordinate has low expected information gain.
3. **Repair the N65 contact interpretation on existing data.** If the archived
   TSV/baseline factors are available, apply one exposure convention and form
   full, selected, horizontal/gauge, and remainder quantities inside every
   shared delete-one replicate.  This refines contribution coverage without
   repeating MC and without delaying the separate #275 theory line.
4. **Complete the restricted-trace transport law.** For a named vacuum/Ward
   and thermal-Q4/Jordan alternative, derive the same-source rank-0/rank-2
   `delta Z_r, partial_t delta Z_r` vector, rank-1/normalizer action, root
   counterterm, and allowed common amplitudes.  Feed it through the common
   two-coordinate forward map rather than fitting an amplitude to every
   observer.
5. **Score existing observer-specific assets before pooling stories.** Use
   K1/K2 for the global matching-odd/topological activation response, P156 for
   its fixed norm-2 primitive transfer, rho-child data for E_top/primitive
   child rays, and the N65 pair for its signed-real primitive phase.  Compare
   constrained model images inside each covariance block; combine only after
   an explicit cross-observer map exists.
6. **Choose the next acquisition by the remaining nuisance.** If the transport
   law separates existing rows, consume them directly.  If primitive H8 versus
   H0 remains the decisive ambiguity, the already proposed convention-
   invariant three-way angle statistic at a geometry where H4, H8, and H0
   predict separated values is more informative than another copy of the
   nearly H8/H0-aliased N65 angle.  If no candidate supplies a distinguishable
   vector, record the missing coupling or jet as the result and keep alternative
   research lines visible at lower attention.

This ordering preserves the external packages' strongest insight--move from
names and nonzero coordinates to original-observable forward predictions--
without importing their snapshot-stale production instructions or turning
scientific priorities into locks.
