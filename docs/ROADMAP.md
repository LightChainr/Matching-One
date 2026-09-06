# Roadmap

This roadmap optimizes for **information gained per unit effort**. It is not a permission system. Existing-data analysis, exact work, pilots and exploratory production may proceed whenever useful.

While exploring, the whole rule set is `GOVERNANCE.md` §2 — five items, of which three bear on this roadmap: chronology, observable-semantic compatibility for claim-bearing scores, and non-duplication of correlated evidence. Everything else applies when there is a paper.

## Theory-bound questions

Three items below are blocked on theory rather than on compute, and no amount of
sampling moves them. They are packaged as self-contained questions for an external
mathematical model in [`docs/astra/`](astra/README.md): the descendant logarithmic
coupling that #275 needs to reopen, the additive-shape ambiguity that decides whether
item 2 is worth running, and an escape from the unit-automorphism no-go that closed
the Gaussian-cover line. An answer there is a theory input, not evidence; the pack's
README states the ordering that keeps it out of the claim ledger.

## Active — highest information now

### 1. Same-N norm-5 coalescence control — #205

The N325/N425 coalescence design is already frozen and cheap relative to large production:

```text
N325: 5 M_C - 11 M_A + 6 M_B = 0
N425: 20 M_C + 13 M_A - 33 M_B = 0
```

It removes radial exponent, parent amplitude and thermal metric. The new C nodes also change Smith class, so one block tests H4 interpolation, conjugation and quotient sensitivity at once.

**Action:** run/score the fixed 10M common-field block if the C targets are still unrevealed. If it fails specifically at C, move quotient/deck structure ahead of extra RG fields. If it passes, norm-4 becomes a cleaner radial/Jordan test.

### 2. Thermal Q4/Jordan modulus fingerprint

Scale-log behavior is no longer enough to identify the Q4 Jordan module. Use exact shape assets instead:

- `B_logN/A_epsilon = -493/192 * g2(tau)` in the frozen module normalization;
- rectangular/CM `11/4` relation;
- hexagonal degree-2 E4 child phase projector `(1,zeta,zeta^2)` which cancels a common scalar mode.

**Run, and negative** — see `notes/modulus-fingerprint-n290-result-20260905.md`. The N=290 square-vs-rectangular ratio is `1.880 +/- 0.177`, which excludes `11/4` at 4.9 sigma along with every other prediction on the frozen list. The hypothesis tested was a conjunction — weight-4 shape *and* a normalization that removes the same block — so what died is the cheap version of the fingerprint, not the module. The known spin-8 systematic would need `A8/A4 = 3.44` to explain the gap and is ruled out by the committed H4-beats-H8 results.

**Designed and frozen, not yet scored — `predictions/aspect_ladder_n580_20260905.yaml`.** The design is the *ladder* `r = 1, 2, 4` at one site count, and it lives at **N=580** — not the N=3380 this item first named, nor the N=1300 of the two superseded prediction files.

Each rung is one paired run: two orientations of the same norm, giving `(O₁−O₂)/Δcos4 = A4 + A8·(Δcos8/Δcos4)`. Three runs in all. N=580 is picked by a stated objective over the sixteen site counts up to 1000 that carry the ladder — **maximum shared angular leverage first, then minimum spin-8 leakage** — and it wins on both: `Δcos4 = 8064/4205` in all three rungs, the same maximum the N=290 design had, and leakage `1148/21025 ≈ 0.055` against `196/625 ≈ 0.31` for the runners-up. The search is in `results/aspect-ladder-design/latest.json`, not asserted here.

The property that makes it better than a bigger, more expensive design: **the r=1 and r=4 rungs carry the *same* leakage**, so the spin-8 bias cancels to leading order in `A4(4i)/A4(i)` — the ratio that discriminates. The N=290 pair could not do this; its two families had equal and *opposite* leakage, so the systematic entered the score twice. It still does not cancel in the `r=2/r=1` entry, and the frozen file says so.

`r=4` is where the live hypotheses separate: weight 4 predicts `10.99`, the bare aspect ratio `4.00`, area `16`, none `1`. At `r=2` the first two are `2.75` and `2.00`, which the existing 9 % measurement cannot split — which is why the ladder goes to 4. The `r=2` rung is separately a replication of the N=290 number.

The linear-in-`r` law is named in the frozen file as `bare_aspect_ratio` **before** any block runs, which is the whole point: it is a post-hoc reading of the N=290 point and this is its one chance to lose.

**Why not N=1300, measured rather than guessed.** A 1M pilot there returned a per-difference noise of `0.0131` against `0.0065` at N=290 — a factor 2.0 for a 4.5× larger torus — while the amplitude falls roughly as `N^-5/4`. A decisive ratio at N=1300 is about three orders of magnitude beyond what we can spend. The same pilot found that **the analysis path returned zeros at N=1300**: the binomial tail's recurrence starts at `(1−p)^N`, which underflows to exactly zero near **790 sites** at the percolation threshold and then stays zero, silently. That bound had never been noticed, and it capped every future large-`N` plan in this repository. It is fixed (`analyze_p48_retrospective` now anchors the recurrence at the mode) and `N ≳ 4000` is analyzable.

**Piloted 2026-09-05** — `notes/aspect-ladder-n580-pilot-20260905.md`. Measured throughput 33 s/M samples and per-difference noise 0.0018–0.0024 per 10M, comparable across all three rungs. **200 M samples per rung, three rungs, about 5.5 hours** puts ~20 % on the denominator of both score entries, which separates 4.00 from 10.99 at r=4 with room to spare. The pilot's central values are noise (33 %, 1416 %, 22 % relative) and may not be read or pooled.

**Action:** run the three rungs — ticket #567.

**Run, 2026-09-05 — underpowered** — `notes/aspect-ladder-n580-result-20260905.md`. `A4(4i)/A4(i) = 4.58`, scored by Fieller contrast because the denominator `A4(i)` is only 3.6σ from zero (the ratio z is recorded but not used). `no_modulus_dependence` `1.00` is excluded at `9.5` sigma; the bare aspect ratio `4.00` (z=+0.50), the weight-4 shape `10.99` (z=−2.08) and area scaling `16` (z=−2.56) all survive — underpowered, three survivors. The 3σ Fieller interval is `[2.40, 27.47]`. The cross-rung covariance was measured by the #575 deterministic replay (`ρ=−0.1648`, reconstructed was −0.1526) and moves no verdict. The r=2 entry (`3.23 +/- 0.93`) cannot split `2.75` from `2.00`, which is why the ladder went to r=4. No optional stopping: any further run is a new frozen design.

Aspect ratio **3 is arithmetically impossible** here: `N = 3|w|^2` and 3 is inert in `Z[i]`, so no 3:1 rectangle shares a site count with a square torus. The reachable ladder is `r` that are themselves sums of two squares.

### 3. Norm-4 dyadic closure with deck characters — #154

The N260/N340 pilot and production scorer are ready. The existing allocation reaches about three-sigma expected q2/Jordan separation at the tested cost.

Exact quotient arithmetic adds information rather than a gate: norm-4 `2i` has `Z2 x Z2` deck group and `(1+i)^2=2i` has an exact coarse/detail Hadamard split.

**Action:** production may run whenever compute is available. If cheap, record character-resolved sufficient statistics in the same run; do not postpone production merely to perfect that extension.

### 4. New local pivotal/RG readout on an injective geometry — #155

The microscopic second direction exists exactly, but N130/N170 response matrices remain nearly rank-one. The multiradius prototype also shows `R=8` is non-injective on those tori and the observed shells do not support a simple constant log-flow.

**Action:** stop adding replicas to the same N130/N170 rows. Choose the smallest larger geometries where `R=2,4,8` are injective, or introduce a genuinely different local/sublattice perturbation. Freeze the observable first, then run a modest covariance pilot.

## Ready — useful parallel work

### Primitive square-bond spin-4 sector — #156

Two prospective norm-2 generations already establish repeated negative H4 phase transfer while positive-phase adversaries fail. Vacuum-KdV gives an excellent zero-new-compute geometry ratio.

**Do not run a third same-purpose norm-2 generation just to show another sign flip.** The next useful target must distinguish finite-size corrections or the KdV/identity-family shape, for example an amplitude-free modulus ratio or a new character projector.

### Multi-u / intrinsic coordinate — #119

The N145->290 quantile-center `N^-3/4` transfer passed while the width metric drifts precisely. Multi-u work is useful if it separates coordinate nonlinearity from S-prime/Jordan dynamics; it is not a new independent evidence block when built from the same histograms.

### Boolean/noise and energy-log-pair exact programs — #227/#234

The exact/no-new-compute programs in open PRs #245/#246 can proceed in parallel. Treat them as mechanism-discovery tools. They do not block the active compute choices above.

### Publication track P2 — the algebraic exclusion manuscript

`docs/manuscripts/p2-algebraic-exclusion/` is drafted end to end, with every number
generated from committed artifacts. It needs no new compute.

The blocker is cleared. [#574](https://github.com/LightChainr/Matching-One/issues/574)
returned the primary reading: Ziff, Phys. Rev. E **73**, 016134 (2006) prints the "A
lattice" quintic and its own statement that the method does not reach square site, both
now quoted verbatim in §1.1; the Suding–Ziff sentence is quoted as printed rather than
as the shortened paraphrase we had. That reading also corrected §4.2 — the historical
exact-bond record reaches height **6** (Wierman 1984), not 4, so the class the paper
exhausts is a choice we defend rather than a bound the literature hands us, and §4.2 and
§8.1 now say so.

**One `[LIT]` marker remains:** Scullard 2006 (martini lattices), which #574 did not
cover. It supports a background sentence, not a result.

**Action:** read Scullard 2006, or weaken the sentence that cites it; then submittable.

### Publication track P3 — denominator-free projective inference

`docs/manuscripts/p3-projective-inference/` is drafted end to end from committed
artifacts and needs no new compute. Target: *Physical Review E*.

The methodological claim is that a model predicting proportions predicts a **ray**, so
the test is the covariance-weighted distance to that ray and no coordinate should be
nominated as a denominator; for two entries and one ray this is exactly Fieller's `z`
squared, verified to `1.5e-15` on the real N=580 covariance, so every verdict change is
attributable to the third rung and not to a change of statistic. Using all three rungs
flips two of eight frozen verdicts from compatible to excluded, at 7.0σ and 7.1σ.

The physical finding is that no competitor's class can produce the measured concavity:
`f[1,2,4] = -4.66e-04 ± 1.53e-04`, `z = -3.05`, against zero or strictly positive for
every frozen ray. Reconciling the rung the design dropped needs `|A8/A4|` between 3.2 and
17.5 against the `<< 1` under which it was dropped -- solved exactly from the per-rung
leakage signs, not bounded.

**Action:** the draft's own §6 —
[#583](https://github.com/LightChainr/Matching-One/issues/583), N=650 with three
orientations per family, which measures `A8` instead of assuming it. One deterministic
N=580 replay would also settle the single undetermined verdict of §5.2.

**#583 is now blocked on an identifiability problem, and the cheap fix has been tried
and failed.** [#589](https://github.com/LightChainr/Matching-One/issues/589) showed the
N=650 design moves the lattice angle and the Smith/quotient class together, so its third
fitted coefficient absorbs any noncyclic offset;
[#591](https://github.com/LightChainr/Matching-One/issues/591) made the leakage exact
(`0.32494 delta` into `A8`, `±0.77170 delta` into `A4`, same-sign `A8` in both modulus
families). #589's option 2 was to calibrate that offset from #205's archived histograms
in the derivative channel — no new sampling.

That attempt is complete and negative (`notes/p205-derivative-smith-loading-20260906.md`,
`results/p205-derivative-smith-loading/latest.json`). Rescoring the frozen #205 block in
`M`/`S`/`D`/`Sp`/`Dp` at the same `p_ref` with the same aligned jackknife — the `M`
channel reproduces the published `analysis/score.json` exactly, as the control that this
is the same pipeline — gives a noncyclic offset consistent with zero everywhere
(`max |delta_z| = 1.83`) but an angular amplitude that is *also* barely resolved
(`|A4_z|` from 0.06 to 2.18). The exported ratio `rho = delta / A4` is therefore
unconstrained: in **all ten** channel×size cells the 2σ interval for `rho` reaches the
pole `rho = ∓1.29584` at which the quotient response annihilates N=650's fitted `A4`
entirely. This is a power failure, not a detection, and it holds in the published `M`
channel too — #205 bounded an offset, never an offset *relative to* the amplitude.

Cost of fixing it by this route, assuming a zero true offset: `Sp` is the only viable
channel (~3× the existing 10M/node for a banded reading, ~26× for a clean one); `Dp`
needs 2.6e4–1.5e6×. So the live options are #589's 3 and 4 — the all-cyclic N=2210
design, or a crossed design — and one cheap prospective question that changes the cost by
an order of magnitude: **which channel maximizes `A4` resolution at fixed sample count?**
`Sp` beats `M` by 4.5× on exactly that figure of merit and nothing has yet been optimized
for it.

### State-object contraction — #582 / #581 / #580

Three zero-new-sampling tests of what the finite "state" actually is, ahead of any
further large production.

**#580 is complete** (`notes/p398-intervention-transport-20260906.md`). On P398's exact
noncrossing process a frozen rank-6 observable-Krylov span, with its readout and source
coordinates and its reduced tangent `B` all fixed at `eta = 0`, predicts the intervened
response at `eta = ±1/4` at a cost of `+0.006` relative error over its own `eta = 0`
truncation, stably from width 4 (14 states) to width 8 (1430). It survives an
intervention deliberately chosen **outside** the operator pencil that built it, so the
pass is not an algebraic tautology. A random span of the same rank sits at 0.997.

The failure is elsewhere and it is the useful part: the same frozen realization does not
represent the five readouts it was not built from — `covering_depth` 0.73,
`halves_linked` 0.52 — and that gap is already there at `eta = 0`. Additive local counts
transport; nesting and long-range boundary linking are not represented at all. The
low-dimensional state is real and transportable **relative to a declared observable
dictionary**, and is not a state of the system.

**#588 Phase A is complete** (`notes/p398-projected-memory-20260906.md`). The exact
Mori-Zwanzig kernel of that same frozen span answers the fork #580 opened. On the
declared dictionary the residual is 96% projected memory, of effective order 3-4,
decaying in `tau ~ 0.14` — a small state plus three or four poles is a real candidate
description. On the held-out readouts memory and unrepresented-readout are jointly
required, ~40% and ~43% of the closure error each, stably from 14 states to 1430.
Neither repair alone is the object.

Two qualifications carried in the note. The memory order is bounded only relative to
an accuracy target: four poles reproduce 99.9% of the kernel at every width, but the
exact numerical degree goes 4, 9, 12, 13, 14 — sublinear and apparently flattening, not
flat. And the small leading order is partly arithmetic: `rank C = 3` uniformly, because
a rank-6 Krylov prefix over a 3-dimensional seed span contains its own first level, so
`rank K(tau) <= 3` for structural reasons. What is not an artifact is the pencil
comparison: per unit of perturbation the out-of-pencil intervention moves
the kernel 1.80 against the declared one's 2.58, where #580's exact lumping went from
750 blocks to the identity partition under the same probe. The three descriptions of
one object -- `r_transport`, `r_positive`, memory -- respond differently to the same
intervention, which is the concrete reason a single "state dimension" cannot be
reported for it.

**#588 Phases B and C are complete** (`notes/p398-memory-closure-20260906.md`), and
they revise the headline above. #580's held-out failure was a property of the *span*,
not of the dictionary. A finite-horizon balanced realization of **order 4**, built from
baseline data only and frozen at `eta = 0`, predicts all eight readouts' contrast
responses at `eta = ±1/4` to 5.9% at width 8 — where the rank-6 Krylov span sits at
41.8% on the held-out block and never passes at any rank up to 12, and where rank 6
plus the *exact* memory kernel reaches only 0.245. A four-dimensional Markov state with
no memory at all beats both, on more readouts, with fewer coordinates.

The Hankel spectrum barely moves with width (D0: `1, 0.45, 0.09` at width 4 and
`1, 0.58, 0.14` at width 8), so the task exposes 3-4 Markovian input/output coordinates
independent of state count. That compact object is a state of the triple
`(sources, readouts, horizon)` — change any of the three and it changes — but it is
genuinely compact and it transports.

Two secondary findings. The Krylov prefix ladder is **not monotone**: rank 7 is worse
than rank 6 on the held-out block at widths 5-8, because a rank cutting inside a Krylov
level adds half a level that the Galerkin closure then misuses. And the memory closure
saturates at **three poles** — real, exact, and dominated.

**Action:** this is the input #275 needs. Score its mechanism classes on dictionaries
that are known to separate, not on a shared one; two models sharing a dictionary sharing
an image is the expected outcome and reads `UNIDENTIFIABLE_WITH_CURRENT_ASSETS`. It also
gives #581 an independent bulk-versus-topology datum from an exactly solvable process.

Not opened by this: any square-site application, which needs a declared intervention with
common observables before and after and a supplied map between microscopic state spaces.

**#598 is complete, and it revises #580's headline** (`notes/p398-reflection-parity-20260906.md`).
The exact positive quotient was never a fragile emergent state. It **is** the orbit
partition of one reflection — checked as the same partition, not merely the same block
count, by exact refinement at every width 4-8 (`10, 26, 76, 232, 750` = R-orbits exactly).
The reflection is identified by the task, not assumed: seven of the eight declared
readouts are fully dihedral-invariant and `wrap` fixes the remaining freedom, so
`R: i -> (w-1) - i` is the unique reflection preserving `D0`. The orbit-count formula
`[Catalan(w) + C(w,floor(w/2))]/2` matches all five widths but identifies nothing on its
own — every reflection has the same fixed-point count, and that caveat is part of the
result.

The collapse under a localized tilt then needs no explanation, but the *survival of the
task* does, and parity supplies it exactly. `H_single` splits into `H_even + H_odd` with
`||H_odd||/||H_single|| = 0.56` — the odd part is over half the perturbation — yet for
`R`-even sources and readouts the first-order integrand
`mu^T e^{(t-s)G} H_odd e^{sG} f` vanishes **pointwise**, to `6.6e-15` over 288 declared
pairs at width 4 and `5.5e-16` over 288 at width 8. The consequence checks out on a
declared epsilon ladder: log-log slopes `2.0006 -> 2.0010` for the odd direction against
`0.980 -> 0.983` for the even one, the latter being the control that a first-order effect
is visible at all. At `eta = 1/4`, width 8, the odd direction moves an invariant task 26x
less than the even one.

The test is falsifiable and fires: `halves_linked` is not `R`-even at odd widths, so the
rule *requires* a nonzero first-order response there, and one appears at `0.292` (w=5)
and `0.107` (w=7). Demanding the lumping carry all eight readouts at w=5 collapses it to
the identity, as it must.

**Action:** #593 changes character — it is no longer "two more rows" but a falsification
of a closed form declared in advance, `r_positive(9) = 2494` and `r_positive(10) = 8524`.
And #594's Q1 no longer needs the codimension conjecture in its body: the closed condition
is a symmetry and the surviving descriptions are invariant. The general statement is the
`C2` case of a character selection rule and is the same shape as #244's deck rule — an
analogy in representation theory, not an identification of the two systems.

## Completed high-information blocks

- **#50 N145->290 full curve:** complete. Corrected slope/root structure survives; a single three-level multiplier shape does not.
- **#57 norm-5 N325/425:** complete. Frozen H4 beats H12/H8; child block alone remains compatible with zero.
- **#212 independent matching-odd synthesis:** complete. Global zero strongly disfavored; fixed H4 compatible.
- **#155 current N130/N170 tangent gate:** complete negative decision for these readouts; do not buy more identical samples.
- **#156 two primitive norm-2 generations:** complete for the sign/phase question.

## Existing-data work — analyze freely

Useful no/low-new-compute work includes low-rank full-curve transfer, covariance-aware thermal-jet mixing, metric-free ratios, standardized profiles, intrinsic/multi-u coordinates, pivotal normalization and exact deck-character projections.

Do not turn multiple derived views of one raw block into extra evidence votes.

## Low-information loops to stop

Not forbidden, simply poor use of time/compute now:

- more N290 replicas repeating completed scores;
- more N130/N170 replicas with the same two self-matching tangent rows;
- a third primitive norm-2 generation whose only purpose is another sign flip;
- another scalar width/boundary correction fit to P57;
- another free exponent fit before testing shape/modulus information;
- (withdrawn 2026-09-05) *"rectangular tori beyond r=2 separate nothing"* — true only against the area competitor, and the N=290 measurement landed on a third hypothesis where longer rectangles separate a great deal. See `notes/modulus-fingerprint-n290-result-20260905.md`;
- large production that stores only final scalars instead of reusable sufficient statistics;
- treating registry/doc synchronization as a prerequisite for science.

## Decision logic

Choose the next experiment by which ambiguity it can kill:

```text
same-N H4 vs quotient/conjugation ambiguity -> coalescence #205
q2 vs generic log vs Q4-Jordan identity     -> scalar-cancelled modulus shape
scale composition + deck arithmetic          -> norm-4 #154
microscopic second RG direction               -> larger injective local pivotal geometry
primitive x≈4 correction structure            -> new KdV/character shape, not more sign tests
```

A failed discriminator is a successful result if it removes a mechanism class.
