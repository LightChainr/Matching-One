# P5 — What a bounded state can and cannot carry: exact continuation, branching, and a no-go

Draft manuscript notes for **Paper 2** of [issue #644](https://github.com/LightChainr/Matching-One/issues/644).
Portfolio track **P5**, following the `p2`/`p3`/`p4` manuscript convention.

| File | Role |
|---|---|
| `README.md` (this file) | the draft: results, boundaries, limitations, omissions |

**Status of this document.** A first-pass writing draft only. It promotes nothing,
computes nothing, and touches `docs/STATUS.md` not at all; levels quoted below are current
`main` content, cited as-is. Issue #644 stays open. No venue selected; the ticket defers
that and this draft does not pretend otherwise. Unlike `p2`/`p3` there is **no generated
`tables.md` and no evidence-assembly script yet** — every value below is carried with the
committed file path (and field, for JSON) where its full-precision value lives. All paths
are relative to the repository root.

## 1. Abstract (draft)

Three exact results about the limits of bounded state descriptions of continuation
processes. (i) On a family of two-terminal cut networks, the complete unbranched survival
law can be *identical* while branching behaviour separates `k+1` exact predictive classes,
and the response-matrix rank achievable by experiments of depth `≤ d` is exactly
`min(d+1, k+1)`: richer branching futures require predictive classes that no fixed depth —
and no bounded-radius terminal summary — supplies. (ii) A constructive no-go: two analytic
families of local, irreducible-at-finite-`L` Markov generators with a *common declared
task* have identical finite-horizon responses to machine precision and closing-gap
thresholds at `1/2` and `1/3` respectively; bounded finite-horizon task rank therefore
carries no threshold information unless critical-sector completeness is separately
assumed. (iii) For a finite connectivity process (P398), the coarsest admissible exact
positive lumping is *provably equal, block for block*, to the orbit partition of the one
reflection its declared readout dictionary preserves — a first-order response selection
rule makes the symmetry-breaking direction invisible to the invariant task (the Duhamel
integrand vanishes pointwise on protected pairs, while the deliberately exposed channel
fires at odd widths) — and the task's balanced realization factors exactly through the
quotient because the odd sector is inert. The through-line: compressibility, symmetry
structure, and threshold content are separate properties of a state description, and
conflating them is the failure mode each result closes. The reflection-fixed counts in
(iii) are classical (Callan–Smiley; Burnside; Ding) and are cited, not claimed; the
quotient identity is verified at five widths and stated as an instance, not a theorem —
coarsest lumping need not equal an automorphism-orbit partition (D'Angeli–Donno,
Prop. 13). P398 is a calibration process, not a percolation model; the result is a
statement about the frozen finite-horizon input/output object.

## 2. In its own voice: what this paper is not

Per the ticket, this belongs in the body, not a footnote:

- **P398 is not percolation.** It is a finite positive width-`w` process program used as a
  calibration and mechanism-de-identification object (`docs/RESEARCH-ATLAS.md` §11). Part
  III is an exact finite symmetry theorem about that declared input/output process. The
  committed claim boundary already says so in the working-note voice —
  `notes/p598-phase-c-balanced-quotient-20260906.md`: "this is a statement about the frozen
  finite-horizon input/output description of one exactly known process, not about the state
  space, and not a percolation result."
- **The cut networks are proof devices in a declared gauge, not physical marks.** The
  occupied-bridge convention is stated as such in
  `notes/no-bounded-radius-quotient-20260906.md` §3 ("the network is a proof device in a
  gauge, not a physical mark"), after #491/#550.
- **No continuum claims.** Nothing here bounds the state dimension of a continuum field
  theory, identifies an operator, or touches `p_c` numerically. Part II is precisely the
  warning against reading threshold content out of state-counting.
- **Levels.** `docs/STATUS.md` line 38 carries Part III's Phase C at level
  "exact (bounded to the frozen I/O object)". Parts I and II have **no `docs/STATUS.md`
  lines** on `main`; they are indexed in `docs/RESEARCH-ATLAS.md` §9 (integration
  `open_pr` for #491/#492) and the #435/#549 rows of its integration table. This draft
  reports that as-is and does not treat writing as promotion.

## 3. Part I — Branching continuations need unbounded predictive classes (#435 / #491 / #549 / #550)

**The fact that starts it.** Two N16 torus configurations with the *identical complete
unbranched survival vector* `(1, 7/8, 9/14, 5/14, 4/35, 0, 0, 0, 0)` — one `L_0`
predictive class — where a single declared fork separates them
(`notes/memory-branching-and-moment-algebra-20260906.md`, FACT #435). The exact fork
probability on the parallel-gadget family (#549) is affine,
`F_{k,a} = [343k^3 − 182k^2 + 25k + 4a] / [8k(8k−1)^2]`, strictly increasing in the hidden
class count `a` with consecutive gap `1/[2k(8k−1)^2]` — `1/98` at `k = 1`, reproducing
#435's witness (`notes/probe-predictive-rank-separation-20260906.md`,
`notes/branching-hankel-rank-lower-bound-20260906.md` §4).

**Theorem A (rank upgrade).** Grouped-fork experiments `E_1..E_k` of depth `1..k` give a
`(k+1)×(k+1)` response matrix of exact rank `k+1` for every `k`; machine-verified by exact
rational elimination at `k = 2..8` — committed rank/class table in
`results/cutnetwork-rank-vs-depth/latest.json` (keys `2..8`, each `rank = classes = k+1`);
statement and proof sketch in `notes/branching-hankel-rank-lower-bound-20260906.md`
(Theorems 1–2). Script `scripts/cutnetwork/rank_lower_bound.py` (note: the working note's
Files section names `scripts/rank_lower_bound.py` / `results/rank-vs-depth/` — the
committed locations differ and are the ones cited here).

**Theorem B (depth is the resource).** The span of all experiments of depth `≤ d` has
dimension `r_d(k) = d+1` for `d ≤ k`: the class count grows like `k+1` while every fixed
depth sees a constant-size response space. The unbranched language (`d = 0`) sees
essentially nothing; #549's single-fork language is exactly `r_1 = 2`.

**Theorem C (no bounded radius suffices).** For every radius `r ≥ 0` and every `k`, `k+1`
hidden classes admit planar two-terminal realizations sharing the unbranched survival law
and the radius-`r` terminal neighbourhood, with fork gap `1/[2k(8k−1)²]` — the bank is
pushed beyond the radius on an occupied bridge. Hence no bounded-radius quotient of the
cut network preserves all branching futures.
`notes/no-bounded-radius-quotient-20260906.md`; script
`scripts/cutnetwork/radius_insufficiency.py`; data
`results/cutnetwork-radius-insufficiency/latest.json` (radii `1,2,4,8`). This upgrades
#550's `r = 1` certificate at every radius simultaneously, and already for depth 1.

**Theorem D (ordinary rank vs nonnegative rank).** An explicit family (slack matrices of
the regular `n`-gon) with ordinary rank `≤ 3` and nonnegative rank `n` (4 at `n = 4`) —
so even positivity, the natural constraint on probability responses, does not rescue a
bounded signed summary. `notes/positive-vs-signed-slack-separation-20260906.md`; script
`scripts/cutnetwork/slack_family.py`; data
`results/cutnetwork-positive-vs-signed-slack-family/latest.json`.

**Background the paper must cite but not claim:** #491's reduction of embedded rank-one
continuation to two-terminal vertex connectivity, with the repository's own boundary —
"an exact finite continuation representation, not a proof of bounded continuum state
dimension" — quoted from `docs/RESEARCH-ATLAS.md` §9.

## 4. Part II — Bounded task rank carries no threshold content (#615)

**Setting.** Declared task `(G_L(p), B_L, C_L)`, finite-horizon response
`R_L(t,p) = C_L e^{tG_L(p)} B_L`; exact task order = Kalman controllable-observable
dimension; closing-gap threshold defined analytically (gap positive and analytic at every
finite `L`, vanishing in the `L → ∞` limit) — all in
`notes/bounded-task-rank-threshold-no-go-20260907.md` §0.

**Theorem E (no-go, constructive).** Two analytic Markov families
`G_L^{(j)} = G_vis ⊕ H_L(p − p_c^{(j)})`, `j = 1,2`, share source/readout: visible sector
a fixed 2-state chain; hidden sector a biased nearest-neighbour walk on `{1..L}` with
reflecting boundaries, rates `a = e^{p−p_c}`, `b = e^{−(p−p_c)}`, exact gap
`γ_L = a + b − 2√(ab)·cos(π/L) → 4·sinh²((p−p_c)/2)`. The hidden chain is connected,
positive, local, irreducible at every finite `L`. The responses are *identical* for all
`L, p, t`, while `p_c^{(1)} = 1/2 ≠ 1/3 = p_c^{(2)}`; both have `sup_L r_ex = 2`.

Committed witness: `results/theory-no-go/latest.json`
(`witness.family1_pc = 0.5`, `witness.family2_pc = 0.3333…`, per-`p` response vectors
byte-identical between families and equal to the visible-only response); generator
`scripts/theory/no_go_theorem.py`. (The note's Files section lists
`scripts/no_go_theorem.py` / `results/no-go/latest.json`; the committed paths above are
what this draft cites.)

**Why the counterexample is not a cheat (and where it is genuinely limited).**
Proposition 1 (hidden ⟺ off-Kalman-CO) is the certificate that the *missing hypothesis*
is completeness, not block-diagonality. The paper then states the three-level honesty
ladder from the note §3: (a) the direct sum; (b) symmetry-protected hiding, which is
necessarily reducible (an irreducible generator admits only trivial commuting finite-group
actions with >1 isotypic type) — P398's reflection-odd sector is the worked instance;
(c) irreducible exact masking, which *is* fine-tuned: committed numbers in
`results/theory-irreducible-masking/latest.json` (`partial_fraction_coeff`: generic
`0.48371173070873924`; readout-masked `4.1291408707408205e-17`; source-masked
`−2.312360816083912e-17`; both `−1.9739160635202627e-33`; script
`scripts/theory/irreducible_masking.py`). What survives generically in the irreducible
class is *weak coupling*, not an exact dark sector — the note says so and the paper keeps
that qualifier.

**Corollary / scope sentence (paper-ready, from the note §5):**

> "A bounded finite-horizon balanced order of a declared task certifies the
> compressibility of that task's input–output map and the existence of an exact or
> symmetry-protected uncontrollable/unobservable sector; it carries no information about
> any thermodynamic threshold unless the critical modes are separately shown to be
> controllable and observable."

**Escape conditions (note §4), stated as the theorem's other half:** critical-sector
completeness is *necessary* and not sufficient; a horizon growing with the correlation
length amplifies but cannot create coverage (exact order is horizon-independent by
analyticity — Lemma 0); a signed-observable zero with a monotonicity argument (the
repository's own `M_L = P_2 − P_0` route) is an *independent* mechanism that never goes
through rank at all. Task compressibility ≠ threshold identifiability: statements about
`(G,B,C)` vs about `spec(G_L)` plus a sign-change argument.

## 5. Part III — The exact reflection quotient (#598 / #600 / #601)

**The object.** P398's positive lumping at widths `w = 4..8` on noncrossing partitions.
The reflection named by the *task*: `R = (i ↦ (w−1) − i)` — of the eight declared readouts,
seven are fully dihedral-invariant and `wrap = [state[0] == state[w−1]]` is preserved
exactly by `R`, the unique reflection preserving the declared `D0` dictionary at every
width (`notes/p398-reflection-parity-20260906.md`).

**Theorem F (quotient identification, at five widths — an instance, §7.2).** By exact
partition refinement, the coarsest `D0`-admissible strong lumping (stable for the declared
`J`, `D` moves) **equals the `R`-orbit partition, block for block** — not merely in
cardinality — at `w = 4,5,6,7,8` with `10, 26, 76, 232, 750` blocks; and the test bites:
at odd widths `halves_linked` is not `R`-invariant, so demanding all eight readouts
collapses the coarsest admissible partition to the identity at `w = 5` (all 42 blocks).
Reproduce: `scripts/p398_reflection_parity.py` →
`results/p398-reflection-parity/latest.json` (fields `by_width[*]`, coarsest-vs-orbit
equality; 17 tests in `tests/test_p398_reflection_parity.py`).

**The counts are classical.** `[Catalan(w) + C(w,⌊w/2⌋)]/2 = 10,26,76,232,750` is Burnside
for `{1,R}` once the reflection-fixed count `C(w,⌊w/2⌋)` is known; that count is
Callan–Smiley (arXiv:math/0510447, Thm 1, 2005) with Ding 2016 (Miami PhD, Thm 2.1.2,
Kreweras anti-isomorphism of the two even-`n` reflection classes) for the even-`n`
second class; equivalently OEIS A007123(`w+1`). **The paper cites; it claims nothing here.**
Full bibliographic detail, verbatim quotes and the explicit do-not-cite list (Reiner–
Stanton–White cyclic sieving is the *rotation* count; type-B NC counts a different object)
are in `notes/literature-officer-20260906-issue601.md` §Q1. The out-of-sample predictions
`r_positive(9) = 2494`, `r_positive(10) = 8524` are stated in the notes (as
A007123(10)/(11)); **widths 9–8…10 have not been run in this repository** — §8.2.

**Theorem G (parity selection rule — pointwise).** Decompose the localized tilt into
`R`-even and `R`-odd halves; over half the breaking direction is odd
(`||H_odd||/||H_single|| = 0.55277, 0.55902, 0.56026, …` for `w = 4..8`; committed at
`results/p398-reflection-parity/latest.json`, fields `by_width[*].parity_split.odd_share_of_norm`).
For an invariant baseline and `R`-even source/readout, the first-order Duhamel integrand
`μᵀ e^{(t−s)G} H_odd e^{sG} f` vanishes **pointwise**, not merely in integral — checked on
a declared grid for all four sources against every readout: protected-pair worst relative
integrand `~6.6e-15` or below; **positive control:** `halves_linked` is odd-containing at
`w = 5, 7` and the rule *requires* it to fire there — observed largest exposed relative
integrand `0.292092871260188` (w=5) and `0.10691283943577032` (w=7), same artifact
(`by_width[*].selection_rule.largest_exposed_relative_integrand`). "A run in which
everything vanished would have been a broken test, not a stronger theorem."
Phase B (second order): the odd direction moves the invariant task at slope `2.0006–2.0010`
in `ε`, the even direction at `0.98`, a factor 26 at `ε = 1/4`, `w = 8`
(`notes/p398-reflection-parity-20260906.md` Phase B).

**Theorem H (factorization).** The task's balanced realization factors through the
quotient: microscopic vs quotient constructions agree with `max |A−B| ≤ 2.7e-15` in every
resolvable Hankel direction at `w = 4..8`, balanced order and numerical rank preserved;
the odd sector (`dim = 4/16/56/197/680`) is **exactly inert**, reach/observe energy
`1.69e-17 … 5.69e-17` — which is why the reduction is exact, not merely accurate. The
declared trap fires: the exposed dictionary breaks at exactly `w = 5, 7`
(`0.017944…` / `0.008650…`) and only there, recovering to `1e-16` under symmetrization.
`notes/p598-phase-c-balanced-quotient-20260906.md`;
`scripts/p598_phase_c_balanced_quotient.py` → `results/p598-balanced-quotient/latest.json`
(fields `blocks[*].dictionaries.*`, `odd_sector.reach_odd_energy_relative`,
`decision.*`); `tests/test_p598_phase_c_balanced_quotient.py`. STATUS level: line 38,
"exact (bounded to the frozen I/O object)" — the parenthetical is the claim ceiling and the
paper keeps it.

**Prior-art position (#601).** Round 1 and round 2 of the literature pass
(`notes/literature-officer-20260906-issue601.md`,
`notes/literature-officer-20260906-issue601-round2.md`) establish: the abstract vanishing
is Schur / Wigner–Eckart; the *statement* for equivariant Markov generators with the
Duhamel integrand vanishing pointwise "is still unnamed" (Hänggi I/II,
Golubitsky–Stewart–Schaeffer, Antown–Dragičević–Froyland 2018, Diaconis 1988 Ch. 3E,
Santos Gutiérrez–Zagli–Carigi 2025 all checked and are adjacent, not it). The paper's
contribution there is **language and the pointwise verification**, claimed as a remark
citing Schur, not a theorem. Phase H's factorization is the Kalman decomposition applied
to a `C2`-equivariant realization — cite Wonham, per round 2's own instruction; the
odd-width exposed control is what makes the test non-vacuous.

## 6. The through-line

1. **What a bounded state carries about *futures*** (Part I): unbranched memory can be
   complete and still blind; branching depth, not width, is the resource; positivity of
   responses does not help (Thm D).
2. **What a bounded task carries about *thresholds*** (Part II): nothing, without
   critical-sector completeness; the witness is analytic, local, irreducible at finite
   size, with the response identical between families to machine precision.
3. **What a bounded state carries about *symmetry*** (Part III): exact quotients *are*
   automorphism orbits when they exist — and here the coarsest lumping provably equals one,
   while first-order response selection rules say precisely which symmetry-breaking
   structure any such description is blind to, and the factorization theorem says what
   survives the blindness intact.

Together: compressibility, branching expressiveness, threshold content and symmetry
structure are four separate properties. Each part is a theorem or construction about one
implication failing to run to the others.

## 7. Boundaries the prose must respect (settled by the ticket; restated here)

1. Reflection-fixed counts: **cite, do not claim** (Callan–Smiley Thm 1 + Burnside; Ding
   2016 even-`n`; OEIS A007123). Cardinality matching identifies nothing — every
   reflection fixes `C(w,⌊w/2⌋)`; the load is carried by the partition equality and the
   task that names `R` (already correct in the working note).
2. "Coarsest lumping = automorphism-orbit partition" is **not a theorem in general**
   (D'Angeli–Donno 2013, Prop. 13 gives a lumping induced by no subgroup; their Thm 12 is
   the version *with hypotheses*, which P398 does not satisfy). The `w = 4..8` agreement
   is **an instance** and the abstract above calls it one. Gate 1's branch-A resolution
   must never be written as "by the theorem relating coarsest lumpings to Aut orbits."
3. P398 is a calibration process, not percolation — §2 in the paper's own voice.
4. Parts I and II: the constructions are about the *declared* protocol (the published
   #549 fork composition, the occupied-bridge gauge, the declared `(G,B,C)` triples); the
   full-N16 lift is explicitly not claimed (`notes/branching-hankel-rank-lower-bound-20260906.md` §5).
5. #636 territory and any `p_c` inference are out of scope entirely.

## 8. Limitations

1. **Five widths.** The quotient identification is proved by exhaustive exact refinement
   at `w = 4..8` only; it is a finite instance, and the generic-rates lemma that would
   turn it into a theorem with hypotheses has *not* been written (round 2 checked the
   neighbourhood — Watanabe–Wolfer 2024/2026 is adjacent and does not contain it).
2. **Widths 9–10 untested.** `r_positive(9) = 2494`, `r_positive(10) = 8524` are
   predictions against A007123; the partition equality re-derivation commits at
   `w = 4..8`; no committed artifact in this repository extends it (or the cardinality) to
   widths 9–10. Reportable as a finding
   rather than produced: per this ticket's rule the gap is stated, not computed.
3. **Part I's theorems are family statements.** They bound what the declared fork/radius
   languages can carry; they do not exhibit a percolation observable whose branching
   futures *require* unbounded classes at any fixed depth (the gap from proof device to
   lattice observable is open and the note's §5 says so).
4. **Part II's strongest irreducible form is fine-tuned or weak-coupling** (§4c): exact
   hiding in the irreducible/local class survives only as single-mode orthogonality; a
   generic irreducible task has only `O(ε)` hiding. The no-go's force is for *declared
   tasks and symmetry-protected sectors*, which is what the paper claims.
5. **Phase C's agreement is bounded to the frozen I/O object** (STATUS line 38
   parenthetical): the `Catalan − r_reflect` "extra" states are not asserted unphysical.
6. **The #601 gap is a negative literature result**, bounded to the searches recorded in
   the two packets; a selection-rule statement could exist under language neither pass
   indexed.
7. **No STATUS promotion.** Parts I–II have no claim-ledger rows; that is the honest state
   of integration (`open_pr` snapshot per atlas), and the draft treats it as such.

## 9. What had to be left out (for the next draft)

1. `references.bib` entries for this cluster — Callan–Smiley 2005, Ding 2016,
   D'Angeli–Donno 2013, Kemeny–Snell 1960, Buchholz 1994, Godsil–Royle 2001,
   Hänggi 1978 I/II, Wigner–Eckart (textbook cite), Diaconis 1988, Antown–Dragičević–Froyland
   2018, Wonham, OEIS A007123 — **are not in `references.bib`** (10 keys, none matching;
   a finding). Deliberately not added here: the full bibliographic detail lives in the two
   #601 packets and the bib patch should be their transcription with `PRIMARY_TEXT_READ`
   marking per convention, done in a venue-prep pass, not silently inside a no-compute
   writing PR.
2. Generated tables / evidence-assembly script (`p2`/`p3` style) locking every Part III
   number to its JSON field with a drift test.
3. Widths 9–10 exact refinement run (new computation — separate ticket).
4. The generic-rates lemma (§8.1) — optional theory that would upgrade instance→theorem.
5. Related-work prose on lumpability/bisimulation (round 1 §Q2 "Bisimulation / coalgebra")
   and the model-reduction Krylov-prefix material (round 1 §Q5) — adjacent to Part III/II
   framing but beyond this paper's claim budget.
6. Venue. Ticket says the owner picks once drafts exist; this draft neither proposes nor
   pretends one was chosen.

## 10. Evidence index

| # | statement / value | committed artifact (field) |
|---|---|---|
| 1 | #435 identical unbranched vector; single `L_0` class | `notes/memory-branching-and-moment-algebra-20260906.md` (Program H, FACT) |
| 2 | `F_{k,a}` closed form; gap `1/[2k(8k−1)²]`; `1/98` | `notes/probe-predictive-rank-separation-20260906.md`; `results/cutnetwork-radius-insufficiency/latest.json` |
| 3 | rank `= classes = k+1` at `k = 2..8` | `results/cutnetwork-rank-vs-depth/latest.json`; `scripts/cutnetwork/rank_lower_bound.py` |
| 4 | no bounded-radius quotient (all `r`) | `notes/no-bounded-radius-quotient-20260906.md` §2; `scripts/cutnetwork/radius_insufficiency.py` |
| 5 | ordinary rank ≤ 3 / nonnegative rank `n` | `notes/positive-vs-signed-slack-separation-20260906.md`; `results/cutnetwork-positive-vs-signed-slack-family/latest.json` |
| 6 | #491 cut-network theorem + boundary | `docs/RESEARCH-ATLAS.md` §9 |
| 7 | Theorem E families; `1/2`, `1/3`; identical responses | `notes/bounded-task-rank-threshold-no-go-20260907.md` §2; `results/theory-no-go/latest.json` (`witness.*`); `scripts/theory/no_go_theorem.py` |
| 8 | masking coefficients `0.4837…/4.13e-17/−2.31e-17/−1.97e-33` | `results/theory-irreducible-masking/latest.json` (`partial_fraction_coeff`); `scripts/theory/irreducible_masking.py` |
| 9 | scope sentence (paper-ready) | note §5 quoted verbatim above |
| 10 | `R` unique to dictionary; block-for-block equality; `10,26,76,232,750`; `w=5` identity collapse | `notes/p398-reflection-parity-20260906.md`; `results/p398-reflection-parity/latest.json`; `scripts/p398_reflection_parity.py` |
| 11 | classical count citation + do-not-cite list; A007123; 2494/8524 | `notes/literature-officer-20260906-issue601.md` §Q1 |
| 12 | odd share `0.55277…`; exposed `0.29209…/0.10691…`; protected `~1e-15` | `results/p398-reflection-parity/latest.json` (`parity_split.odd_share_of_norm`, `selection_rule.largest_exposed_relative_integrand`) |
| 13 | Phase B slopes `2.0006–2.0010` vs `0.98`; factor 26 | note Phase B table; same artifact |
| 14 | Phase C `≤ 2.7e-15`; inert odd `1.7e-17…5.7e-17`, dims `4/16/56/197/680`; breaks `1.8e-2/8.7e-3` | `results/p598-balanced-quotient/latest.json`; `scripts/p598_phase_c_balanced_quotient.py`; `tests/test_p598_phase_c_balanced_quotient.py` |
| 15 | STATUS line 38 level | `docs/STATUS.md:38` (read-only) |
| 16 | #601 gap statements (Q3 round 1/2) | `notes/literature-officer-20260906-issue601.md`, `...-round2.md` |
