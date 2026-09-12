# #679 — Does digital Alexander give r_b + r_w = 2 on every honest square-cell torus?

2026-09-08. WorkBuddy CLI (local machine). Notes/proof only — no Monte Carlo, no new
census, no Huawei, no STATUS edit, does not close #613/#632/#668. Parent #650.

**Verdict in one line: yes — `r_NN + r_NN+NNN = 2` is a theorem, configurationwise,
on every honest square-cell torus; it is proved in PR #271
(`notes/digital-alexander-duality-proof.md`) and the site L=3,4 / bond L=3 census
rows are finite *checks* of it, not its evidence base. The 16-pattern certificate
is only the local half; the "`= 2`" is the global half and needs one further
lemma, which we make explicit below.**

---

## 0. Honest tori: `ell_N > sqrt(2)` iff four distinct corners per cell

Let `T^2 = R^2/Lambda` with `Lambda <= Z^2` an index-`N` sublattice, so the
standard unit-square grid descends to a periodic square-cell decomposition with
`N` cells and `N` vertices. The four corners of the cell `[0,1]^2` are the
cosets of `(0,0), (1,0), (0,1), (1,1)`, and the differences between two corners
of (possibly different) cells range exactly over

```text
(+/-1, 0), (0, +/-1), (+/-1, +/-1),
```

which are precisely the eight vectors of `Z^2` with Euclidean length `<= sqrt(2)`.
Hence two corners are identified in the quotient iff `Lambda` contains a nonzero
vector of length `<= sqrt(2)`, i.e. iff `ell_N <= sqrt(2)`:

```text
four distinct corners per embedded cell   <=>   ell_N > sqrt(2).
```

Call such quotients **honest**. Axis `L x L` tori are honest iff `L >= 2`.
Degenerate (self-identifying) quotients are *excluded* by the theorem, not
refuted; see section 6.

## 1. What the census establishes, and what it cannot

The committed exact censuses confirm, with zero failures:

| source | scope | rank pairs `(0,2)/(1,1)/(2,0)` |
|---|---|---|
| site L=3 (512 configs) | honest axis | 259 / 162 / 91 |
| site L=4 (65 536 configs) | honest axis | 36 559 / 19 932 / 9 045 |
| bond L=3 (262 144 configs, after the #653 repair) | honest axis | 75 460 / 111 224 / 75 460 |

(`results/homological-balance-exact-torus/latest.json`, independently rerun in
`notes/literature-officer-20260906-homological-balance-verify.md`; bond row from
PR #653.) A null census through any finite `L` cannot distinguish "structural"
from "accidental" — that is exactly the worry PR #676 raised about its own
Conjecture 3.2. The question is therefore whether a proof exists in scope.
It does, and it predates the census-based uses.

## 2. (a) The local statement: complementary neighbourhoods, face by face

Fix any configuration on an honest torus. Let `G_B` be the black
nearest-neighbour graph (embedded, because cells are honest) and let `U` be a
closed regular neighbourhood of `G_B`; put `V = closure(S \ U)`. Let `W` be the
white matching graph: white NN edges plus both diagonals of every square.

The **local** claim is facewise, over the 16 black/white corner patterns:

> In each face, every active matching diagonal either (i) has its endpoints
> exactly the two white corners of the face — retained (this happens for
> exactly the two opposite-white-pair patterns, masks 5 and 10); or (ii) has
> its endpoints joined by a white NN path along the boundary of that same
> face — replaced by that path; and no pattern retains crossing diagonals.

Replacing each redundant diagonal by its boundary path changes every cycle by a
chain contained in a single contractible face, hence preserves classes in
`H_1(S)`. Call the result `G_W`. Then `G_W` is embedded, `G_W subset V`, and
cell by cell `G_W` is a spine of `V` (each face's part of `V` is a disk or an
annulus-adjacent piece that deformation-retracts onto its `G_W` part; glued
over faces, `V` retracts to `G_W` up to disks). Consequently

```text
im[H_1(W)   -> H_1(S)] = im[H_1(G_W) -> H_1(S)] = im[H_1(V) -> H_1(S)].   (*)
```

This is the precise role of the 4/8 complementary-adjacency convention, and it
is *purely local*: nothing about ranks, and nothing about the number 2, has
been used. The machine certificate
(`scripts/digital_alexander_local_bridge.py` against
`analysis/digital_alexander_local_bridge_manifest.json`, merged in PR #271)
checks all 16 patterns; re-run on 2026-09-08 it reports 16/16 pass, zero
replacement failures, retained diagonals exactly `{mask 5: (1,3), mask 10:
(0,2)}`, `all_local_cases_pass = true`, bit-identical to the committed
`results/digital-alexander-local-bridge/latest.json`.

## 3. (b) The global statement: the one-sentence bridge, made explicit

The global input is the complementary-subsurface duality lemma, quoted from
`notes/digital-alexander-duality-proof.md` section 2 (PR #271):

> Let `U, V` be complementary compact subsurfaces of a closed oriented surface
> `S`, meeting on their common boundary. Over `Q`, put
> `A = im[H_1(U) -> H_1(S)]`, `C = im[H_1(V) -> H_1(S)]`.
> Then `C = A^perp` for the nondegenerate intersection pairing on `H_1(S)`.

(Proof chain: Poincare duality identifies `A^perp` with
`ker[H^1(S) -> H^1(U)]`; the long exact sequence of `(S,U)` identifies that
kernel with the image of `H^1(S,U)`; excision gives `H^1(V, boundary V)`;
Poincare-Lefschetz gives `H_1(V)`; naturality matches the maps.)

The bridge from (a) to (b) — the sentence the ticket demands — is exactly:

> `r_NN = rank A` because `U` deformation-retracts to `G_B`;
> `r_NN+NNN = rank C` by the local bridge (*);
> `A` and `C` are orthogonal complements by the lemma;
> `dim H_1(T^2; Q) = 2`;
> therefore `r_NN + r_NN+NNN = rank A + rank A^perp = 2`, for every
> configuration on every honest square-cell torus. `[]`

The only torus-specific input is the dimension of `H_1`. On a closed oriented
surface of genus `g` the *identical* local certificate would give
`rank A + rank C = 2g`: the local statement constrains the two images, and the
global statement fixes their sum. Neither half alone yields the identity.

Allowed rank pairs, exhaustively: `(0,2), (1,1), (2,0)`. In particular no
configuration has `(0,0)` or `(2,2)`, which is the structural half that PR
#676's Conjecture 3.2 leans on (see section 5).

## 4. Bond is not site (#646/#653 vs #271)

The bond-side identity at L=3 rests on a **different mechanism**: geometric
dual transport `T` (`notes/square-bond-transport-parity-theorem.md`,
`notes/square-bond-duality-tiny-torus.md`) is a bijection on bond
configurations that swaps primal and dual wrapping, so the odd combination `D`
satisfies `E[D] = 0` at the self-dual point `p = 1/2` — duality-oddness under a
measure-preserving involution. It is structurally tied to the square-bond torus
being self-matching (the dual grid *is* the primal grid).

The square-site pair (NN vs NN+NNN) is not self-matching — exactly,
`M_L(1/2) = -21/64` at L=3 and `-13757/32768` at L=4 — admits no such transport
involution, and its rank-sum is proved by the digital-Alexander route of
sections 2-3, not by oddness. Conversely, the #271 theorem says nothing about
bond observables. Neither direction imports. Both identities are
configurationwise and deterministic; the probability statements (such as
`E[D] = 0` at `p = 1/2` for bond) are additional and separate.

## 5. Scope against the literature (#613 item 4 convention)

- **Duncan-Kahle-Schweinhart (arXiv:2011.11903, AIHP 2025).** Supplies the
  ambient-`H_1` "giant-cycle" observable framework and the sharp-threshold
  mechanism for the square-bond torus. It neither states nor needs the site 4/8
  matching rank-sum; the proof obligations do not overlap. The rank identity
  here is a finite-`L`, configurationwise input; DKS-style asymptotics sit on
  top of it and are out of scope.
- **Classical topology.** The subsurface duality lemma is classical
  Poincare-Alexander duality for complementary regular-neighbourhood
  subsurfaces. The repository claims no novelty for it. The repository-specific
  obligation — that the 4/8 digital matching complement has the correct image
  in `H_1(S)` despite crossing diagonals and boundary-redundant diagonals — is
  discharged by the 16-pattern certificate (PR #271).
- **Cote-Uzcategui-Aylwin (arXiv:2503.17861).** Modern digital-connectivity
  treatment of 4/8 adjacency; cited as consistent context, not as a dependency
  of the proof.
- **van den Berg caveat (recorded in #670 item 3).** Concerns
  matching-lattice critical-point *relations* under dependence or
  generalization. Irrelevant to the configurationwise deterministic identity
  proved here — but equally, the identity must not be levered into a `p_c`
  statement without the #613/#614 probability inputs (exponential decay for
  site on both graphs; the amenable matching relation via Grimmett-Li). No
  such lever is pulled in this note.
- **Verdict per #613 item 4: repository lemma, not a new theorem.** Proved in
  PR #271 (merged 2026-08-29); the census rows are checks.
- **Effect on PR #676.** Conjecture 3.2's structural half — "`r_b + r_w = 2`
  at every `L`, hence no rank-2 x rank-2 and no rank-0 x rank-0 cell" — is
  exactly the #271 theorem on honest tori and travels to every honest `L`. The
  per-`L` wrap-cell support claim (exclusive crosses at each `L`; spiral
  cancellation `a_k = Delta #(both-same)`) remains genuinely conjectural: it is
  a statement about wrap-cell *counts*, which the rank theorem does not see.
  A-prime and the proposed `F = (1+M)/2` identity may use
  `r in {0,1,2}` with complementary ranks at every honest `L`; any remaining
  gap in those proposals lives in the count-level claims (and, per #668, in
  sector rank-purity), not in the rank-sum.

## 6. Degenerate quotients: excluded, not refuted

The hypothesis "honest" excludes short-period/self-identifying quotients
(`ell_N <= sqrt(2)`; 47 of the 140 HNF representatives through index 13). For
these the cellwise proof does not apply, because some face fails to embed with
four distinct corners and the regular-neighbourhood/spine construction breaks
down. The frontier oracle
(`notes/digital-alexander-short-period-frontier.md`,
`results/digital-alexander-quotient-frontier/latest.json`) exhausts all 140
representatives through index 13 — 101,140,028,118 complete filtrations — with
**zero** `rank_sum` failures across honest and self-identifying geometries
alike. That is evidence, not proof: a degenerate-quotient counterexample, if
one exists, would delimit the honest-cell theorem rather than contradict it.

## 7. Non-claims

No Monte Carlo; no new census; no threshold value, rate, or scaling claim; no
CFT field identification, no `V_(2,2)` selection rule; no STATUS edit; no
ticket closed. The identity is deterministic and probability-free; the #276/#613
qualitative-convergence targets remain exactly where they were.

## References

- PR #271 (merged): `notes/digital-alexander-duality-proof.md`,
  `scripts/digital_alexander_local_bridge.py`,
  `analysis/digital_alexander_local_bridge_manifest.json`,
  `results/digital-alexander-local-bridge/latest.json`.
- #613 (item 4) and PR #670: repository-lemma vs new-theorem convention.
- #676: Conjecture 3.2 and the structural/per-`L` split.
- #653: repaired bond L=3 census; #646: bond rank defect fix.
- `notes/digital-alexander-rank-oracle.md` (weak vs strong identity, finite
  oracles); `notes/digital-alexander-short-period-frontier.md` (degenerate
  frontier); `notes/square-bond-duality-tiny-torus.md` and
  `notes/square-bond-transport-parity-theorem.md` (bond transport parity);
  `notes/literature-officer-20260906-homological-balance-verify.md`
  (independent L=3,4 rerun).
- Duncan-Kahle-Schweinhart, arXiv:2011.11903 (AIHP 2025); Cote and
  Uzcategui-Aylwin, arXiv:2503.17861; Mertens and Ziff, arXiv:1603.07289.
