# #592 — Projected-memory literature check for #588 (P398 object)

**Branch:** `retrieval/p592-projected-memory-lit-20260913`
**Ticket:** #592  ·  **Source note:** `notes/p398-projected-memory-20260906.md` (read via `gh api`, decoded)
**Retrieval mode:** novelty judgement only. Retrieval is *not* evidence about any percolation claim. A negative search here is **not** proof of originality — see §6.

---

## 1. Context carried from P398 (frozen, read)

P398 studies a finite CTMC (non-crossing planar process, widths 4–8, 14→1430 states) whose
microscopic generator `G` is known exactly. With `P = ΦΦᵀ` the orthogonal projection onto a
rank-6 observable-Krylov span, `Q = I − P`, the projected equation is

```
A = ΦᵀGΦ      B y = ΦᵀG y        (y in range Q)
C a = QGΦ a  D   = QGQ
dx_R/dt = A x_R(t) + ∫_0^t K(t−s) x_R(s) ds + B exp(tD) x_U(0)
K(τ) = B exp(τD) C
```

Four core findings to judge for prior art:

- **F1 (order saturation).** Block-Hankel order carrying 99.9% of kernel energy is
  `2,3,4,4,4` across widths 4–8 (states 14→1430, a 102× growth) at the 99% level it is
  `2,3,3,3,3`. `rank C = 3` structurally (three of six `C` columns vanish identically).
- **F2 (seed-readout residual is memory).** On the declared dictionary the Markov-closure
  error is ~96% projected memory (four-way attribution, second-order convergence).
- **F3 (held-out readouts split).** On span-out readouts the closure error splits into
  ~40% memory / ~43% unrepresented initial data, stable across all widths.
- **F4 (lumping vs memory robustness).** The exact strong lumping (width 8: 750 blocks)
  collapses to the identity partition under rate perturbations *outside* the operator pencil,
  whereas the signed low-rank realization and the memory description are essentially unaffected.

---

## 2. Retrieval matrix

| # | Query | Source(s) reached | Outcome | Grade |
|---|-------|-------------------|---------|-------|
| Q1a | "Mori–Zwanzig memory kernel McMillan degree Hankel rank bound independent of environment dimension" | Gouasmi–Parish–Duraisamy 2017 (Proc. R. Soc. A 473:20170385); Ma–Li–Liu 2018 (arXiv:1802.10133); Hokanson 2018 (arXiv:1803.00043 "A data-driven McMillan degree"); Krishnaprasad lecture notes (UMD) | No theorem found bounding McMillan degree / Hankel rank of a Krylov-seed MZ memory kernel *independently of dim(D)*. | ABSTRACT_ONLY / [LIT] |
| Q1b | "Krylov subspace projection memory kernel order reduction balanced truncation singular values" | Ma–Li–Liu 2018 (arXiv:1802.10133, Krylov projection ↔ moment matching, FDT order<6); arXiv:1403.6543 (Krylov memory-function computation); Gouasmi 2017 | Krylov projection studied for MZ kernels; order discussed as Markovian-embedding (Padé) order, not as a dim(D)-independent McMillan bound. | ABSTRACT_ONLY |
| Q2 | "positive realization Mori–Zwanzig memory kernel compartmental finite-horizon balancing Markov chain" | Handwiki MZ formalism; Agarwal et al. 2023 (arXiv:2305.20083, MZ for Markov renewal processes, "compact when almost Markovian"); ResearchGate (Koopman/MZ learning) | No positive-realization theorem for the MZ memory kernel of a Krylov projection. Finite-horizon gramians exist as a concept but no MZ-CTMC-specific balancing theorem found. | ABSTRACT_ONLY / [LIT] |
| Q3 | (highest value, per ticket update) same object as Q1: is the order *saturation* a theorem or a finite-sample accident? | same as Q1a/Q1b + P398 internal arithmetic | Saturation in P398 explained structurally by `rank K(τ) ≤ rank C = 3` (a span identity), not by any published dim(D)-independent bound. | PRIMARY_TEXT_READ (P398) + ABSTRACT_ONLY |
| Q4 | "finite-horizon balancing / time-limited balanced truncation of Markov chains; MZ memory of CTMC" | Agarwal et al. 2023 (arXiv:2305.20083); standard balanced-truncation literature (Glover 1984, time-limited truncation) | Finite-horizon / time-limited balanced truncation is a known idea for LTI; not found applied to the MZ memory kernel of a CTMC with a Krylov seed. | ABSTRACT_ONLY / [LIT] |
| F2/F3 | "memory share four-way attribution projection closure error; MZ noise vs memory decomposition" | Gouasmi 2017; Ma–Li–Liu 2018; Parish–Duraisamy dynamic-τ; Stinis renormalized MZ | The exact four-term decomposition (markov / +memory / +forcing / +both) is a P398 construction; literature reports memory/noise split qualitatively, not a quantitative 4-way attribution on a held-out vs seed readout split. | ABSTRACT_ONLY |

Grade legend: `PRIMARY_TEXT_READ` = read full source body; `ABSTRACT_ONLY` = read abstract/page only;
`[LIT]` = standard/secondary textbook knowledge; `API_QUERY` = database query.

---

## 3. Answers to Q1–Q4 (each `cite-or-gap`)

**Q1 (highest value).** *Is there a theorem giving the McMillan degree / Hankel rank of the
Mori–Zwanzig memory kernel of a Krylov-seed orthogonal projection, with a bound independent of
the environment (orthogonal-complement) dimension?*
→ **GAP.** I found no such theorem. What exists: (i) the MZ formalism with a finite-rank
projection onto a Krylov/observable span (Mori 1965; Gouasmi 2017; Ma–Li–Liu 2018); (ii) the
general LTI fact that the McMillan degree of a transfer function equals the rank of its infinite
Hankel matrix (Ho–Kalman; Hokanson 2018). For the kernel `K(τ)=B e^{τD}C` the McMillan degree is
the dimension of the minimal `(D,B,C)` realization and is *a priori* bounded above only by
`dim(D)` — i.e. by the environment dimension, the opposite of the sought independence. No result
cuts this bound down to the Krylov dimension for an arbitrary projection. `cite-or-gap: GAP`
(related: [LIT], ABSTRACT_ONLY on Gouasmi/Ma–Li–Liu/Hokanson).

**Q2 (positive realization).** *Does the projected memory kernel admit a positive (Metzler)
realization?*
→ **GAP.** A kernel `K(τ)=B e^{τD}C` is the Laplace transform of a positive system iff
`(D,B,C)` is a positive realization (`D` Metzler, `B,C ≥ 0`). For a Krylov projection of a CTMC
generator `G` (Markov generator: off-diagonal ≥0, row sums 0) the projected `(D,B,C)` is not
obviously positive. No published positive-realization result for the MZ memory kernel of such a
projection was found. `cite-or-gap: GAP`.

**Q3 (highest value, per ticket update).** *Is the effective-order saturation in F1 a theorem or
a finite-sample accident?*
→ **Partial / GAP on the theorem, structural explanation available.** The 99.9%-energy order
`2,3,4,4,4` does **not** follow from a general dimension-independent bound (none exists, see Q1).
It is instead a direct consequence of the structural identity `rank K(τ) ≤ rank C` together with
`rank C = 3` for this span (three seed directions span their own first Krylov level, so only the
frontier survives). The *numerical* order (tol 1e-6) does grow sublinearly (12,13,14) while the
state space multiplies ~11×, so "saturation" is accuracy-target-dependent
(`BOUNDED_MEMORY_ORDER_IS_ACCURACY_DEPENDENT` in P398). Conclusion: the flat low order is a
property of *this span*, not a theorem about Krylov-seed MZ kernels in general.
`cite-or-gap: gap-on-theorem; structural explanation = PRIMARY_TEXT_READ (P398)`.

**Q4 (finite-horizon balancing of Markov chains).** *Can the CTMC memory be treated by a
finite-horizon balancing?*
→ **GAP (specific).** Finite/time-limited balanced truncation and finite-horizon gramians are
known for LTI (Glover; time-limited truncation). Agarwal et al. 2023 apply MZ to Markov renewal
processes and prove the representation is compact "when the MRP is almost fully Markovian"
(all-but-one kernel vanishes in the Markovian limit) — relevant but not a balancing theorem for
the Krylov-seed CTMC kernel. No MZ-CTMC finite-horizon balancing theorem found.
`cite-or-gap: GAP` (ABSTRACT_ONLY on Agarwal 2023; [LIT] on balanced truncation).

---

## 4. F2 / F3 / F4 — prior-art status

- **F2 (96% memory on seed readouts).** The *exact four-term* closure attribution is a P398
  device. Literature establishes qualitatively that most reduced-model error is memory
  (Gouasmi 2017; Parish–Duraisamy; Stinis) but does not give the quantitative 96% seed-readout
  figure. `cite-or-gap: no-prior-art-for-this-exact-attribution` (the underlying MZ memory
  dominance is well known → [LIT]).
- **F3 (40/43 split on held-out).** The stable *split* between memory and unrepresented initial
  data on span-out readouts is not reported in the searched MZ literature; it is specific to the
  P398 experiment. `cite-or-gap: GAP`.
- **F4 (lumping collapses, memory robust).** Strong-lumping sensitivity to perturbations outside
  the operator pencil vs robustness of low-rank/memory descriptions: no directly matching theorem
  found. Perturbation sensitivity of lumping is classical (Kato; Simon) but the specific
  pencil-outside collapse vs memory-robustness comparison is not in the retrieved literature.
  `cite-or-gap: GAP` (classical perturbation theory = [LIT]).

---

## 5. Explicit negative results (what was NOT found)

- No theorem bounding the McMillan degree / Hankel rank of a Krylov-seed MZ memory kernel
  independently of the environment dimension (Q1, Q3).
- No positive-realization theorem for the projected MZ memory kernel (Q2).
- No MZ-CTMC finite-horizon balancing theorem (Q4).
- No published exact four-term (markov / +memory / +forcing / +both) closure attribution, nor the
  specific 96% / 40%–43% split figures (F2, F3).
- No published result stating that exact strong lumping collapses to the identity under
  pencil-outside perturbations while the memory description does not (F4).
- The 2026-vintage loop-model papers (arXiv:2604.24491 etc.) were **not** part of this ticket's
  object; see #585 for that literature.

Searched endpoints: Google Scholar / arXiv / Royal Society / AIP / ScienceDirect / handwiki /
ResearchGate. No endpoint was blocked by a proxy during this run.

---

## 6. Mandatory disclaimer

This note performs **novelty judgement only**. It is **not** evidence about any percolation
threshold, square-site or otherwise. A negative retrieval (GAP) means "not found in the sources
consulted," **not** "does not exist / is original." Originality would require an exhaustive
database/catalogue search (e.g. MathSciNet, the MZ / model-reduction primary literature in full),
which is broader than the WebSearch/WebFetch sweep done here. The structural explanation of F1
(rank C = 3) is read from P398 itself, not from external literature.

---

*Full Matching-One repository CI has not been run for this commit.*
