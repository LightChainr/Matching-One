# Literature retrieval: double-pulse / bilinear kernels against the P398 odd-sector certificate

**Date:** 2026-09-12
**Ticket:** #713 (parent #650). Context: draft #709 on `analysis/p398-double-pulse-20260912`.
**Type:** retrieval only. No new production, no enumeration, no transfer matrix, no hardware.
**Does not:** merge, close, or edit `docs/STATUS.md`. No novelty claim is made anywhere below.

Object under discussion, quoting draft #709's own definitions (not re-derived here):

```text
G = sum_j (J_j + D_j),   R : j -> w-1-j,
H = J_0 - J_(w-2),       RGR = G,  RHR = -H,
K(tau) = S H exp(tau G) H F = C_- exp(tau G_-) B_-.
```

The four questions in #713 are about what earlier literature already contains for this shape.
Every source is marked **PRIMARY_TEXT_READ**, **ABSTRACT_ONLY**, or **[LIT]** (bibliographic
pointer; body not retrieved). Quotations below are verbatim from the retrieved text only.

## Tripwire answers (stated first, because they bound the rest)

1. **Bilinear/Hankel prior art claimed as PRIMARY?** No. The only bilinear-Hankel text I read in
   full is Petreczky, ESAIM: COCV **17** (2011) 446–471, Part II (theorem numbers below). Fliess
   1974, Isidori 1973 and Ho–Kalman 1966 are **not** upgraded to PRIMARY here: they are
   ABSTRACT_ONLY or [LIT]. Arbib–Manes 1980 is ABSTRACT_ONLY; it names the Fliess/Isidori Hankel
   matrices but I did not open those originals.
2. **Does the P398 double-pulse certificate itself appear in print?** **No.** Nothing retrieved
   matches the specific scalar observer `Khat(z) = 2(z^2+11z+27)/[(z^2+10z+23)(z^2+11z+26)]`, the
   4×4 determinant `−16`, or the width-4/5 odd-sector rank attainment on `G = sum_j (J_j+D_j)`
   with the original `S`, `F`. The expected answer holds. Absence of a retrieved match is **not**
   an originality certificate — it is the only statement the search supports.

## Evidence legend

| Mark | Meaning in this note |
|---|---|
| **PRIMARY_TEXT_READ** | I retrieved and read the full body (arXiv HTML, or publisher/Numdam PDF text extracted with `pdftotext`). Section/equation/theorem numbers are from that body. |
| **ABSTRACT_ONLY** | Only the abstract and landing-page metadata were retrieved. Equations attributed to the paper are quoted from its own abstract, not reconstructed. |
| **\[LIT\]** | Bibliographic pointer only (title/venue/year). The body was not retrieved; no statement of its content is made here. |

---

## Q1. Two-pulse / Volterra kernels recovering a hidden sector after a quotient

**What was retrieved and read.**

Lucarini, *Interpretable and Equation-Free Response Theory for Complex Systems*,
**arXiv:2502.07908v2** HTML, §II.2 (§IV.1.3 for the Volterra remark).
**PRIMARY_TEXT_READ.** Published version: *Phil. Trans. R. Soc. A* **384** (2026) 20250081,
DOI 10.1098/rsta.2025.0081 (**ABSTRACT_ONLY**; the DOI landing page was not opened here).

This is the closest retrieved match to Q1's "two-pulse kernel for a Markov chain". For a
finite-state Markov chain with perturbation `M -> M + eps f(n) m`, the paper's second-order
response is a genuine two-time kernel. Verbatim from the retrieved HTML:

- second-order measure, Eqs. (13)–(14):
  `nu^(2)(n) = sum_k sum_p Theta(k)Theta(p) M^k m M^p m nu_inv f(n-k-p-2) f(n-k-1)`,
  with the text "the Θ's ensure the correct time ordering of the perturbation";
- second-order Green's function, Eq. (18):
  `G^(2)_{m,Psi}(k,p) = Theta(k)Theta(p) <m^T (M^T)^p m^T (M^T)^k Psi, nu_inv>`,
  called a "double convolution sum";
- **two-mode spectral decomposition**, Eq. (19):
  `G^(2)_{m,Psi}(k,p) = Theta(k)Theta(p) sum_{i,j=2}^{N} alpha_{ij} lambda_i^k lambda_j^p`,
  `alpha_{ij} = <m^T Pi_j m^T Pi_i Psi, nu_inv>`;
- §IV.1.3: "The nonlinear Green functions can be formally seen as **Volterra kernels**".

**What this does and does not supply for P398.**

- It supplies the general, already-published statement that a Markov chain's second-order
  response is governed by a two-time (Volterra-type) kernel with two-mode spectral weights.
  P398's `K(tau) = S H exp(tau G) H F` is an instance of this language, **not a new response
  formalism**. This matches #709's own disclaimer.
- It does **not** single out a parity/odd sector. Lucarini's decomposition is by Koopman modes
  `lambda_i`, not by a `C2` (reflection) irrep. The P398 selection rule that an invariant task is
  blind to `H` at first order is a representation-theoretic statement (trivial irrep in a tensor
  product, as in #598/#244), and I did not retrieve a primary paper that phrases a two-pulse
  hidden-sector recovery as a parity quotient of a Markov generator.
- Mueller, Basu, Sollich, Krueger, *Coarse-grained second-order response theory*, Phys. Rev.
  Research **2** (2020) 043123, DOI 10.1103/PhysRevResearch.2.043123, remains **ABSTRACT_ONLY**
  (already so in #709). Its equilibrium hypotheses are not imported here.

**Finding for Q1.** The *machinery* is published (two-time second-order kernels, spectral
two-mode decomposition, explicit Volterra framing). A primary source that uses a two-pulse kernel
to recover a *parity-odd hidden sector after an even quotient* in a Markov chain or interacting
particle system was **not found**. Treat that as a gap in my retrieval, not as a no-go.

---

## Q2. Bilinear realization, Fliess kernels, and Hankel rank of `C A^k B`

### 2a. The verified primary text

Petreczky, *Realization theory for linear and bilinear switched systems: a formal power series
approach. Part II: Bilinear switched systems*, **ESAIM: COCV 17 (2011) 446–471**, DOI
10.1051/cocv/2010015. **PRIMARY_TEXT_READ** (Numdam PDF, `pdftotext -layout`).

Word coefficients and the Hankel matrix, §2.1 (verbatim):

```text
c_f(eps) = C x0,   c_f(j1 j2 ... jk) = C B_jk B_j{k-1} ... B_j1 x0     (2.6)
```

> "columns and rows of which are indexed by sequences v in Z*_m. The p × 1 block entry of H_f
> lying on the intersection of the row indexed by v and the column indexed by w equals c(wv). It
> turns out that the generating series c has a representation of the form (2.6) if and only if the
> column rank of H_f is finite. That is, f has a realization by a bilinear system if and only if f
> has a Fliess-series expansion and the column rank of its Hankel-matrix H_f is finite."

Definition 2.8 defines the Hankel matrix `H_Phi` with block entry
`(H_Phi)_{(v,i),(w,f)} = (T_{f,sigma_{K+1}}(wv))_r`, and rank = dimension of the column span.
The minimal-realization characterization, **Theorem 2.3(iii)** (verbatim):

> "The dimension of Sigma_min equals the rank of the Hankel-matrix H_Phi of Phi, i.e.
> dim Sigma_min = rank H_Phi."

**Theorem 2.7** (existence, constrained switching): finite Hankel rank plus generalized
Fliess-series expansion gives a realization. **Remark 2.5**: the proofs construct the realization
from the columns of the Hankel matrix.

Petreczky, Wisniewski, Leth, *Moment matching for bilinear systems with nice selections*,
**arXiv:1605.04414v1**, IFAC-PapersOnLine 49(18):838–843 (2016), DOI
10.1016/j.ifacol.2016.10.270. **PRIMARY_TEXT_READ** (arXiv HTML). This is the paper #709 already
cites. It gives the bilinear system (1a)–(1b), the matrix words `A_w = A_{q_k}...A_{q_1}` with
`A_eps = I`, the coefficient identity `c_f(w) = C A_w x0`, the Fliess operator/series (Eq. (4)),
and the `gamma`-partial-realization construction (Def. 2, Theorems 3–4).

Important negative reading, verified in the same body: **this 2016 paper does not state a
Hankel-rank criterion.** Its minimality test is algebraic — "Sigma is a minimal realization of f
if and only if Sigma is a realization of f, and it is span-reachable and observable", with
observability `intersection_w ker C A_w = {0}` and span-reachability `Span{A_w x0 | w in Q*} = R^n`.
So it is not the place to attribute the rank theorem to; the 2011 Part II is.

### 2b. The classical linear case (`K(tau) = C exp(tau A) B`)

The shape `H_i^j = C A^{i+j} B` and the equality "minimal order = Hankel rank" are the classical
linear realization theorem.

- Arbib, Manes, *Generalized Hankel Matrices and System Realization*, **SIAM J. Math. Anal. 11(3)
  (1980) 405–424**, DOI 10.1137/0511038. **ABSTRACT_ONLY.** The abstract states verbatim: "Our
  definition of the Hankel matrix unifies the familiar `H_i^j = CA^{i+j}B` of linear system theory
  ... with the bilinear Hankel matrix of A. Isidori ... and the Hankel matrix of M. Fliers
  (Matrices de Hankel, J. Math. Pure Appl., 53 (1974), pp. 197–224)." It also states a
  realisability theorem, a partial-realization theorem, and a canonical-realization theorem for
  finite Hankel blocks.
- Fliess, *Matrices de Hankel*, J. Math. Pures Appl. (9) **53** (1974) 197–222 (the Arbib–Manes
  abstract prints the range as 197–224). **[LIT]** — body not retrieved.
- Isidori, *Direct construction of minimal bilinear realizations from nonlinear input-output maps*,
  IEEE Trans. Autom. Control **AC-18** (1973) 626–631. **ABSTRACT_ONLY** (IEEE/scilit landing
  pages; the bilinear Hankel matrix is attributed to it by Arbib–Manes, not read here).
- Ho, Kalman, *Effective construction of linear state-variable models from input/output functions*,
  Regelungstechnik **14** (1966) 545–548. **[LIT]** — not retrieved; named for attribution only.

### 2c. What is *not* claimed as new if we only certify a finite instance

A finite 4×4 or 16×16 nonzero Hankel minor is an **instance** of Theorem 2.3(iii)/Theorem 2.7 of
Petreczky Part II, specialised to `C A^k B` with `A = G_-`, `B = B_-`, `C = C_-`. If #709's
deliverable is read as "certified a finite instance", then all of the following are prior art and
must not be presented as new:

- Fliess-series / word-coefficient representation `c(w) = C A_w x0`;
- the Hankel-matrix construction and the finite-rank ⇔ realizable equivalence;
- the identity "minimum realization order = rank of the Hankel matrix";
- partial-realization-from-Hankel-blocks.

The only candidate for a new statement is the **model-specific** one: that P398's *particular*
`S`, `F`, `G`, `H` at `w = 4, 5` attain the parity upper bounds with a nonzero integer minor, and
that the **equal-time** mixed derivative is zero while the delay kernel has full odd rank. Per the
tripwire, that specific certificate is not in print (Q4/tripwire).

---

## Q3. Noncrossing partitions / Temperley–Lieb / join-split generators as CTMCs

Two named families were retrieved; neither is P398's generator.

### 3a. Temperley–Lieb stochastic processes (closest named relative)

Pearce, Rittenberg, de Gier, Nienhuis, *Temperley-Lieb Stochastic Processes*,
**arXiv:math-ph/0209017v2**, J. Phys. A **35** (2002) L661–L668, DOI
10.1088/0305-4470/35/45/105. **PRIMARY_TEXT_READ** (arXiv HTML).

The abstract generator is Eq. (2.1): `H = sum_a c_a (1 - w_a)`, `c_a >= 0`, an intensity matrix
satisfying the master equation (2.2) `dP_a/dt = -sum_b H_ab P_b`. For the Temperley–Lieb algebra
the generator is Eq. (2.6) `H = sum_{j=1}^{L-1} (1 - e_j)`; the **cylindrical** version (Eq.
(2.15)) adds the closed bond, `H = sum_{i=1}^{L} (1 - e_i)`. The state space is the set of
**link patterns / connectivities** (words in an ideal of the TL algebra), dimension given by (2.9)
with the Catalan count (2.13).

Crucially, verbatim: "**the terms in the Hamiltonian may connect disconnected lines but it is not
possible to have the reverse process**". So this generator is a **join-only** (irreversible)
dynamics on planar connectivities.

**Resemblance and difference for P398.** The *cylindrical* TL generator `sum_{i=1}^{L}(1-e_i)` is
structurally the nearest published object to a sum of local generators on *cyclic* planar
connectivities, which is why it is worth naming. But it is not P398's `G`:

- TL acts on link patterns (noncrossing **matchings**), P398 on noncrossing **partitions** of `w`
  cyclic points;
- TL's `sum (1-e_j)` is **join-only**, while P398's `G = sum_j (J_j + D_j)` contains the detach
  half `D_j` as well. The two are therefore different generators, and the TL result does not
  transfer as a statement about `G`.

### 3b. Exchangeable fragmentation–coagulation (both directions, but exchangeable/nonlocal)

Bertoin, *Two-parameter Poisson-Dirichlet measures and reversible exchangeable
fragmentation-coalescence processes*, **arXiv:0704.3122v1** (2007). **PRIMARY_TEXT_READ**
(arXiv HTML, §2.2).

An EFC process is a `P_N`-valued exchangeable Markov process whose restriction to `P_[n]` "only
evolves by fragmentation of one block or by coagulation"; the allowed first jumps are "obtained
from `pi^[n]` by **splitting exactly one of its blocks** ... or by **merging at least two of its
non-empty blocks**". This *does* have both join and detach, and can be reversible. But it is
**exchangeable** (non-local, all partitions, invariant under relabelling), so it is not the
planar/local object P398 uses.

### 3c. Answer: is `G = sum_j (J_j + D_j)` a named process?

**No named process matching P398's generator was found.** The retrieved named relatives are:
TL stochastic processes (join-only, matchings, elliptic/local but irreversible) and EFC /
coalescent–fragmentation (join+detach and reversible, but exchangeable and nonlocal). P398's
combination — **adjacent join plus point-detach, on cyclic noncrossing partitions, both
directions, local** — was not identified with a name in the retrieved literature. Pitman and
Diaconis appear in the surrounding noncrossing-partition and coalescent literature, but no primary
noncrossing-partition CTMC with this join+detach generator was retrieved; those remain **[LIT]**
pointers only. **Do not call `G` "the Temperley–Lieb process" or "an EFC process."**

Consistent with repo policy: an unretrieved name is not evidence that the process is unnamed, and
it is not an originality claim.

---

## Q4. A published "zero delay misses a delay-mode" example of the same shape

**Linear-systems side (verified shape, not a special example).** For `K(tau) = C e^{tau A} B`, the
first Hankel block is `C B` and later blocks are `C A^{i+j} B` (Arbib–Manes **ABSTRACT_ONLY**,
quoted above). `K(0) = 0` with nonzero Hankel rank is the ordinary **strictly proper / relative
degree >= 1** case: minimal order equals the full Hankel rank regardless of the vanishing of the
first Markov parameter (Ho–Kalman **[LIT]**; Petreczky Part II **PRIMARY_TEXT_READ** for the
bilinear version). In other words, "zero at zero delay, nonzero rank afterwards" is **generic**, not
a named counterexample, and certifying one such instance carries no novelty by itself.

**Two-pulse / Markov-chain side.** I found **no published example** of exactly the P398 shape — a
two-pulse Markov-chain delay kernel with `K(0) = 0` and Hankel rank 4 on an odd sector — in the
retrieved literature. The nearest physical analogue I noticed in search (a two-photon/dark-state
transition with zero one-photon matrix element, i.e. a selection-rule analogue) was **not retrieved
as a primary text** and is listed here only as a **[LIT]** pointer, not a claim.

**Finding for Q4.** The *linear* shape is standard and prior; the *P398-specific* two-pulse
certificate is not in print.

---

## Sources, with marks

| # | Source | Mark | Retrieved |
|---|---|---|---|
| 1 | Lucarini, *Interpretable and Equation-Free Response Theory for Complex Systems*, arXiv:2502.07908v2; Phil. Trans. R. Soc. A 384 (2026) 20250081 | **PRIMARY_TEXT_READ** (arXiv HTML; §II.2 Eqs. 10–19, §IV.1.3) | https://arxiv.org/html/2502.07908v2 |
| 2 | Lucarini, journal version, DOI 10.1098/rsta.2025.0081 | **ABSTRACT_ONLY** | (via search result; not opened) |
| 3 | Petreczky, *… Part II: Bilinear switched systems*, ESAIM: COCV 17 (2011) 446–471 | **PRIMARY_TEXT_READ** (§2.1 Eqs. 2.6; Def. 2.8; Thm 2.3(iii); Thm 2.7; Rem. 2.5) | https://www.numdam.org/item/COCV_2011__17_2_446_0.pdf |
| 4 | Petreczky, Wisniewski, Leth, *Moment matching for bilinear systems with nice selections*, arXiv:1605.04414v1 / IFAC-POL 49(18):838–843 | **PRIMARY_TEXT_READ** (§2.1–2.2, §3 Def. 2, §4 Thms 3–4; no Hankel criterion) | https://arxiv.org/html/1605.04414v1 |
| 5 | Pearce, Rittenberg, de Gier, Nienhuis, *Temperley-Lieb Stochastic Processes*, arXiv:math-ph/0209017v2; J. Phys. A 35 (2002) L661 | **PRIMARY_TEXT_READ** (Eqs. 2.1–2.3, 2.6, 2.9, 2.13, 2.15) | https://arxiv.org/html/math-ph/0209017v2 |
| 6 | Bertoin, *Two-parameter Poisson-Dirichlet measures and reversible exchangeable fragmentation-coalescence processes*, arXiv:0704.3122v1 (2007) | **PRIMARY_TEXT_READ** (§2.2) | https://arxiv.org/html/0704.3122v1 |
| 7 | Arbib, Manes, *Generalized Hankel Matrices and System Realization*, SIAM J. Math. Anal. 11(3) (1980) 405–424 | **ABSTRACT_ONLY** | https://epubs.siam.org/doi/10.1137/0511038 |
| 8 | Isidori, *Direct construction of minimal bilinear realizations…*, IEEE TAC AC-18 (1973) 626–631 | **ABSTRACT_ONLY** | https://ieeexplore.ieee.org/document/1100424 |
| 9 | Fliess, *Matrices de Hankel*, J. Math. Pures Appl. (9) 53 (1974) 197–222 | **[LIT]** | (cited via #7's abstract) |
| 10 | Ho, Kalman, *Effective construction of linear state-variable models from input/output functions*, Regelungstechnik 14 (1966) 545–548 | **[LIT]** | https://ntrs.nasa.gov/citations/19670049337 (pointer) |
| 11 | Mueller, Basu, Sollich, Krueger, *Coarse-grained second-order response theory*, PRR 2 (2020) 043123 | **ABSTRACT_ONLY** (as in #709) | https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.043123 |

## What this note does not do

- No novelty claim, in either direction. A gap in my retrieval is not a no-go theorem, and an
  unretrieved name is not an originality certificate.
- No new computation, no width campaign, no transfer matrix, no hardware.
- No `docs/STATUS.md` edit; no merge; no close of #713 or #709.
- The unique P398 science (the certificate itself) stays in the #709 GitHub PR, not here.
