# Forward analysis plan, 2026-09-06

Written **before** the next implementation block, at the owner's instruction, so that
what follows is executed against a declared plan rather than assembled after the fact.
Everything below is prospective: the decision rules are fixed here, before the numbers
that would decide them exist.

Scope of this plan: the next four work blocks, plus the paper track. It does not
re-open anything already decided, and it does not schedule new large production.

---

## 0. Where the research object actually stands

Two things changed in the last two blocks, and together they change what the next
question should be.

**#582/#581 (state contraction, threshold side).** The threshold law does not move by
center and width; the Q tangent splits at Boolean degree one; what is left is a small
*structured* remainder. #584 asks whether that remainder carries a discrete index.

**#580 (state contraction, exact-process side).** On P398 — an exactly known finite
noncrossing generator — a rank-6 observable-Krylov span frozen at `eta = 0` transports a
declared intervention to `eta = ±1/4` for `+0.006` excess relative error, stably from 14
states to 1430, and survives an intervention chosen deliberately outside the operator
pencil that built the span. But the same frozen realization does not *represent* five
readouts it was not built from, and that failure is already present at `eta = 0`.

The refinement pass added the part that matters here: the three rank notions do not
agree, and their disagreement grows.

```text
width   states   r_linear(D0)   r_positive   r_transport
  4       14         10             10            4
  5       42         26             26            4
  6      132         72             76            6
  7      429        218            232            6
  8     1430       >=150           750            6
```

`r_transport = 6` is flat. `r_positive = 750` at width 8. A factor of 125, widening.
And the exact lumping that certifies `r_positive` **collapses to the identity partition
under the out-of-pencil intervention at every width** — no nontrivial exact positive
realization survives it.

So the compact object is real, transportable, and *dictionary-relative*. It is a state
of its readouts, not a state of the system.

### The one question the current language cannot answer

The held-out failure has two completely different possible causes, and #580 cannot
separate them:

```text
(i)  the dictionary is missing an instantaneous state coordinate
        -> add state; a bigger Markov model is the right object;

(ii) the exact Markov process, projected onto an observable-relative subspace,
     generates memory
        -> the residual is history, not a hidden coordinate; adding state is
           the wrong repair, and any "fiber label" found elsewhere may be
           memory wearing a label.
```

This is the fork. Every downstream item — #584's index, #275's mechanism verdict,
#249's Hankel rank, and eventually any square-site reading — inherits whichever answer
is true. That makes it first in the queue, and it is the owner's #588.

---

## 1. Priority queue

Ordered by (information about the fork) / (compute + risk). Items 1 and 2 need **zero
new sampling and zero new width**.

### P1 — #588 Phase A: the exact Mori–Zwanzig kernel on the already-frozen span

The whole point is that P398 is exact, so `K(t)` is not estimated, it is computed.

With `P` the frozen #580 projection and `Q = I - P`, in a `(P,Q)`-adapted orthonormal
basis,

```text
L = [ A  B ]        K(t) = B exp(tD) C
    [ C  D ]

d x_R/dt = A x_R(t) + int_0^t K(t-s) x_R(s) ds + B exp(tD) x_U(0)
```

`K(t)` is *exactly* the dynamics the Markov closure threw away. Deliverables, in order:

- **A1 controls.** Full-rank `P` must give `K = 0`; an exactly invariant `P` must give
  `K = 0`; a synthetic block with planted `B,C,D` must reconstruct. And the block
  formula must be verified against the corresponding block of `exp(tL)` rather than
  trusted from notation — the repo stores observables in one convention and the
  literature in the other.
- **A2.** `||K_eta(t)||` on the existing lag grid, integrated memory weight, decay
  time / tail mass, kernel pole count, Markov-only response error, and
  Markov+exact-memory response error, for widths 4–8 on the declared `eta` ladder.
  The last is the closing control: exact memory must remove the projection error down
  to the separately tracked unresolved-initial-condition term.
- **A3.** The same for the out-of-pencil `single_point_join`. The question is not
  whether `K_eta` changes — it must — but whether the *baseline memory structure*
  transports, or whether the memory spectrum changes in kind.

**Reuse, not rebuild.** The span, projection, readouts, sources, lag grid, `eta`
ladder, uniformized propagator, modular rank, exact lumping and the readout-balanced
metric all already exist in `scripts/p398_intervention_transport.py` and are tested.
Phase A adds the block decomposition, the kernel, and the memory-corrected solve.

**Pre-registered decision rule** (fixed now):

```text
declared readouts short/weak memory AND held-out strong but LOW-ORDER memory
    -> the compact state is genuinely dictionary-relative; prefer state+memory
       over new latent coordinates, and warn #584 explicitly;

few memory poles Markovianizable by few COMMON auxiliary coordinates
    -> strongest finite version of base+fiber; freeze it and ship it to a new
       intervention before attaching any meaning;

memory rank/tail GROWS with width or dictionary
    -> the low-rank state is a compression of selected outputs; connect to
       predictive noncompression, stop enlarging the operator model;

memory small BUT held-out error still large
    -> the fault is the observable map / projection semantics, not history;
       fix the readout before touching state dimension.
```

**Phase B** (nested dictionaries D0/D1/D2, already declared in #580 — no new observable
may be invented after seeing kernels) reports `r_linear`, `E_markov`, `M_memory`,
`r_memory`, `r_transport` **separately per dictionary**. The scientific object is the
growth profile under refinement, never one number called "state dimension".

**Phase C** compares the two repairs at equal cost per degree of freedom: `S` = smallest
baseline-only Krylov/leakage augmentation, `M` = smaller state + fixed low-order rational
kernel learned from baseline data only. Plus the owner's added adversary — a
finite-horizon balanced reduction on the zero-sum contrast subspace,

```text
W_c(T) = int_0^T exp(tG) B_src B_src^T exp(tG^T) dt
W_o(T) = int_0^T exp(tG^T) C_out^T C_out exp(tG) dt
```

with the stationary mode removed explicitly (never fed to an infinite-horizon Lyapunov
solve), balancing transform frozen at `eta = 0` and transported exactly as in #580. This
is the strongest deterministic Markov competitor; if the balanced state stays compact and
transports, #580's Krylov choice was not special.

The output metric is frozen before balancing, and reported under **both** the pooled and
the readout-balanced norm — #580 already showed the pooled Frobenius norm is
magnitude-weighted and can rotate a verdict.

Risk: low. Cost: one script, no new sampling. This is the block to do next.

### P2 — #589/#591: derivative-channel Smith rescore of the existing #205 block

Zero new sampling. The preflight note establishes the confound exactly: for a square
Gaussian torus from `w = a + ib`, `d1 = gcd(a,b)`, `d2 = N/d1`, and at N=650 each
three-angle family carries one noncyclic row with `Smith = (5,130)`, so a
three-coefficient harmonic fit absorbs the quotient loading exactly into `(C, A4, A8)`.
Until that is calibrated, **an N=650 third coefficient cannot be read as pure `A8`.**

The cheapest calibration is to rescore #205's *archived* threshold-rank histograms in
the `Sp`/`Dp` derivative channels with the aligned 3×3 delete-one covariance. The exact
H4 nulls are already frozen in `notes/norm5-conjugate-coalescence-tomography.md`:

```text
N65 parents 8+i, 7+4i;  observed 2-i -> A=(17,6),  B=(18,1);  missing 2+i -> C=(15,10)
N85 parents 9+2i, 7+6i; observed 2+i -> A=(16,13), B=(19,8);  missing 2-i -> C=(20,5)
C has Smith invariants (5,65) and (5,85); A and B are cyclic.
N325:  5 M_C - 11 M_A +  6 M_B = 0
N425: 20 M_C + 13 M_A - 33 M_B = 0
```

**This is explicitly a re-view of the same evidence block, not a new independent vote**
(GOVERNANCE §2E). It is reported as a nuisance calibration attached to #205, and it
cannot by itself move the H4 verdict in either direction. Its only job is to say whether
the noncyclic row loads the derivative channels enough to contaminate an `A8` reading.

Blocking relation: **#583 (N=650) must not be interpreted before this lands.**

### P3 — #584: index or falsify the #582 remainder

Candidate indices, all predeclared: Smith type, deck group, Gaussian cover word, parent
lineage, primitive homology class, period/orientation class. Permutation controls on
every one; a label that survives no permutation control is not a label.

**#580 hands this a warning that must be carried in the note**: a structured residual
does not imply an omitted instantaneous coordinate. If #584 finds a plausible fiber
label, the very next question is whether it is a static typed coordinate or history
exposed by that context — which is exactly the discrimination P1 calibrates. So P3
runs *after* P1 by design, not by convenience.

### P4 — #275, rescoped

#580 changed the design: mechanism classes must be scored on dictionaries **known to
separate them**, not on a shared one. Two models sharing a dictionary sharing an image
is the expected outcome and reads `UNIDENTIFIABLE_WITH_CURRENT_ASSETS` — which is what a
shared-dictionary design would have produced regardless of the truth. Column-space /
profile-rank verdict follows the separating-dictionary construction, not before it.

### Standing, unscheduled

- #581 first empirical control: a square-bond block with 3–5 predeclared Euclidean
  scales and a predeclared `O_H4`.
- Merge chain #564 → #573 → #575 → #577; #576.
- Reference [6] still `[LIT]`; one deterministic N=580 replay settles `bare_aspect_ratio`.
- Scullard 2006 unread.

---

## 2. What this plan will not do

Carried forward and still binding:

- No top-up of #205.
- No further N130/N170 identical local pivotal rows.
- No new sizes for #537.
- No new free-exponent or scalar-width closure.
- No new width beyond 8 for P398 in #588 Phase A/B.
- No monitoring, watching, subscription or scheduled check-in without the owner
  saying so first.
- No repeated re-audits of settled blocks.

And one added here: **no new observable may be introduced into the P398 dictionary
after a kernel has been looked at.** D0/D1/D2 membership is closed as of #580.

---

## 3. Paper track

The point percolation ("点渗流") goal is the frame, and it is worth saying plainly where
the current results stand relative to it: none of #580/#582/#588 is a threshold result.
What they are is a *methodological* result about what a finite state means under partial
observation, on an exactly solvable process, with the arithmetic done exactly.

That is publishable on its own terms, and it is what the paper should claim. Working
title and spine:

```text
"A transportable state is a state of its readouts":
    exact intervention transport, positive realization, and projected memory
    on a finite noncrossing process.

1. The object: P398's exact generator, the affine family G_eta = G_0 + eta H.
2. Three rank notions that do not agree, computed exactly
   (modular elimination for r_linear, certified lumping for r_positive,
    frozen-span transport for r_transport), and their widening gap.
3. Transport is not an algebraic tautology: the out-of-pencil control.
4. But the positive realization sees the pencil: identity-partition collapse.
5. Representability, not dictionary membership, predicts transport.
6. Mori-Zwanzig: is the held-out failure missing state or projected memory?
7. Boundary: what this does not say about percolation.
```

Sections 1–5 are done and in tree. Section 6 is P1. Section 7 is the honest part and
must be written before, not after, the abstract.

Journal choice is mine per the owner's grant; deferred until §6 lands, because the
answer to the fork decides whether this is a methods paper or a negative-result paper,
and those are not the same venue.

---

## 4. Order of execution

```text
1. #588 Phase A          (this block)
2. #588 Phase B + C      (next block; C includes the balanced adversary)
3. #589/#591 rescore     (zero sampling; unblocks #583)
4. #584 with the P1 warning attached
5. #275 rescoped
```

Reported to the owner at each boundary, in Chinese, with the decision rule quoted from
this file rather than restated from memory.
