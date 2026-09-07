# C7: exact threshold-identifiability no-go — a finite even task cannot constrain the hidden sector

**Status:** theorem + exact numerical certificate (two-copy C2 construction).
Answers #599 program F and #594's comment-1 no-go target in the sharpest form
the probe can make: a *declared finite task* is a functional of one symmetry
sector only, and everything outside it — including a critical point placed
there — is invisible by construction.

Script: `scripts/probe/c7_no_go.py`; data:
`results/probe-c7-no-go/latest.json`.

---

## 1. Construction (exact)

Let `X = {0,1} x V` with `|V| = m`, and let `R` swap the two copies.  Even
sector = orbit functions (dim `m`); odd sector = sign functions (dim `m`).
A generator commuting with `R` is block diagonal on even/odd; even-task data
live entirely in the even block.  For a coupling `c in [0,1]` (the cross-copy
jump rate; same-copy rate `1-c`) the even block has off-diagonal rate
`(1-c)+c = 1`, hence is *fixed*, while the cross-copy (hidden) sector is
controlled by `c`.  For every `c in [0,1]` the generator is row-stochastic
and commutes with `R` (verified to 1e-12).

**Theorem (certified numerically).**  For `m in {2,3,5}`, `c in [0,1]`:

1. every even task (any even sources, any even readouts, any finite horizon)
   has response data whose task spectrum is independent of `c` to machine
   precision (max relative difference over the table: 7e-16);
2. the full chain's slowest nonzero mode (spectral gap) is a nontrivial
   function of `c` (table below), so the long-time behaviour of the full
   chain changes with the hidden sector even though the task does not;
3. at `c = 0` the two copies decouple and the chain acquires an extra zero
   eigenvalue (ergodicity loss in the hidden/cross-copy direction), while
   every even task is *still identical* to its value at all other `c`.

Hence: **a finite even task cannot even detect where the hidden sector's
critical point sits.**  The parameter `c* = 0` (ergodicity loss of the full
chain) can be moved or removed by re-defining the cross-copy sector without
changing any even task datum.

**Numbers** (`results/probe-c7-no-go/latest.json`).  Even-task spectra agree
to ~1e-16 across `c in {0, .25, .5, .75, 1}` at every `m in {2,3,5}`.
Full-chain spectral gaps at those `c` values: m=2: 2.0, 0.5, 1.0, 0.5, 2.0;
m=3: 3.0, 1.0, 2.0, 1.5, 1.0; m=5: 5.0, 2.0, 4.0, 3.5, 3.0.  Zero-eigenvalue
multiplicity (loss of irreducibility in the hidden direction): 2 at `c = 0`
(for m = 2, 3), 1 elsewhere — the hidden-sector ergodicity loss sits at
`c* = 0`, invisible to every even task.

## 2. Reading for the repository

* **Kalman-decomposition flip side.**  #600/#603 used the even/odd split
  (C2-equivariant Kalman) as a *positive* fact: quotient balanced data equals
  full-space data.  The same split is a *negative* fact for identifiability:
  everything the task can know is in the even sector, so any claim of the form
  "bounded finite-horizon task order constrains an infinite-volume threshold"
  must prove that the threshold lives in the task's own sector.  This is
  exactly #594's no-go target, now with a minimal exact certificate.
* **Where the positive bridge would have to go.**  A positive threshold
  statement needs one of: critical-sector completeness (the singularity is in
  the even/task sector), horizon growing with the correlation length,
  uniform spectral approximation, or an explicit map from task data to the
  singularity.  The two-copy family shows what happens when none of these is
  imposed: the gap is moved by an invisible parameter.
* **Relation to the probe's other no-go material.**  Round-1 note
  (`probe-threshold-identifiability-and-smooth-null`) gave the *unbranched-
  vs-fork* certificate (same survival law, different fork behaviour).  This
  note gives the *sector-surgery* certificate (same task, different
  hidden-sector dynamics).  They are the two complementary ways a bounded
  task description can fail to identify the object behind it.

## 3. Boundary

* Exact finite statement about the two-copy construction; it does **not**
  claim that site percolation's `p_c` is invisible to the repository's tasks
  — it only fixes the burden of proof: a task that is a functional of a
  sector that does not contain the singularity cannot constrain it.
* `c` here is an abstract coupling, not a percolation probability; no
  continuum statement is made.

Files: `scripts/probe/c7_no_go.py`,
`results/probe-c7-no-go/latest.json`.
