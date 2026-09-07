# Probe PROGRESS (one page)

## Round 0–2 (completed)
- Independent P398 reimplementation; Gate-1 reproduced (orbit == lumping); Phase C quotient-balanced factorization verified to 1e-16 + inner-product diagnosis 37–61%; finite-group selection rule (order-ell) with C3 demo; predictive rank lemma (affine response rank <=2 vs classes k+1); no-go scope; atlas slice.
- Deliverables under `/workspace/matching-one-round1/notes/` (probe-*.md).

## C1 (completed, 2026-09-06)
- Snapshot v2 written: sibling front #600/#601/#602/#603 absorbed; no conflicts; round-1 confirmed by #603; round-2 positioned as the #601 Q3 gap-filler (markov-generator pointwise form), Theorem 2(b) = Kalman/Wonham (cited); experiment-language complexity territory still empty.

## C2 (completed)
- Protocol reconstructed from the #549 branch verifier; closed form re-verified exactly k=1..12.
- Proposition (proved): single-root fork language is affine in a => response rank <= 2 regardless of clone width.
- Abstract exact d-root model: degree(P)=d; tests d=0..k give exact rank k+1; h_A=h_B degenerates to rank 1.
- Answers #599 Program B Q1/Q2/Q3 at protocol level: rank growth requires language DEPTH (independent roots), not clone width; depth d caps rank at d+1; separating k+1 classes needs depth k.
- Cheapest lift to real N16: two-root protocol on actual #435 pair (documented, not run).

## C4 (completed)
- Formalised experiment languages (static/test; kappa_L, r_L), minimal distinguishing depth d*.
- New P398 numbers: D0/D1/D2 static signature classes at w4-8 (e.g. w8: 32/156/209 vs 750 orbits vs 1430 states) -> full declared dictionary never resolves C2 orbits at w>=5; at w=4 D2=10=orbits.
- Clean hierarchy 1430 > 750 > 209 > 32 > 3-4, each term attached to a distinct formal object.
- d*(#549 k+1 classes) = k (single-root cap rank 2, depth-k reaches rank k+1).
- Next: dynamical kappa for D2 on P398; compositional-depth theorem; N16 two-root lift.

## C3 (completed, acquisition)
- Cost probe (numpy sparse): w9/w10 fully feasible (one e^{tG}v = 0.06 s / 0.24 s; nnz ~11-12/row).  Repo "pure-Python unreachable" is a code cap (width<=8), not a barrier.
- Real gate = porting the frozen #580/#588 projected-memory construction + continuity gate.  Decision table written for the w9 saturation run.  (probe-w9-memory-acquisition-20260906.md)

## C5 (completed)
- Formal Markov-generator selection-rule statement (citation-disciplined; fills the #601 Q3 gap: Schur/Kalman/Wonham cited as prior, claimed content = Markov/Duhamel pointwise form + order-ell tensor criterion + C2-accident).
- Numerics #603 did not run: (a) P398 L2(pi) full-vs-quotient factorization holds at all widths 4-8 to ~1e-16 (spectra convention-dependent: w8 lead sv 0.1197 vs 0.6406 counting; identity convention-robust); odd-width teeth reproduce under the second product.  (b) C3 isotypic-fibre realisation: W-task full vs fibre 1.4e-15; joint spectrum = union of sector spectra.  (probe-markov-selection-rule-and-l2pi-fibre-20260906.md)

## C6 (deferred)
- #582 discriminator empirical run needs the archived Q_N data / owner authorisation; not run to avoid document-without-data.  Discriminator protocol + decision criteria already live in round-1 threshold note (Part 2) and are cited in the ladder.

## C7 (completed)
- Exact threshold no-go on a two-copy C2 family: every even task datum independent of g to ~1e-16 while full-chain gap = g*m -> 0 at g=0 (ergodicity loss in the hidden sector).  A finite even task cannot constrain where a hidden-sector critical point sits.  (probe-threshold-no-go-theorem-20260906.md, scripts/probe/c7_no_go.py)

## C8 (completed)
- Invariant/observer-change catalogue over all probe material: selection zeros, quotient identity and per-sector inertness are the robust invariants; all positive magnitudes are dictionary/inner-product/language dependent.  (probe-invariant-catalogue-20260906.md)

## C9 (completed)
- Complexity-ladder theorems v1: every arrow labelled P/F/C/U with the counterexample archive (X1..X9) and ranked unknown arrows with cheapest falsifications.  (probe-complexity-ladder-v1-20260906.md)

## C10/C11/C12 (completed)
- C10: three cheap tests executed (T1: even-language dynamic resolution saturates at the orbit quotient, w4-8; T2: marked odd channel sees the C7 hidden coupling, even tasks do not; T3 = C2 protocol algebra).  Ranked bridge conjectures 1-8 with cheapest falsifications.  (probe-bridge-conjectures-ranked-20260906.md, scripts/probe/c10_cheap_tests.py)
- C11: full submission pack = submit-pack-v2/ (branch research/probe-full-round1-20260906) + matching-one-probe-full-20260906.bundle; open PR, no auto-merge.
- C12: mission update - north-star proposition supported with sharpened boundaries; durable deliverables listed; continuation = r_lin(w8) definition, w9 memory port, two-root N16, then C6.

## Open items ledger
- w8 exact r_linear (blocked on #593 definition).
- r_positive(9/10)=2494/8524 falsification run (needs w9/10 compute; #593).
- #582 discriminator run (planned C6).
