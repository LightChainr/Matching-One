# Issue #576 — wrapping channel on published ground

See `notes/issue576-wrapping-grounding.md` for the full report (Parts 1–3).
This REPORT reflects the 2026-09-13 resolution pass (compute1-retry).

## Part 1 — Akhunzhanov–Eserkepov–Tarasevich torus wrapping polynomials
- Retrieved `torus.txt` from arXiv ancillary (https://arxiv.org/src/2204.01517v1/anc/torus.txt).
- Two independent exact enumerations at L=2..5:
  - doubled-grid criterion (`exact_wrapping_enum.cpp`): 7, 175, 19571, **8853291**;
  - displacement-union-find (`exact_wrapping_enum2.cpp`): 7, 175, 19571, **8853301**.
- **RESOLVED:** the published Σ c_k = 8853301 (L=5) is **correct**; per-coefficient
  match at L=3,4,5 against the displacement-DSU enumerator. The doubled-grid
  criterion misses exactly **10 spiral configurations** (cluster winds x AND y;
  no same-row (r,0),(r,L) pair in the doubled strip; all 10 in the c_15 bin;
  listed in `derived/spiral_discrepancy_resolution.json`).
- **RETRACTED:** the earlier "L≥9 blocks appear corrupted" claim. It compared
  P_L(1/2) against 0.1694 — a category error, since P_L(1/2) → 0 (p=1/2 < p_c).
  Correct anchor: P_L(p_c) → R^v = 0.521058290… (Mertens–Ziff). All blocks
  L=3..12 pass the paper's proven divisibility test (c_k ≡ 0 mod L²/gcd(k,L²)),
  and P_L(p_c) converges smoothly to 0.521058… (table in the note).
- Only genuine file defect: L=10 block has 100 of 101 coefficients; the missing
  trailing coefficient is necessarily c_100 = 1 (truncation, not corruption).
- The main (NN square site) channel is on published closed-form ground; the
  finite-L polynomials and the continuum anchors (R^v at p_c, Pinson/Arguin
  class probabilities) are mutually consistent. NN+NNN (Sq8) matching remains
  out of the paper's scope and unvalidated against it.

## Part 2 — Pinson/Arguin continuum wrapping at r=1,2,4
- π({1,0})(i)  = 0.16941543532134688938260796919875445000145337645375
- π({1,0})(2i) = 0.50303589769463904028462207879367242007166407775608
- π({1,0})(4i) = 0.85596932105387227330836216781610004621209130143655
- Recomputed from scratch on the resolution pass; direct vs theta forms agree to
  ~1e-110…1e-116 at 100 dps. Ratios: 2.9692447842221780796 (r=2), 5.0524872154079183417 (r=4).
- Pruessner–Moloney (cond-mat/0310361) report agreement with Pinson at r=2 to < 1e-8.
- Decision: (ii) explicit non-claim — matching-odd slope is NOT claimed to be Pinson wrapping.

## Part 3 — three distinct 11/4 objects
- Ê4(2i)/Ê4(i) = 11/4 (modular weight-4 amplitude ratio)
- Newman–Ziff L^{-11/4} (finite-size convergence rate, 1/ν+θ = 11/4)
- Mertens–Ziff matching root ~ L^{-4} (third, different exponent)
Not conflated; no claim-ledger entry added.

CI: Full Matching-One repository CI has not been run for this commit.
