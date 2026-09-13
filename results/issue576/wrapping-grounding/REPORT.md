# Issue #576 — wrapping channel on published ground

See `notes/issue576-wrapping-grounding.md` for the full report (Parts 1–3).

## Part 1 — Akhunzhanov–Eserkepov–Tarasevich torus wrapping polynomials
- Retrieved `torus.txt` from arXiv ancillary (https://arxiv.org/src/2204.01517v1/anc/torus.txt).
- Independent exact enumeration (`scripts/exact_wrapping_enum.cpp`) for L=2..5:
  L=2: 7/16; L=3: 175/512; L=4: 19571/65536; L=5: 8853291/33554432.
- Validated published coefficients for L=3 (175) and L=4 (19571) exactly.
- L=5 discrepancy: published Σ c_k = 8853301 vs independent 8853291 (diff 10) — flagged, unresolved.
- L≥9 (and malformed L=10) blocks in torus.txt give P_L(1/2) below the continuum limit → appear corrupted; not used as ground truth.
- Main (NN square site) channel converges to the repo's committed Pinson/Arguin continuum baseline π_i({1,0}) = 0.1694154… .

## Part 2 — Pinson/Arguin continuum wrapping at r=1,2,4
- π({1,0})(i)   = 0.16941543532134688938260796919875445000145337645375
- π({1,0})(2i)  = 0.50303589769463904028462207879367242007166407775608
- π({1,0})(4i)  = 0.85596932105387227330836216781610004621209130143655
- Two independent evaluations (direct Gaussian sum vs theta form) agree to 1e-82…1e-90.
- Decision: (ii) explicit non-claim — matching-odd slope is NOT claimed to be Pinson wrapping.

## Part 3 — three distinct 11/4 objects
- Ê4(2i)/Ê4(i) = 11/4 (modular weight-4 amplitude ratio)
- Newman–Ziff L^{-11/4} (finite-size convergence rate, 1/ν+θ = 11/4)
- Mertens–Ziff matching root ~ L^{-4} (third, different exponent)
Not conflated; no claim-ledger entry added.

CI: Full Matching-One repository CI has not been run for this commit.
