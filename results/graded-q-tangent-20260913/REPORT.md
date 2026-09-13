# #586 Graded Q→1 tangent benchmark — REPORT

**Issue:** #586 (graded tangent benchmark: separate Betti-even and ambient-homology response)
**Branch:** `theory/p586-graded-q-tangent-20260913`
**Method script:** `scripts/graded_q_tangent_benchmark.py`
**Reused exact gate:** `scripts/qtangent_scale_decomposition.py` (#581)

## Headline

- **Phase A continuum internal check passes:** the Q-tangent of a representative
  (2,1)-degenerate BPZ ODE satisfies the differentiated inhomogeneous ODE with
  residual `max|R|/||G|| = 4.4e-6` (threshold 1e-3).
- **Phase A/B lattice control exact:** `Cov(O,T) = ½Cov(O,B_even) + ½Cov(O,X)`
  holds **exactly** (`split_is_exact = True`) for all six observables at L=2 (256
  configs) and L=3 (262144 configs); identity gate `T − T* = X` has **0 failures**.
- `open_edges` loads **purely on X** (B_even = 0); `wrap_either`/`wrap_cross`
  carry **opposite-sign B_even** pieces (∓427/65536) with identical X — confirming
  B_even and X are physically distinct, not aliased.
- **Phase B bulk-logarithmic numeric declared a buy-back** (VJS / Camia-Feng
  generic-Q data not in tree). Spin-4 tangent reopening of #263 **not yet earned**.

## Two measure-score contributions reported separately (never recombined first)

| observable (L=2) | Cov(O,T) | B_even | X |
|---|---|---|---|
| open_edges | 27/64 | 0 | 27/64 |
| wrap_either | 8405/65536 | −427/65536 | 69/512 |
| wrap_cross | 9259/65536 | +427/65536 | 69/512 |
| components | 2911/65536 | 7903/65536 | −39/512 |
| cycle_rank | 30559/65536 | 7903/65536 | 177/512 |
| wrap_direction_0 | 4335/32768 | −81/32768 | 69/512 |

(Identical exact structure at L=3; see `derived/graded_q_tangent.json`.)

## Decision

Promote the graded tangent machinery as an **exact method-level diagnostic** for
#581's scale tomography. The narrow spin-4 tangent reopening of #263 remains
**conditional on the Phase B bulk-log positive control**, which is a declared
buy-back. No result is transported to square site.

Full Matching-One repository CI has not been run for this commit.
