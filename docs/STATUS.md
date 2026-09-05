# Project Status and Claim Ledger

**Status date:** 2026-09-05

`main` is the shared research line. Claim strength follows evidence and chronology, not PR state. `docs/ROADMAP.md` ranks information gain; it is not a permission system.

## What is required while exploring

`GOVERNANCE.md` §2 is the whole list, and it is short: don't fool yourself about a
number, don't destroy data, don't misdate a freeze, say which observable, and count
one random block once. Nothing else in this repository gates exploratory work.

The three that bear directly on this ledger: frozen predictions and result history are
preserved rather than overwritten; a claim-bearing comparison uses identical observable
semantics or an exact registered map; and correlated views of one raw random block are
counted once, not as independent primary evidence.

Everything else — digests, provenance chains, second implementations, preregistration,
power — is publication-time work and lives in `docs/PUBLICATION-CHECKLIST.md`.

## Strongest current evidence

| Statement | Level | Current evidence |
|---|---:|---|
| Square-site matching-odd orientation signal exists | C3 | Independent P43+P57 primary synthesis rejects global zero: `chi2=31.1857355515/4`, `p=2.81e-6`; fixed H4 predictions give `3.4622795373/4`, `p=.484` |
| Central square-site odd sector is compatible with `DeltaCos4*N^-13/8` | C3 | P31/P32/P37/P43/P50/P57 |
| Frozen norm-5 H4 transfer beats H12/H8 aliases | C3 | H4 `0.4163/2`; H12 `35.1931/2`; H8 `16.0120/2` |
| P57 child block alone rejects zero | negative refinement | No: zero `1.77635/2`; its value is harmonic/transfer discrimination |
| N145->290 full curve is one scalar multiplier | C3 negative | No: three-level transfer `9.3520/2`, `p=.0093`; the resolved shape mode fails while the common amplitude direction remains viable |
| Frozen finite-size center-slope correction predicts N290 | C3 | corrected slope residual `z=-0.666`; bare `2^(3/8)` gives `z=-22.690` |
| Pure `P4[S'] ~ N^-5/4` is sufficient | C3 negative | Prospectively falsified; `52.71634/2` on P48 new geometry |
| One scalar width explains the higher thermal jet | C2 negative | No: full covariance norm-5 width score `24.5004/10`; width-corrected q2 `22.2386/10` |
| Rank-2/Jordan is uniquely established | C2 | No. Jordan/log is compatible (`17.0513/10`) and now has a precise Q4-module origin, but scale-log behavior alone is not module identification |
| Intrinsic quantile-center transfer obeys `N^-3/4` on N145->290 | C3 | frozen ratio observed `0.59584549` vs `2^-3/4=0.59460356`, `z=-1.033` |
| Square/rectangular spin-4 ratio carries the area-normalized weight-4 shape | C3 negative | Prospectively falsified at N=290: ratio `1.880 +/- 0.177` excludes `11/4` at `4.9` sigma, area scaling `4` at `12.0`, and no-dependence `1` at `5.0`. Tests the fingerprint as constructed, not the module: the normalization assumption is a second conjunct (`docs/astra/Q2`) |

## Exact semantics and controls

The Issue #43 even-sector channel correction remains

```text
DeltaS_cross = -DeltaS_either
```

with corrected score `0.5700315436/2`.

Finite Russo/chain rule is exact:

```text
M'(p) = pivotal_mass_primal(p) + pivotal_mass_matching(1-p).
```

The N=26 frozen finite families remain falsified:

```text
Beta(5,5): first k=5 difference = -96
Beta(7,7): first k=5 difference = +156
```

## Square-site thermal spin-4 sector

The durable empirical picture is a leading matching-odd H4-like sector with `x=21/4` candidate scaling, plus non-scalar finite-size mixing in the derivative/full-curve state.

The exact LCFT bridge is now sharper: the percolation energy Jordan pair can be lifted by the repository Q4 descendant to a rank-2 `x=21/4`, spin-4 pair. In the repository normalization `<Q4|Q4>=4930`. The resulting logarithmic slope has the exact module relation

```text
B_logN(tau) = -(lambda_top/2) A_q(tau),
A_q/A_epsilon = (493/96) g2(tau),
```

so the frozen module coefficient is `-493/192`. This supplies a representation-theory origin for Jordan/log scaling; it does **not** prove that the lattice `P4[S']` overlaps that module.

The next identifying evidence must use shape/modulus information. Exact assets include the rectangular/CM `11/4` ratio and the hexagonal degree-2 E4 phase projector. A scalar-cancelled modulus fingerprint is more identifying than another radial exponent fit.

## Pivotal and self-matching mechanism

Pivotal normalization gives two stable archived relations:

```text
N * P4[D']/Mbar'                 chi2 = 8.793/7
[P4[S']/Mbar']/P4[D]             chi2 = 9.458/7
```

while `N^(13/8)P4[S']/Mbar'` is nonconstant (`117.880/7`). A genuinely local landing-marked pivotal H4 observable is measurably orientation-sensitive.

Microscopically the N=10 local odd tangent has two independent response rows,

```text
[[15/8, 5/4],
 [-3/64, 11/64]],
```

but at N130/N170 the second singular direction is unresolved (condition numbers about 1687 and 608). More samples of the same two rows are therefore low information. The multiradius N130/N170 prototype also rejects a simple constant shell-log story and shows `R=8` is geometrically non-injective there. Future local tomography should change geometry/readout, not merely add replicas.

## Rank gap versus local thermal jet

The exact neutral-area covector maps the full Krawtchouk expansion to `E[K_plus-K_minus]/(N+1)`, but the expansion is severely ill-conditioned when truncated around the intrinsic center. The global rank gap is therefore **not** a redundant low-order thermal-jet coordinate. This closes the scalar-width/common-state shortcut.

## Distinct primitive square-bond spin-4 sector

Primitive homology characters form a separate mechanism from the square-site thermal Q4 candidate.

The continuum-subtracted non-scalar C3 character is directly observed with a passing reflection null. The simple scalar `C proportional E4(tau)` Pell phase bridge fails, so this sector should not be identified with the thermal Q4 one-point function.

Two prospectively frozen norm-2 generations instead select the negative rank-4 H4 phase:

```text
first generation:  H4 -1/2  chi2=5.4171/2
second generation: H4 -1/2  chi2=1.7077/2, p=.426
```

Frozen positive-phase alternatives are strongly excluded in both generations. Individual lineages do not converge monotonically to `-1/2`, so monotonic convergence is not claimed.

A zero-new-compute vacuum-KdV calculation predicts `C30/C56=1.99068780`; observed is `1.99360564`, with the C-only score essentially exact. The scalar S residual is a separate direction. Current interpretation: a distinct `x≈4`, spin-4 identity/vacuum-family response with finite-size corrections.

## Norm-4 quotient structure

The general integer-period backend is production-ready. Exact Gaussian-cover arithmetic shows norm-4 `2i` has deck group `Z2 x Z2`, and `(1+i)^2=2i` has an exact coarse/detail Hadamard character decomposition. There is no Gaussian scalar norm-4 cyclic `Z4` comparator. Therefore quotient dependence should be tested with character-resolved readouts when cheap, not treated as an unspecified nuisance.

## Current interpretation

1. **Signal existence is no longer the bottleneck.** Independent primary square-site blocks strongly reject global zero while remaining compatible with fixed H4 predictions.
2. **The central/derivative state is not scalar.** N290 shape, norm-5 thermal jet, rank-gap and width analyses all point toward compact mixing/transfer rather than another free correction exponent.
3. **Jordan has a concrete module origin but is not identified by scale behavior alone.** Modulus/shape is the next orthogonal discriminator.
4. **Local pivotal physics is real, but current N130/N170 readouts are nearly rank-one.** Change the readout/geometry rather than buying more of the same samples.
5. **Primitive square-bond H4 is a separate `x≈4` sector.** Do not fold it into the thermal `x=21/4` story.
6. **Norm-4 has exact deck-character structure.** Use it to sharpen, not delay, the existing production design.

## Explicit non-claims

The project does not claim a closed form for square-site `p_c`, global uniqueness of H4 or `13/8`, a unique q2/Jordan mechanism, a scalar-width explanation, proof of the lattice-to-Q4 overlap, a full matching/OPE automorphism, or a rigorous new percolation bound.
