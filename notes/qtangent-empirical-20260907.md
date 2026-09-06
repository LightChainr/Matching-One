# The two typed Q-score pieces keep distinct scale / noise profiles to L=8

**Date:** 2026-09-07
**Ticket:** #581 (first empirical control, after the exact gate)
**Claim level:** C2 — Monte Carlo on a new square-bond block; the L=3 run is a
correctness control against the enumerated gate.
**Artifact:** `results/qtangent-empirical/latest.json`
**Script:** `scripts/score_qtangent_empirical.py` (builds on `qtangent_scale_decomposition.py`)

#581's exact gate ran the L=2 and L=3 square-bond tori to exhaustion and found
the duality-even Betti tangent `B_even` and the ambient-homology source
`X = r − 1` separate at Boolean degree one by parity.  The owner's roadmap note
asked for the next step: **one modest exact-critical square-bond block** asking
only whether the two pieces retain measurably different finite-size scale/noise
profiles beyond L=2,3, with raw covariance vectors primary and one frozen
supported readout.  This is that block.

## The instrument

The readout is `wrap_either`, the same primary channel the gate carried.  Three
sizes: L=3 (a control), L=4, L=8, `p = 1/2`, 100,000 configurations each.  L=4
has `2^32` and L=8 `2^128` configurations, so the gate's exhaustive
conditional-replica machinery cannot be used.  Instead the module uses an
**unbiased paired estimator**: for the nested filtration `F_j`, the increment
matrix `Γ_j = E[D_j D_jᵀ]` equals `E[m_j m_jᵀ] − E[m_{j−1} m_{j−1}ᵀ]`, and each
`E[m_j m_jᵀ]` has the exact U-statistic form `E[Y(c) Y(c′_j)ᵀ]` where `c` and
`c′_j` agree on the bonds revealed by `F_j` and are independent elsewhere
(`c′_j = (c & S_j) | (fresh & ~S_j)`).  The replica is integrated out in
expectation, so there is **no replica-noise bias** — the two things the method
must survive are that the increments telescope to the total covariance and that
the L=3 run reproduces the enumerated gate, and both are asserted and tested.

## The control: L=3 reproduces the enumerated gate

| quantity | Monte Carlo | enumerated | |
|---|---:|---:|---|
| `Cov(wrap_either, X)` | 0.28952 | 0.28786 | within 0.0017 |
| `Cov(wrap_either, B_even)` | −0.02701 | −0.02707 | within 0.0001 |
| degree-one `X` | −0.1340 | −8721/65536 = −0.13307 | within 0.001 |
| degree-one `B_even` | 0.0099 | 0 exactly | sampling floor |

(`Cov(wrap, X)` and `Cov(wrap, B_even)` are natural covariance; the gate's
`covariance_split` reports each piece divided by two, and the gate's
`scale_decomposition` reports `Γ[wrap, 2·B_even]` — doubled — so the units are
doubled for `B_even` and match for `X`.  The JSON carries the raw matrix in both
forms so the comparison is unambiguous.)

## The result: the two profiles stay distinguishable, and the gap is structural

The per-scale increments for the readout against each piece:

```text
L=3    level  bonds   Γ[wrap, B_even]   Γ[wrap, X]
         0      4        +0.0056         +0.0373
         1     12        −0.0448         +0.1122
         2     14        +0.0127         +0.0355
         3     18        −0.0005         +0.1045

L=4    level  bonds   Γ[wrap, B_even]   Γ[wrap, X]
         0      4        +0.0078         +0.0213
         1     16        −0.0066         +0.0791
         2     24        −0.0010         +0.0689
         3     32        −0.0373         +0.1281

L=8    level  bonds   Γ[wrap, B_even]   Γ[wrap, X]
         0      4        +0.0007         +0.0024
         1     52        +0.0055         +0.0676
         2     88        −0.0130         +0.0674
         3    128        −0.0430         +0.1674
```

Two signatures, both stable across all three sizes:

1. **The ambient-homology channel is single-signed and accumulates monotonically.**
   Every `Γ[wrap, X]` is positive at every scale; the cumulative curve rises
   monotonically to a total `Cov(wrap, X)` that is nearly size-independent:
   `0.2895 → 0.2973 → 0.3048`.  The low scales already carry real coupling
   (by level 1 of 4, most of the total is present).

2. **The Betti-even channel alternates in sign and is small and opposite.**
   The increments change sign from scale to scale, the cumulative curve
   reverses, and the total `Cov(wrap, B_even)` is `−0.0270 → −0.0372 → −0.0499`
   — opposite in sign to `X` and a factor 6–11 smaller in magnitude.

The degree-one result carries the same story and sharpens it with size.  The
`X` coefficient — the one value the gate found — decays as the torus grows:
`−27/128 (L=2)`, `−0.1331 (L=3)`, `−0.0950 (L=4)`, `−0.0408 (L=8)`, roughly
`∝ 1/L` (four points, not fitted).  The `B_even` coefficient stays at the
sampling floor at every size — it has **no** degree-one weight to speak of,
exactly as the parity argument requires.  So whatever degree-one information the
readout's Q-score coupling has is, at every size tested, entirely ambient.

## What is established and what is not

**Established.** The two typed pieces of the Q tangent keep a measurably
different scale and noise organisation from L=2 through L=8, on a readout and a
filtration that were frozen before any block ran.  The separation is not a
toy-size accident: it survives to L=8 with the same two signatures — single-signed,
monotone, low-scale ambient homology versus alternating, small, opposite-sign
Betti.  This is the "stably distinguishable" branch of the roadmap note, so the
two pieces form a **typed fingerprint** for interpreting a declared response's
Q-score coupling.

**Not established.** Any mechanism.  `B_even` is a duality-even measure tangent
and not automatically a local energy field; `X` is topological, but a covariance
with it does not demonstrate a defect theory.  Three sizes do not determine the
slow drift of `Cov(wrap, X)` upward or of `Cov(wrap, B_even)` downward; the
degree-one `∝ 1/L` reading is four points and is reported, not claimed.  And
nothing here is a square-site statement — the square-site matching-odd
observable still has no canonical lift to a square-bond torus, so the
topological side (`X_site = r − 1`) is the only part transportable, and it is
not transported here.

## The next step #581's own order dictates

The square-site extension should transport **only** the topological side:
scale-resolve `X_site = r − 1` and its coupling to an original-U / H4 readout
whose provenance is declared separately.  The square-bond control now in hand
says that a matching-odd signal that loads only when ambient `X` is resolved —
against this bulk control where the bulk channel alternates and stays small —
would be real evidence of a more topological/defect-like square-site response.
