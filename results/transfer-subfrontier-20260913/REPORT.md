# REPORT — P11 Transfer-matrix sub-frontier reproduction and resource model

Date: 2026-09-13. Ticket: issue #11 (`T07`). Branch:
`analysis/p11-transfer-subfrontier-20260913`. All numbers in this report are
either exact, or measured, or model-based extrapolations that are explicitly
labelled as such. No width beyond what we actually computed is claimed.

## 1. Sub-frontier reproduction (Mertens 2022, arXiv:2109.12102)

Definition chain reproduced exactly: `R_{n,m}(p) = Σ_k A(k) p^k (1-p)^{nm-k}`
on the free-boundary n×m grid, spanning = one 4-neighbour cluster touching the
virtual all-occupied row 0 and row m; `p_med`: R_{n,n}=1/2; `p_cell`:
R_{n,n}=R_{n-1,n-1}.

### 1.1 Digit agreement with the published table (data/mertens_2022_square_site_estimators.csv)

| n | p_med digit match | p_cell digit match | path |
|---|---|---|---|
| 2..10 | 30/30 printed | 30/30 printed (n=2..5 checked) | exact DP (path A, python) |
| 11 | 32 (all 30 printed) | — | CRT container |
| 12 | 32 (all 30 printed) | 32 (all 30 printed) | CRT container |
| 13 | 32 (all 30 printed) | 32 (all 30 printed) | CRT container |
| 14 | 32 (all 30 printed) | 32 (all 30 printed) | CRT container |
| 15 | in flight at delivery, not claimed | in flight | CRT container |

"32" means our 60-digit mpmath root agrees with the published 30-digit string
and continues to match beyond it (32 leading digits); the published table
stops at 30.

Exact reproduced values (leading digits):
- p_med(11) = 0.587545601376707076865096376747781238324764
- p_med(12) = 0.588212470606443263171973079741831708330245
- p_cell(12) = 0.598724257102302868949743766602534710805687
- p_med(13) = 0.588753953651382767097855532073340370474632
- p_cell(13) = 0.59802106388297966577943835943963309465498
- p_med(14) = 0.589200171193723644344640059478678608017814
- p_cell(14) = 0.597439041437080848283968950089603109720644

### 1.2 Independent validation of the state encoding

- Path A (label-based connectivity-state DP, exact integer polynomials,
  `scripts/p11_square_dp.py`) vs path B (exhaustive 2^(nm) enumeration with
  union-find, `scripts/p11_bruteforce.{py,cpp}`): A(k) vectors bit-identical
  for every width both reach (n=2..5).
- Path A vs published p_med: 30/30 digits at n=2..10.
- C++ cell-sweep engine (`scripts/p11_exact_dp.cpp`) vs path A: A(k) mod
  independent 30-bit primes identical for n=4..10 (CRT mode).
- State counts equal Mertens' S_n (Motzkin triangle) exactly: 29964 (n=11),
  83304 (12), 232323 (13), 649845 (14), 1815992 (15) — the encoding
  generates precisely the paper's state space.

### 1.3 Binary128 digit audit

On aarch64 Linux `long double` is IEEE binary128. Newton root-finding in
binary128 (with derivative DP and final bracket verification at 1e-25)
reproduces the published 30 digits with margin at every audited width:

| n | f128 p_root | digits matching published | bracket_ok |
|---|---|---|---|
| 5 | 0.575810073211627653605032974314577674... | 32 | true |
| 6 | 0.579702757132443521419439978330558427... | 32 | true |
| 7 | 0.582351295080082980073474691830408397... | 32 | true |
| 8 | 0.584241466489847673860351132398167166... | 32 | true |
| 9 | 0.585641556861396511416995666354642547... | 32 | true |
| 10 | 0.586710034053406359804690473124405622... | 32 | true |

Caveat carried forward: this audit certifies ~32 correct digits only where an
exact cross-check exists (n≤10, 13 via CRT). The f128 runs at n=16..18 have
**no exact cross-check**; we report them with a bracket verification at
1e-25 but make no certified digit claim beyond what the audit supports.

## 2. Jacobsen 2015 cylinder sequence — NEGATIVE result

Our independent strip-sector implementation of the open/closed eigenvalue
identity gives p_c(1)=1/2 (matches Table 2) but p_c(2)=0.52002090625...
vs published 0.5651977173836393964375280132470308160984. Diagnosis: the
paper's sectors are blocks on the periodic Temperley–Lieb s=0 reduced states
(two directions glued), not strip-frontier connectivity states. We claim
**no** cylinder width beyond n=1. Details in
`notes/p11-transfer-subfrontier-20260913.md` §1.3 and
`scripts/p11_cylinder_sectors.py`.

## 3. Measured scaling and resource model for n=25..28

Measured on ARM containers (16 vCPU / 32 GiB cgroup limit), exact CRT mode:

| n | peak states (map size) | peak RSS (KiB) | wall s (primes) | wall s/pass |
|---|---|---|---|---|
| 11 | 29964 | 62992 | 19.9 (10) | 4.0 |
| 12 | 83304 | 194920 | 89.0 (10) | 17.8 |
| 13 | 232323 | 609484 | 394.4 (10) | 65.7 |
| 14 | 649845 | 1953088 | 1980.3 (7) | 282.9 |

Log-linear fits y = a·b^n over measured widths (max relative residual):

- states: b = 2.7887, resid 1.6e-3 (cross-check: exact S_{n+1}/S_n ≈ 2.80)
- peak RSS: b = 3.1401, resid 1.1e-2
- wall per pass: b = 4.0982, resid 6.0e-2
- wall all primes: b = 4.6169, resid 5.0e-2 (not width-scalable — prime
  counts differ per width; reported for completeness)

Extrapolation to n=25..28 (**model-based, not measurements**; the fit uses
only n=11..14 and residuals are small, but widths 25..28 are 11+ widths
beyond the fitted interval, so treat as order-of-magnitude):

| n | peak states | peak RSS | wall per pass |
|---|---|---|---|
| 25 | 5.1e10 | 566 TiB | 1.5e9 s (~50 y) |
| 26 | 1.4e11 | 1.78 PiB | 6.3e9 s |
| 27 | 4.0e11 | 5.58 PiB | 2.6e10 s |
| 28 | 1.1e12 | 17.5 PiB | 1.1e11 s |

Conclusion: exact CRT to n=25..28 is unreachable on 32 GiB containers by 4–5
orders of magnitude in memory alone; even n=16 (extrapolated ~20 GB, i.e.
S_16=5.08e6 states × 17-coefficient big vectors) is at the edge. This is why
the deliverable is a sub-frontier reproduction plus this resource model, not
widths 25–28.

Bounded f128 observation: a binary128 evaluation run at n=18 was observed at
peak RSS 7.75 GB and 617 s elapsed while still running
(`raw/probe_f128_med_n18_partial.json`); f128 RSS at n≤10 is ~15 MB. This is
a partial record of a still-running process, not a completed measurement.

## 4. Checkpoint / restart

`p11_exact_dp.cpp crt --dump-row R --map-out F` writes the row-R state map in
a fixed binary layout (sorted keys); its SHA-256 is a deterministic
state-enumeration hash. Verified: three consecutive runs byte-identical
locally (n=8, row 5); container XPk2PZ n=14 dump row 7 (`ck_n14.p0`,
SHA-256 `8ec26ebf6c42fe07334025ddb5e77e1cfc79164930b7332df6571630806cc4f3`)
+ restart produced prime-0 A(k) **bit-identical** to the uninterrupted run
(`raw/crt_n14_restart.json`, machine-checked equality).

## 5. What is NOT covered (read before citing)

- Cylinder p_c(n) n≥2 (Jacobsen) — negative result, §2.
- Exact CRT beyond n=15 — memory-extrapolated beyond containers.
- f128 n=16..18 — computed, but digit accuracy there is not certified by an
  exact cross-check (see §1.3 caveat).
- Widths 19..24 — not attempted; resource model says out of budget.
- Widths 25..28 — model-based extrapolation only, §3.

## 6. Artifact index

- `raw/roots_n{11,12,13,14}.json` — exact CRT roots + digit match report
- `raw/crt_n{11..14}.json`, `raw/crt_rows_n{11..14}.rows` — container outputs + per-row telemetry
- `raw/crt_n14_parsed.json`, `raw/crt_n14_restart.json`, `raw/ck_n14.sha256` — checkpoint/restart evidence
- `raw/probe_crt_n{11..14}.json` — probe timing/RSS records
- `raw/f128audit_med_n{5..10}.json` (+probes) — binary128 digit audit
- `raw/probe_f128_med_n18_partial.json` — bounded partial observation
- `raw/resource_model.json` — machine-readable fits + extrapolation
- `commands.txt` — exact command lines for every artifact
- `metadata.json` — machine, toolchain, versions
- `raw/resource_model.json` — machine-readable fits + extrapolation
- `commands.txt` — exact command lines for every artifact
- `metadata.json` — machine, toolchain, versions
