# P11 — Transfer-matrix sub-frontier reproduction and resource model (2026-09-13)

Ticket: Issue #11 / `T07` (reproduce cylindrical/helical transfer-matrix
sequences through the published frontier, then build an honest resource model
for `n=25..28`). Branch: `analysis/p11-transfer-subfrontier-20260913`.

## 1. What was reproduced

### 1.1 Mertens 2022 exact spanning probabilities (arXiv:2109.12102)

Definitions verified against the paper: `R_{n,m}(p) = sum_k A_{n,m}(k) p^k
(1-p)^{nm-k}` on the free-boundary `n x m` grid, spanning = one occupied
4-neighbour cluster touching the virtual all-occupied row 0 and row `m`;
`p_med`: `R_{n,n} = 1/2` (paper Eq. 14a); `p_cell`: `R_{n,n} = R_{n-1,n-1}`
(Eq. 14b).

`A_{n,m}(k)` are exact integers computed by transfer-matrix DP carrying
integer polynomials; roots computed with mpmath at 60 digits.

Two independent paths validate the state encoding (acceptance item):

| path | implementation | widths |
|---|---|---|
| A | label-based connectivity-state DP, integer polynomials (`scripts/p11_square_dp.py`) | 2..10 |
| B | exhaustive `2^(nm)` enumeration with union-find spanning test (`scripts/p11_bruteforce.{py,cpp}`) | 2..5 |

`A(k)` vectors agree **bit-for-bit** between paths A and B on every width both
can reach (n=2..5), and path A reproduces the published `p_med(n)` to all 30
printed decimals for n=2..10 (digit-by-digit against
`data/mertens_2022_square_site_estimators.csv`; p_cell checked n=2..5).
A third implementation (C++ cell-sweep engine, `scripts/p11_exact_dp.cpp`,
CRT mode) reproduces path A's `A(k)` mod independent 30-bit primes for
n=4..10 after two bugs were found and fixed (see §4).

State-space cross-check: our full-row state counts equal Mertens' `S_n`
(second column of the Motzkin triangle) — 9, 25, 69, 189, 518, ... at
n=3..8, and 29964, 83304, 232323, 649845 at n=11..14 — i.e. our
label-based encoding generates exactly the same state space as the paper's
parenthesis-signature encoding, independently.

### 1.2 Container runs (exact, CRT) — n=11..15

(Numbers from the 10 ARM containers, telemetry via
`scripts/transfer_resource_probe.py`; raw files under
`results/transfer-subfrontier-20260913/raw/`.)

- n=11: p_med reproduced to 30/30 published decimals (ours: 32 leading
  digits, 0.587545601376707076865096376747781238324764). Peak RSS 62992 KiB,
  19.9 s for all 10 CRT primes.
- n=12: p_med 0.588212470606443263171973079741831708330245 and p_cell
  0.598724257102302868949743766602534710805687 — 32/30 digits each.
  Peak RSS 194920 KiB, 89.0 s.
- n=13: p_med 0.588753953651382767097855532073340370474632 and p_cell
  0.59802106388297966577943835943963309465498 — 32/30 digits each.
  Peak RSS 609484 KiB, 394.4 s.
- n=14: p_med 0.589200171193723644344640059478678608017814 and p_cell
  0.597439041437080848283968950089603109720644 — 32/30 digits each.
  Peak RSS 1953088 KiB (~1.95 GB), 1980.3 s for 7 CRT primes
  (282.9 s/pass). Deterministic checkpoint/restart verified on the same
  container: row-7 dump (`ck_n14.p0`, SHA-256 8ec26ebf…cc4f3) restarted and
  reproduced prime-0 A(k) **bit-identically** (`raw/crt_n14_restart.json`).
- n=15: CRT exact run on ZyTrST, in flight at delivery time (~1.8e6 states
  at row 8 of pass 1; total runtime estimate several hours). Reported as
  in-flight, not as a result; its telemetry is not included in the fits.

Binary128 (long double on aarch64 = IEEE binary128) Newton root-finding
cross-checks the exact roots to ~31-32 significant digits at every width
where both run, with final bracket verification `|f(p±10^-30)|` sign change.

### 1.3 Jacobsen 2015 cylinder eigenvalue identity (arXiv:1507.03027) —
###     attempted, NEGATIVE result, not reproduced beyond n=1

The published cylinder sequence `p_c(n)` (Table 2, 40 digits, n<=21) is
defined by the eigenvalue identity `Lambda_open = Lambda_closed` between the
two topological sectors of the **periodic Temperley-Lieb** transfer matrix.
We implemented an independent strip-sector formulation (closed sector =
no transverse wrap ever; open sector = a wrapping cluster that keeps frontier
presence forever; `scripts/p11_cylinder_sectors.py`):

- n=1: `p_c = 1/2` exactly — matches Table 2 row n=1.
- n=2: our formulation gives `0.520020906250797790345118234280788...`
  vs published `0.5651977173836393964375280132470308160984`. **MISMATCH.**

Diagnosis: the open/closed sectors of arXiv:1507.03027 are blocks of the
transfer matrix on the s=0 *reduced states of the periodic TL algebra*
(both lattice directions glued; states are annular link patterns with a
nontrivial gluing), not connectivity states of a strip frontier. The naive
sector split is not equivalent; reproducing Table 2 faithfully requires
reimplementing that machinery (est. multi-day effort). Per the dispatch
rule ("if your independent check contradicts your main result, that is the
deliverable") we report the mismatch with the full 3x3 sector matrices in the
script, and we do NOT claim any cylinder width beyond n=1.

## 2. Measured scaling and the resource model for n=25..28

All fits are log-linear `y = a·b^n` over **measured widths only**
(`scripts/p11_resource_model.py`; machine-readable output in
`results/transfer-subfrontier-20260913/raw/resource_model.json`):

| quantity | points (n: value) | growth b | max rel resid |
|---|---|---|---|
| peak map states | 11: 29964, 12: 83304, 13: 232323, 14: 649845 | 2.7887 | 1.6e-3 |
| peak RSS (KiB) | 11: 62992, 12: 194920, 13: 609484, 14: 1953088 | 3.1401 | 1.1e-2 |
| wall s/pass | 11: 4.0, 12: 17.8, 13: 65.7, 14: 282.9 | 4.0982 | 6.0e-2 |

The state growth base 2.7887 matches the exact S_{n+1}/S_n ≈ 2.80 ratio
independently — a consistency check on the telemetry. (Widths used
different prime counts — 10 primes at n≤13, 7 at n=14 — so the per-pass
fit is the width-scalable one; the "all primes" fit is also reported in
raw/resource_model.json.)

Extrapolation to n=25..28 (**model-based, not measurements**; the fitted
interval ends at n=14, so treat these as order-of-magnitude):

| n | peak states | peak RSS | wall per CRT pass |
|---|---|---|---|
| 25 | 5.1e10 | 566 TiB | 1.5e9 s (~50 y) |
| 26 | 1.4e11 | 1.78 PiB | 6.3e9 s |
| 27 | 4.0e11 | 5.58 PiB | 2.6e10 s |
| 28 | 1.1e12 | 17.5 PiB | 1.1e11 s |

Even n=16 exact CRT extrapolates to ~20 GB (S_16 = 5080510 states with
17-coefficient big-integer vectors), at the edge of 32 GiB. The honest
conclusion: **widths 25–28 are 4–5 orders of magnitude out of reach in
memory alone on 10×32GiB containers**; the deliverable for those widths is
this resource model, not numbers.

Bounded f128 observation: a binary128 run at n=18 was observed at peak RSS
7.75 GB / 617 s elapsed while still running
(`raw/probe_f128_med_n18_partial.json`); at n≤10 f128 peak RSS is ~15 MB
(flat in n there — fixed overhead dominates at small widths).

## 3. Checkpoint / restart (deterministic)

`p11_exact_dp.cpp crt --dump-row R --map-out F` writes the row-R map
(sorted keys, fixed binary layout) — its SHA-256 is a deterministic
state-enumeration hash; three consecutive runs produced byte-identical files.
`crt --map-in F` restarts from the checkpoint. Verified locally at n=8
(dump row 5, restart, identical `A(k) mod p`) and on container XPk2PZ at
n=14 (dump row 7; restart continues rows 8..14; prime-0 `A(k)` identical to
the uninterrupted run; SHA-256 recorded in `raw/ck_n14.sha256`).

## 4. Bugs found on the way (each would have silently produced wrong numbers)

1. Path A: non-top clusters whose last frontier site ends were wrongly
   dropped (they are merely finished clusters; only the top cluster may
   never be severed). Found by brute-force mismatch at n=5, k=11 (4 configs).
2. Path A: union-find in the pattern successor did not encode old-row label
   connectivity, splitting a cluster whose old sites bonded at two positions.
3. C++ sweep: emptying a cell can leave a gap in canonical label numbering;
   equivalent states then fail to merge (first fires at n=6).
4. C++ sweep: merging two non-top classes relabelled only old-row cells,
   splitting classes that already owned new-row cells (first fires at n=6).
5. CRT prime list contained composites sharing a factor 7 (1073741837 =
   7x153391691, 1073741963 = 7x153391709) — CRT reconstruction silently
   broke. All moduli are now sympy-verified primes, and the combiner asserts
   pairwise coprimality.

Every fix is pinned by the cross-validation matrix of §1.1.

## 5. What we could NOT do

- Jacobsen cylinder n>=2 (see §1.3).
- Exact CRT beyond n=15 in 32 GiB: peak map ~ `0.35*2.85^15` states x
  (n*m+1)-coefficient vectors is already multi-GB; n=16 CRT is ~20 GB and
  was not attempted. n=16..18 use binary128 evaluation instead (digit claims
  audited against exact results at n<=10; binary128 loses ~2 digits at the
  root, bracket verification at 1e-30 included in raw/).
- Widths 19..24 were not attempted: extrapolated cost exceeds what 10x32GiB
  containers can deliver in this ticket's budget (see the resource model).
  Nothing in this report should be read as covering n>=19.

## 6. Files

- `scripts/p11_square_dp.py` — path A exact DP
- `scripts/p11_bruteforce.{py,cpp}` — path B brute force
- `scripts/p11_exact_dp.cpp` — C++ sweep engine (crt / f128 / res modes)
- `scripts/p11_crt_roots.py` — CRT combine + roots + digit comparison
- `scripts/p11_cylinder_sectors.py` — cylinder attempt (negative result)
- `scripts/p11_resource_model.py` — fits + n=25..28 extrapolation
- `results/transfer-subfrontier-20260913/raw/` — all telemetry and outputs
- `results/transfer-subfrontier-20260913/REPORT.md` — summary tables
