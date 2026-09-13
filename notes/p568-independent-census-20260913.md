# Independent census of C(≤4, ≤100) against the four frozen P2 intervals

Issue #568. Everything here was written from the protocol only: the four interval
strings of `analysis/pslq_search_contract.json`, and the screening statement of
`docs/manuscripts/p2-algebraic-exclusion/manuscript.md` §3.2 (Theorem 2) with its
class definition in §3.1.

## Disclosure, first

The three primary files named in the ticket —
`scripts/degree4_interval_exclusion.py`, `scripts/degree4_fixed_point_screen.cpp`,
`scripts/exact_polynomial_root_certificate.py` — were **never read**. No share of
their code, structure or arithmetic is reused below.

This is nevertheless **not a blind replication**. The primary's committed summary
artefacts `results/pslq-degree4-*/latest.json` were read before the replication was
written, because they define the comparison quantities (per-interval retention,
the closest member, the residual convention). That is recorded here rather than
passed off as something it is not. What follows is an independent *code path*
re-deriving the numbers from the protocol and then checked cell by cell against
the committed ones.

## What was implemented

`scripts/independent_census_screen.cpp` — the screen and the enumeration, in C++.
`scripts/independent_census.py` — weights, bound, exact decisions, artefacts.
`scripts/independent_census_agreement.py` — the cell-by-cell comparison.

### The screen, re-derived

With `S` the scale, `m` the interval midpoint, `w_k = round(S m^k)` and `c_0 = S`,
put `T(a) = Σ_{k=0..d} c_k a_k`. Then `|T(a) − S·P_a(m)| ≤ H·ρ` with
`ρ = Σ_k |S m^k − w_k|`, and if `P_a` has a root in `[l,u]` the mean value theorem
gives `|P_a(m)| ≤ D(u−l)/2` with `D = H·d(d+1)/2`. Hence every root-carrying `a`
obeys `|T(a)| ≤ B = ⌈S·D·(u−l)/2 + H·ρ⌉`.

`w_k`, `ρ` and `B` are computed here in exact rational arithmetic from the decimal
midpoints, independently of the primary, and they **agree identically**:

| interval | ρ (mine = primary) | B (mine = primary) |
|---|---|---|
| jacobsen-2015-eigenvalue | 7686546224338895232565910573032297919/10^37 | 20077 |
| mertens-2022-p-med | 53129090819380871648563294193601/625·10^28 | 3000086 |
| mertens-2022-p-cell | 479982034581456010007919/625·10^21 | 80000077 |
| yang-zhou-2024-corrected | 1574738498838320187782581333490561/2441406250·10^27 | 100065 |

`|C(d,100)|` from Proposition 1 comes out `12175 / 3355121 / 749507743 / 157309446881`
for `d = 1..4`; the formula is checked against direct enumeration at `H = 3,5,7` for
`d = 1,2` before it is used.

### Four search paths, deliberately different code

| tag | structure | scale | where |
|---|---|---|---|
| `local` | nested enumeration, leading coefficient **solved exactly** | 10^15 | workstation |
| `loop` | nested enumeration, leading coefficient **looped** 100× more inner work | 10^15 | Huawei fleet, 60 shards/interval |
| `mitm` | split {0,1,2} \| {3,4}, sorted right half, binary search | 10^15 | Huawei fleet, 20 shards/interval |
| `alt` | same as `local` but on a **different screen lattice** | 2^50 | Huawei fleet, 60 shards/interval |

Shards partition the constant term `a_0 ∈ [−100,100]`, so the union of shards is the
class exactly. 560 work units over 10 containers / 144 vCPU, all completed.

`local`, `loop` and `mitm` return **byte-identical retained sets**; `alt` on the
2^50 lattice returns **exactly the same set** as well — the retained set does not
depend on the screen lattice, which is the property a correct screen must have.

### Exact decisions, own path

Retained tuples are decided in exact integer arithmetic on `P(x)·10^{dE}`, never in
binary floating point. Root counting uses **monotonicity from the second-derivative
bound plus endpoint signs**, not a Sturm chain: `|P'(m)| > (Σ k(k−1)|a_k|)(u−l)/2`
certifies that `P'` has no zero in the interval, hence `P` is strictly monotone and
the root count in `[l,u]` is the number of endpoint sign changes. No retained quartic
has a stationary point.

## Result: exact agreement on every cell

| interval | retention (prim.) mine / primary | near set mine / primary | root-bearing mine / primary | closest member | residual |
|---|---|---|---|---|---|
| jacobsen-2015-eigenvalue | 0 / 0 | 1543 / 1543 | 0 / 0 | agree | agree |
| mertens-2022-p-med | 3 / 3 | 1548 / 1548 | 1 / 1 | agree | agree |
| mertens-2022-p-cell | 127 / 127 | 1660 / 1660 | **15 / 15** | agree | agree |
| yang-zhou-2024-corrected | 0 / 0 | 1543 / 1543 | 0 / 0 | agree | agree |

The closest member is `[-84, 99, -7, 99, 58]` on three of the four intervals and
`[-97, 87, 54, 98, 49]` on `mertens-2022-p-cell`, matching the primary coefficient by
coefficient. Class sizes match (`157309446881`).

On `mertens-2022-p-cell` the 15 root-bearing quartics were additionally checked for
shared roots by exact rational polynomial gcd: **no pair shares a non-constant
factor**, so there are 15 distinct roots, and each contributes exactly one root in
the interval. `distinct_roots_in_interval = 15` therefore stands on its own, not on
the primary's assertion.

## The residual / mean-value check

The ticket asks for the check that rules out two implementations agreeing by both
reporting the same zero. The primary evaluates the closest member at the two
endpoints; the midpoint value is a genuinely different evaluation point, and the
mean value theorem caps the gap at `D_P (u−l)/2` with `D_P = Σ_k k|a_k|` for that
polynomial's own coefficients.

| interval | min over `[l,u]` | `|P(m)|` at the midpoint | gap | cap `D_P(u−l)/2` | holds |
|---|---|---|---|---|---|
| jacobsen-2015-eigenvalue | 1.22096352294980850756812524667737566506094010368e−9 | 1.225831e−9 | 4.867e−12 | 1.284e−11 | yes |
| mertens-2022-p-med | 0 (sign change) | — | — | 1.926e−9 | yes |
| mertens-2022-p-cell | 0 (sign change) | — | — | 5.480e−8 | yes |
| yang-zhou-2024-corrected | 5.930735364635135763e−10 | 6.174104e−10 | 2.434e−11 | 6.420e−11 | yes |

Two things follow. The gaps are **non-zero** and the residuals are **95×** and
**>1×** the mean-value cap, so the agreement is the agreement of two implementations
computing the same non-trivial quantity, not of two implementations both returning
zero. And on the two intervals where the residual is zero, it is zero for a reason
that the independent path also finds: the closest member changes sign across the
interval, i.e. it really does carry a root.

## What this replication does not do

- It does not make the manuscript's §7 "Recommended, not done" entry removable, and
  `Table 9`'s scope line is **not** edited here. Because the primary's committed
  numbers were read first, the honest wording is *"an independent re-derivation
  exists, but it was written with the committed summaries in hand"*, not *"the
  quartic census now has a blind second implementation"*. Suggested wording is given
  in the PR body for the owner to apply or amend.
- No new Monte Carlo, no new `p_c`, no merge, no `docs/STATUS.md` edit, no change to
  any existing result, script or manuscript file.
- Full repository CI has not been run.

## Reproduce

```sh
g++ -O3 -std=c++17 -o census_screen scripts/independent_census_screen.cpp
python scripts/independent_census.py 4 brute <outdir>          # weights, bound, search, decisions
python scripts/independent_census_agreement.py <repo> <outdir> agreement.json
```
