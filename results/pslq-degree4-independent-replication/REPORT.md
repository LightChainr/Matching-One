# REPORT — independent degree ≤ 4, height ≤ 100 census (issue #568)

See `notes/p568-independent-census-20260913.md` for the full write-up.

Four search paths (local solve-based brute force, fleet naive `loop`, fleet `mitm`
split, fleet `alt` on the 2^50 lattice) return identical retained sets, and those
sets reproduce the primary's committed per-interval retention, near-set sizes,
root-bearing counts, closest members and exact residuals.

| interval | retention | near | root-bearing | excluded |
|---|---|---|---|---|
| jacobsen-2015-eigenvalue | 0 | 1543 | 0 | true |
| mertens-2022-p-med | 3 | 1548 | 1 | false |
| mertens-2022-p-cell | 127 | 1660 | 15 | false |
| yang-zhou-2024-corrected | 0 | 1543 | 0 | true |

Machine-readable: `latest.json`, `derived/<interval>.json`, `agreement.json`.
