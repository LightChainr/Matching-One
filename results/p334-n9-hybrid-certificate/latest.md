# P334 N9 hybrid corrected-reservoir certificate

Source commit: `4bb75176c56558084c8397917995026e54420b9f`.

This closes all deterministic N9 rows through two non-interchangeable exact paths: 24 byte-identified direct translation-orbit matchings and four recomputed coarse twin-class capacity flows.

| row | path | HNF | line | flow | deficiency | status |
|---:|---|---|---|---:|---:|---|
| 0 | direct_translation_orbit_matching | `[[3, 0], [0, 3]]` | `[0, 1]` | 6912/6912 | 0 | saturated |
| 1 | direct_translation_orbit_matching | `[[3, 0], [0, 3]]` | `[0, 1]` | 4752/6912 | 2160 | FAIL |
| 2 | direct_translation_orbit_matching | `[[3, 0], [0, 3]]` | `[1, 0]` | 6912/6912 | 0 | saturated |
| 3 | direct_translation_orbit_matching | `[[3, 0], [0, 3]]` | `[1, 0]` | 4752/6912 | 2160 | FAIL |
| 4 | direct_translation_orbit_matching | `[[3, 1], [0, 3]]` | `[0, 1]` | 4032/4032 | 0 | saturated |
| 5 | direct_translation_orbit_matching | `[[3, 1], [0, 3]]` | `[1, 0]` | 6912/6912 | 0 | saturated |
| 6 | direct_translation_orbit_matching | `[[3, 1], [0, 3]]` | `[1, 0]` | 4752/6912 | 2160 | FAIL |
| 7 | direct_translation_orbit_matching | `[[3, 2], [0, 3]]` | `[1, -1]` | 4032/4032 | 0 | saturated |
| 8 | direct_translation_orbit_matching | `[[3, 2], [0, 3]]` | `[1, 0]` | 6912/6912 | 0 | saturated |
| 9 | direct_translation_orbit_matching | `[[3, 2], [0, 3]]` | `[1, 0]` | 4752/6912 | 2160 | FAIL |
| 10 | direct_translation_orbit_matching | `[[9, 2], [0, 1]]` | `[0, 1]` | 11232/11232 | 0 | saturated |
| 11 | coarse_twin_capacitated_hall | `[[9, 2], [0, 1]]` | `[0, 1]` | 45360/45360 | 0 | saturated |
| 12 | direct_translation_orbit_matching | `[[9, 2], [0, 1]]` | `[0, 1]` | 7776/7776 | 0 | saturated |
| 13 | direct_translation_orbit_matching | `[[9, 3], [0, 1]]` | `[0, 1]` | 4032/4032 | 0 | saturated |
| 14 | direct_translation_orbit_matching | `[[9, 3], [0, 1]]` | `[1, -3]` | 6912/6912 | 0 | saturated |
| 15 | direct_translation_orbit_matching | `[[9, 3], [0, 1]]` | `[1, -3]` | 4752/6912 | 2160 | FAIL |
| 16 | direct_translation_orbit_matching | `[[9, 4], [0, 1]]` | `[1, -2]` | 11232/11232 | 0 | saturated |
| 17 | coarse_twin_capacitated_hall | `[[9, 4], [0, 1]]` | `[1, -2]` | 45360/45360 | 0 | saturated |
| 18 | direct_translation_orbit_matching | `[[9, 4], [0, 1]]` | `[1, -2]` | 7776/7776 | 0 | saturated |
| 19 | direct_translation_orbit_matching | `[[9, 5], [0, 1]]` | `[1, -2]` | 11232/11232 | 0 | saturated |
| 20 | coarse_twin_capacitated_hall | `[[9, 5], [0, 1]]` | `[1, -2]` | 45360/45360 | 0 | saturated |
| 21 | direct_translation_orbit_matching | `[[9, 5], [0, 1]]` | `[1, -2]` | 7776/7776 | 0 | saturated |
| 22 | direct_translation_orbit_matching | `[[9, 6], [0, 1]]` | `[1, -1]` | 4032/4032 | 0 | saturated |
| 23 | direct_translation_orbit_matching | `[[9, 6], [0, 1]]` | `[2, -3]` | 6912/6912 | 0 | saturated |
| 24 | direct_translation_orbit_matching | `[[9, 6], [0, 1]]` | `[2, -3]` | 4752/6912 | 2160 | FAIL |
| 25 | direct_translation_orbit_matching | `[[9, 7], [0, 1]]` | `[1, -1]` | 11232/11232 | 0 | saturated |
| 26 | coarse_twin_capacitated_hall | `[[9, 7], [0, 1]]` | `[1, -1]` | 45360/45360 | 0 | saturated |
| 27 | direct_translation_orbit_matching | `[[9, 7], [0, 1]]` | `[1, -1]` | 7776/7776 | 0 | saturated |

## Scientific card

- Question: Does the corrected two-carrier reservoir saturate every deterministic N9 HNF/line row?
- Design: 24 immutable direct orbit matchings plus four exact twin-class capacity flows
- Result: 22/28 saturate; failures are exactly rows 1,3,6,9,15,24, each deficient by 2160
- Heavy gate: rows 11,17,20,26 each saturate exactly at 45360/45360 after proved compression
- Meaning: The N9 obstruction is a six-row topology class, not a missing-heavy-row artifact; the repaired MM output-mark channel remains the minimal known local repair for those six rows.

The certificate is exhaustive only for the frozen deterministic N9 shell. It does not claim the arbitrary-HNF reservoir theorem.
