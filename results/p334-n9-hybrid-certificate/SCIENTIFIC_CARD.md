# P334 N9 hybrid certificate — scientific card

- **Question:** Does the corrected two-carrier reservoir saturate every deterministic N9 HNF/line row?
- **Design:** Close the shell with 24 immutable direct translation-orbit matchings and four independently recomputed twin-class capacitated Hall flows, all anchored to source commit `4bb75176c56558084c8397917995026e54420b9f`.
- **Result:** 22 of 28 rows saturate. The failures are exactly rows `1,3,6,9,15,24`, each with exact Hall deficiency `2160 = 5/16 × 6912`.
- **Heavy gate:** Missing direct rows `11,17,20,26` each saturate exactly at `45360/45360`; they are proved through the demand-9 coarse path, not inferred from neighboring rows.
- **Meaning:** The N9 obstruction is a six-row topology class, not an incomplete-shard artifact. The two-output-mark fixed-base MM repair remains the minimal known local repair for precisely this class; no arbitrary-HNF theorem is claimed.
