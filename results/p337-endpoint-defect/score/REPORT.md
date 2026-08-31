# One saturation defect: the closed source and original global U

Primary decision: **source_independent_endpoint_gain_rejected**.
Secondary fixed question: **interior_thermal_only_mixed_null_rejected**.

| Endpoint quantity | Numerical evaluation of exact coefficients |
|---|---:|
| p0_endpoint | 0.407334460671773 |
| U | 2.71572887734847 |
| U_s | 3.70824092928232 |
| U_t | 0.389147178497717 |
| U_st | 10.7557184075641 |
| Xi_U_t_epsilon | -10.7557184075641 |
| R | 27.7665635812302 |
| d_s_Ut_over_U | 3.76486418289627 |
| source_free_gain_slope | 1.36546801862744 |
| gain_predicted_U_st | 0.531368026777735 |
| U_st_minus_gain_prediction | 10.2243503807863 |

The source, graph pair and original global-U normalization are fixed.
U_s and U_st use s increasing toward full saturation; epsilon=1-s, so
Xi=U_t,epsilon=-U_st. The primary R=U*U_st-U_s*U_t eliminates an unknown
source-independent geometric gain. Decisions use the rational enclosures
of R/A50^2 and Xi/A50, not the rounded displayed values. Every ratio term,
root displacement and slope displacement is included in `score.json`.

One removed A site represents all25 positions by translation. Each geometry
enumerates all2^25 free-B configurations exactly. Its Bernstein degree is25.
The intact endpoint uses the prior N25 coefficients by complement, with no
new baseline enumeration or resampling. Source normalizers are retained
separately for each geometry before the pooled root and angular projection.
The p-dependent defect dose25(1-p) is differentiated, not frozen at the root.

This is a finite mechanistic equality test, not a confidence interval or
a continuum field identification. A nonzero R excludes the local scalar-gain
extension, not the exact saturated identity. No extra source, fitted curve,
defect class, sample extension or old-experiment rescue enters this result.
The original F4/P154/P334 stop decisions remain unchanged.

Frozen specification: `notes/checkerboard-endpoint-defect-decision-freeze.md`
at9024fdbf. See the producer contract and receipts in the parent result folder.
