# A small early weight does not by itself approach the two-center boundary

| N | Fixed rank-two gap Delta2 | LOO SE | early weight | leverage Delta2/w0 | leverage SE |
|---|---:|---:|---:|---:|---:|
| 100 | 0.209677127 | 0.0104424 | 0.180589136 | 1.16107277 | 0.06773752 |
| 400 | 0.173267747 | 0.02874162 | 0.0653619969 | 2.65089432 | 0.5165757 |
| 900 | 0.172723152 | 0.04111921 | 0.032009977 | 5.39591615 | 1.622501 |

Changes, not independent of each other:
- N400_minus_N100: Delta2 change -0.0364093805 +/- 0.03057981.
- N900_minus_N400: Delta2 change -0.000544594328 +/- 0.05016842.
- N900_minus_N100: Delta2 change -0.0369539748 +/- 0.04242444.

Same three source blocks reused; no new MC, old Gaussian-boundary/center refits, or p-value combination. One gap plus its exact factorization into early weight and leverage. Affine invariance does not mean invariance under arbitrary thermal warps.

No Gaussian boundary p-value is assigned at the singular rank-two null.
