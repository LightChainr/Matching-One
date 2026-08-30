# P250 radius-six Level-S elimination certificate

## Certified statement

The locked P250 radius-six dependency block statistically eliminates the
gauge-free endpoint-Hankel classes `rank(H3)<=5`, `<=6`, and `<=7` separately
for plus and R2-gauged minus at the preregistered alpha `0.01`.

The certificate treats the model class invariantly:

```text
rank(H3)<=r
<=> dim ker(H3)>=10-r
<=> every (r+1)-minor vanishes.
```

The numerical statistic uses the pre-reveal maximum-volume Schur chart only as
a locally complete coordinate system for that determinantal variety.  The
excluded class itself does not depend on a fitted eigenbasis, root labels,
diagonalizability, or a plus/minus similarity gauge.

Rank eight is the first class not eliminated: `p=0.1978` for plus and
`p=0.1214` for minus.  The certificate therefore records a statistical lower
bound of eight, not an exact dimension or a physical field count.

## Machine check

Run:

```bash
python3 scripts/verify_p250_radius6_elimination_certificate.py \
  analysis/p250_radius6_level_s_elimination_certificate.json
```

The verifier fails closed unless it can:

1. match the authorized manifest, raw lock, old4/old5/fresh6 batches, exact
   gate, response, scorer, result manifest, and stored score hashes;
2. recompute the frozen score from all three locked batch streams;
3. reproduce the stored score exactly;
4. reproduce every certified T-squared, Hotelling F, covariance rank,
   denominator degree of freedom, alpha-.01 critical F, and p-value;
5. confirm that ranks 5--7 cross the frozen boundary, rank 8 does not, and the
   lower bound is eight in both hands;
6. confirm that the R2 kernel-projector bridge remains locked and untested.

No optimizer, random initialization, post-reveal pivot, or new Monte Carlo is
involved.

## Level and boundary

This is a dedicated **Level-S statistical elimination certificate** under
#370's vocabulary.  It is stronger than reporting a fit failure because the
rank constraints are invariant, the charts and covariance rule were frozen,
and a standalone verifier reconstructs the result from immutable artifacts.

It is not Level E and not an exact SOS/Positivstellensatz certificate.  The
confidence region remains the repository's 400-batch finite-Hotelling
calibration, not a rational interval enclosure of every Monte Carlo mean.  It
excludes only the declared rank classes, observer, geometry, dependency group,
and covariance convention.

The following remain explicitly untested:

- exact physical state dimension;
- rank-eight flatness at the next moment order;
- the full R2-conjugate kernel-projector bridge;
- ordered-path or noncommutative translation relations.

The bridge is not a negative result.  Its rank-five support prerequisite failed,
so the frozen scorer correctly never reached it.

## Scientific card

- Mechanism space changed: commuting endpoint realizations of dimension at
  most seven are statistically eliminated in each hand.
- Surviving class: rank at most eight is compatible only, not certified true.
- Evidence type: one Level-S conclusion from the shared old4/radius5/radius6
  dependency group.
- New production: none.
- Next lift: separately freeze a rank-eight cross-hand projector certificate if
  the R2 relation is the target; use higher moments only for rank-eight
  flatness.
