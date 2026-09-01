# P537 full-T transport quotient

Status: `POST_HOC_SECONDARY__COMPLETE_CANONICAL_QUOTIENT`

The already-produced N65 global sufficient statistics give the complete
canonical-pair thermal response, after restoring the one C4-equivalent NN
source column omitted by the local-carrier producer:

```text
pooled root p        0.5927311266364432
T_t                 -7.43314231e-6  +/- 8.49959721e-7
J_65=A_65*T_t/M_t   -0.00162250989 +/- 0.00018553008
95% interval        [-0.00198614885, -0.00125887093]
```

The exact N25 full response is `J_25=-0.0055194314248394015`.  Therefore

```text
J_65/J_25 = 0.29396323 +/- 0.03361398
95% interval = [0.22807982, 0.35984664]
effective abs(J) power = 1.28130 +/- 0.11967
```

This is a substantial finite contraction of the complete normalized response,
not only of the contact-stage carrier.  The displayed effective power is a
two-point description.  It is not evidence for a fitted exponent and does not
establish `J_N -> 0` or the asymptotic little-o statement.

The frozen canonical selected-carrier total is `5.892%` of full `T_25` and
`2.551%` of full `T_65`.  Thus it does not saturate the finite full response.
The selected cells are not thermal-gauge invariant, so this share is a
coordinate-dependent diagnostic rather than a physical operator fraction.
No cell exponent or CFT interpretation is inferred from two sizes.

## The pooled-root displacement does not generate the finite signal

Reweighting the same N65 sufficient statistics to the prescribed square-site
reference `p_ref=0.592746050790` gives

```text
J_65(p_ref)                 -0.00162276136 +/- 0.00018552167
J_65(root)-J_65(p_ref)      +2.51467029e-7 +/- 9.06647139e-8
absolute transport fraction 1.54986e-4
```

Thus the observed displacement from this reference to the pooled root changes
the response by only about `0.0155%` of its magnitude within this finite
reweighting calculation.  The reference is not a rigorous enclosure of the
mathematical `p_c`; a uniform near-critical transport theorem remains the
exact-`p_c` gap.

## C4 reconstruction and dependencies

The 63 retained directions alone give `J=-0.00150763301`.  The omitted
`y=z=+e1` NN column, reconstructed by the C4-unbiased mean of
`-e1,+e2,-e2`, contributes
`-0.000114876876`, about `7.08%` of the complete value.  It cannot be silently
dropped.

The fill is exact at the expectation level under C4, not samplewise.  The
canonical Bell kernel agrees on all 4,140 partitions and all 16,560 tested
common 90-degree rotations; both N65 quotient tori have the required graph
automorphism.  The three retained NN batch estimates are statistically
compatible with equality (`Q=2.355` on 2 df), and their covariance propagates
through the mean fill.

Uncertainty combines two independent delete-one groups: the new 20M source
block and the P45 100M baseline, each with 100 batches.  Baseline omissions
resolve the pooled root before reweighting the source table.  The root-minus-
reference contrast is paired within each group before their variance
contributions are added.

No configurations were generated or replayed.  This score was chosen after
the N65 primary result and reuses the same dependency block, so it is a
secondary analysis rather than independent evidence.  The scorer, formulas,
input hashes and complete numerical output are in
[`RESULT.json`](RESULT.json).
