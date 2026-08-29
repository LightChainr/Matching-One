# Gaussian x annulus semantic crosswalk

## Verdict

No revealed PR277 row and P253 row has the same source, observer, normalization, and transfer coordinate. There is therefore no numerical 2x2 rectangle to score from the old archives.

| archive row | source/readout | sector | transfer coordinate | cross-context eligible? |
|---|---|---|---|---|
| PR277 scalar U | global threshold-rank thermal derivative | matching-odd scalar side view | cover generation | no |
| PR277 r2-r6 | global Hermite-Krawtchouk derivative jet | matching-odd | cover generation | no |
| P253 A_plus | fixed-p root-toggle landing H4 | ordinary/matching-even | log2(R/2) | target row |
| P253 A_minus | fixed-p root-toggle landing H4 | matching-odd | log2(R/2) | target row |

Scaling or renaming cannot turn a product-measure derivative into a conditional root-toggle observable.

## Frozen row choice

The existing annulus context selected `N425` for the correlated `(A_plus,A_minus)` pair. Its minimum adjacent-lambda Mahalanobis separation is `0.0433206`. This is weak; it freezes the row/geometry, not a mechanism conclusion.

## Missing cells

| context | A_plus ordinary | A_minus matching-odd |
|---|---|---|
| annulus radius doubling | existing P253 N425 | existing P253 N425 |
| Gaussian cover doubling | missing | missing |

The missing cells must use the P253 fixed-p root-toggle and landing-shell definitions verbatim on a norm-4 cover chain. The existing cyclic multiradius runner accepts the primitive N85/N170 parents but rejects the nonprimitive N340/N680 children, so the frozen acquisition requires only a general-period geometry adapter, not a new observable framework.

## Frozen score

The shared model uses one of `lambda={0,1/2,1}` in both contexts and both rows. The minimal context-enriched adversary allows one lambda per context while keeping it shared across A_plus/A_minus. All fixed scores and the predeclared bootstrap comparison use complete within-context covariance; the future Gaussian block is independent of P253.

## Scientific card

- MECHANISM QUESTION: can one low-dimensional generator transport the same local pivotal rows through cover and radius contexts?
- SEMANTIC VERDICT: no PR277 global-rank row is identical to a P253 local landing-H4 row, so no current numerical rectangle exists.
- SELECTED ROWS: N425 A_plus ordinary and A_minus matching-odd, chosen jointly with the full annulus covariance.
- MISSING CELLS: the same fixed-p root-toggle A_plus/A_minus readout along one norm-4 Gaussian cover chain.
- DECISION: compare shared lambda with context-specific lambda only after the missing cells pass the general-period runner preflight.
