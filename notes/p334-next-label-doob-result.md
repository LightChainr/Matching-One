# P334: one next-label drives first-birth and completion responses together

Fresh nested continuations expose a positive next-label cross response between
the paired first-birth and completion observables. With
`B=Dnext` in the `(F1,F2)` basis,

```
Gamma = B12 = (B_AA-B_EE)/4,
A = F1+F2, E = -F1+F2  (paired contrasts).
```

| New conditional readout | N325 mean +/- original-batch SE | N425 mean +/- original-batch SE |
|---|---:|---:|
| canonical Gamma | 0.001381329 +/- 0.000044623 | 0.000997369 +/- 0.000040007 |
| integrated Gamma | 0.000022648 +/- 0.000001509 | 0.000011994 +/- 0.000000823 |

The shared label tends to move the two conditional responses in the same
direction: A adds this innovation and E partially cancels it. `K1<=K2` alone
does not force this sign; responses to different labels have no such total
order, especially after taking the paired orientation difference. This is a
finite-prefix shared-label observation, not identification of a new field or
a claim of path memory.

All of this first/completion cross response comes from the five joint-rank
cells containing R0. If both orientations have rank at least one, their K1 is
already determined by the prefix; its next-label innovation and cross response
are exactly zero by measurability. The score keeps every cell and the original
full-population denominator.

## Where the remaining response randomness lives

01+10 carries about 45.9% of the full K1 suffix variance and 33.1%/33.3% of
integrated E suffix variance at N325/N425. Within this group the next label
accounts for only:

| Coordinate | N325 next-label fraction | N425 next-label fraction |
|---|---:|---:|
| K1 | 14.540% +/- 0.507pp | 11.388% +/- 0.574pp |
| canonical E | 19.023% +/- 0.442pp | 15.974% +/- 0.390pp |
| integrated E | 5.128% +/- 0.489pp | 4.064% +/- 0.479pp |

The remaining suffix therefore dominates conditional variance, even though a
clear common next-label mode is present. These are innovation-variance shares,
not fractions of the signed E mean and not a Markov sufficiency test.

## Fresh-tail mean and dependence boundary

The N425 32-tail estimates are `E_integral=-0.000722835 +/- 0.00031093`,
`K1=-0.435483 +/- 0.21069`, `K2=-0.127555 +/- 0.16095` in the same H4
convention. Their differences from the old baseline and old safe estimator
have shared 20-batch covariance in the score. No independent pooling of these
prefix-reusing estimates is performed.

The new source `e32a85939279b8574278024d647b56d2d1485247` contains 640,000 new
tails per size, from 20,000 original prefixes and eight quartets each. Each
tail is shared by both orientations. Sixteen independent next-label draws,
each with two suffix replicas, leave covariance `B/16+Vafter/32`; the removed
conditional suffix covariance is `(15/16)B+(31/32)Vafter`. The outer product
of old baseline minus fresh finite-fork mean is not removed noise.

Artifacts: `results/p334-next-label-doob-quartets/score.json` and the two
`N*.joint_covariance_factor.json.gz` files. The factor supplies the complete
common covariance, including all nine cells' full matrices and old/new means,
without a dense covariance inversion. Dnext is retained signed, never clipped
to PSD. No new simulations, DP or validation suite were run by the scorer.
The next thin handoff adds the named cross-response magnitude and the already
computed gate readouts on these same batches, not another data block.
