# Identical immediate rank and Euler increments, different future birth responses

The new contact readout finds a clear distinction **inside topology-safe
insertions with identical contact degree**. Replacing component-merging
attachments by contractible-cycle attachments is associated with later
first and second essential births, and with a longer rank-one lifetime.
The direction is resolved at both N325 and N425. This supplies a concrete
connectivity coordinate beyond the current rank/Euler summary.

## Fixed-rank, fixed-Euler discriminator

For an R0 prefix let a next site touch e occupied edges belonging to c old
components. A safe insertion introduces `loop=e-c` contractible graph cycles
and `merge=max(c-1,0)` component mergers. Its Euler increment is `1-e`.
The [named pre-readout discriminator](https://github.com/LightChainr/Matching-One/blob/2b52d337/notes/p334-safe-contact-euler-discriminant.md)
compares independently sampled U,V only when both preserve the selected
orientation's R0 rank and `e(U)=e(V)`.

Define, on this mask, the equal-orientation, full-prefix-weighted statistics

```
G = E[mask * (loop(U)-loop(V))^2 / 2],
H_Y = E[mask * (loop(U)-loop(V)) * (Ybar_U-Ybar_V) / 2],
s_Y = H_Y/G.
```

Each Ybar uses the two existing suffixes under that label. The numerator is
an unbiased estimate of a weighted within-prefix covariance of the true
conditional mean; s is its pooled linear-projection slope. No per-prefix
class frequency or fitted nuisance model is required.

For each prefix, the exact target is a sum over degrees e with weights
`pi_safe,e^2`. The hypothesis `m(Z,u)=f_Z(e(u))` within safe labels gives
H_Y=0 for every prefix, however complicated or prefix-specific f_Z is.
The nonzero observed H thus demonstrates failure of this particular
rank/Euler-only next-response closure in the sampled population.

## Main measured direction

The slopes below are birth-clock steps per unit of loop-count contrast in
that weighted linear projection, not the causal effect of editing a loop
while holding the rest of a configuration fixed.

| Future clock | N325 slope +/- original-batch SE | N425 slope +/- original-batch SE |
|---|---:|---:|
| first birth K1 | +.590092 +/- .046187 | +.710190 +/- .062401 |
| completion K2 | +.885115 +/- .059950 | +.965016 +/- .094448 |
| center C=(K1+K2)/2 | +.737603 +/- .045094 | +.837603 +/- .069171 |
| lifetime W=K2-K1 | +.295024 +/- .057625 | +.254826 +/- .080560 |

The unnormalized same-degree covariance numerators are themselves clearly
nonzero, so the conclusion does not come from a ratio alone:

| Numerator | N325 H +/- SE | N425 H +/- SE |
|---|---:|---:|
| loop versus K1 | .01035703 +/- .00080415 | .01362344 +/- .00128206 |
| loop versus K2 | .01553516 +/- .00105874 | .01851172 +/- .00191199 |

At the fixed `p_ref=.59274605079`, the corresponding F1 slopes are
`-.0140533 +/- .0011013` and `-.0140084 +/- .0012505`; F2 slopes are
`-.0093252 +/- .0007785` and `-.0086506 +/- .0009870`. Both near-critical
birth responses decrease along this attachment direction. The close F1
magnitudes across these two sizes are a useful exploratory fingerprint,
not an estimated asymptotic exponent or a proof of universality.

At fixed e, `delta loop=-delta merge` exactly. The merger slopes are the
opposites of the table and are **the same finding**, not independent
corroboration. Isolation cannot vary within an equal-degree pair; its raw
zero covariance is retained but its undefined 0/0 slopes are omitted.

Without the equal-degree restriction, the loop/K1 slopes are only
`+.11360 +/- .01860` and `+.17767 +/- .01919`, while merger/K1 slopes are
`-.56282 +/- .02726` and `-.66309 +/- .03579`. Contact degree and partition
are different coordinates; combining them in a single untyped attachment
count obscures this distinction.

## Physical interpretation and an independent mathematical witness

Closing a zero-winding cycle spends an insertion without necessarily
connecting previously separate lift footprints. A merger can change which
contact addresses become comparable and hence alter later essential-cycle
formation. The population result supports this direction at a fixed Euler
increment, rather than equating every new graph cycle with progress toward
an ambient torus birth.

The new [exact torus witness](https://github.com/LightChainr/Matching-One/blob/21bdb7b0e59155639452e26f3e75833234bfdaa5/notes/p334-safe-role-innovations-and-contact-mechanism.md)
is complementary: from one R0 tree, a safe extension gives one-further-step
first/completion probabilities `(4/16,1/16)`, while a safe local-loop
attachment gives `(2/16,0)`. That paper example fixes occupied count and
rank but not contact degree. It explains why positive first/completion
cooperation does not by itself mean loop facilitation; the measured
same-degree population comparison provides the additional Euler separation.

## Source and relation to the Gamma result

This reader joins the completed conditional suffix source
`e32a85939279b8574278024d647b56d2d1485247` to checkpoint-only contact marks
`959a7fa26677c416b874d272f1ba66523fb38f73`. It acquires no new suffixes and
solves no conditional network. All 40,000 original prefixes remain in the
denominators; the two orientations are pooled with weight one half. The
same twenty original batches per size provide all errors, C/W contrasts
and saved covariance factors.

The domain here is **each orientation's own R0-safe labels**. This is not
the other analysis's both-orientations-safe mask, nor the paired H4
contrast. It identifies a finite microscopic response mechanism, without
claiming that this coordinate already explains global H4 or the complete
E mean. The new forks, mask decomposition and this contact result share
one original-prefix dependency group and must not be counted as independent
population replications.

Artifacts: `scripts/p334_safe_contact_response.py` and
`results/p334-safe-contact-response/score.json`. The latter retains all
87 raw per-batch coordinates, derived C/W contrasts and joint LOO factors.
The canonical and integrated results are correlated projections of those
same trajectories, with no covariance inverse or PSD clipping.
