# P418 radius-resolved mechanism elimination

## Decision

The large P418 CRT-mask penalty is **not** a rejection within any single
archived radius.  It is a cross-radius transport/coherence obstruction: each
radius separately admits a mask-times-positive spectrum, while one spectrum
shared by radii four, five and six is rejected in every hand/charge channel.
The corresponding shared *raw* spectrum is not rejected.

| channel | raw sharing penalty | raw p | masked sharing penalty | masked p |
|---|---:|---:|---:|---:|
| plus r1 | 15.84 | 0.912 | 527.04 | 0.003984 |
| plus r2 | 30.15 | 0.183 | 617.17 | 0.003984 |
| minus r1 | 32.12 | 0.139 | 1093.32 | 0.003984 |
| minus r2 | 26.91 | 0.231 | 1111.51 | 0.003984 |

`0.003984 = 1/251` is the frozen 250-replicate bootstrap floor.  The sharing
penalty compares the common-spectrum cone distance with the sum of three
independent-radius cone distances and is calibrated under the fitted common
spectrum.

## What happens inside each radius

All twelve independent-radius masked fits survive at `alpha=0.01`.

- Radius four: raw and masked distances agree within `2.7e-11`; the mask is a
  reparameterization on the resolved coordinates for this shell.
- Radius five: the four raw-to-masked distance increments are `+1.95`, `+5.94`,
  `-1.97`, and `+0.52`; masked p values are `0.414`, `0.207`, `0.817`, and
  `0.733`.
- Radius six: both cones interpolate exactly (`d2=0`, `p=1`) because only 14
  covariance modes are resolved against a rank-14 design with 101 spectral
  coordinates.  This is saturation/nonidentification, not affirmative support.

The independent masked squared-distance sums are only `56.98`, `51.48`,
`42.68`, and `41.67`, compared with shared masked distances `584.02`, `668.65`,
`1135.99`, and `1153.18`.  Thus the earlier P418 rejection cannot be assigned
to single-radius stationarity or within-radius translation transport.

## Mechanism update

A fixed CRT gauge can de-gauge every observed shell independently, but the
recovered positive spectrum cannot be held fixed across the three radii under
the frozen observable-transport contract.  At least one of the common-spectrum
assumption, cross-radius observable normalization/landing semantics, or the
fixed de-gauging interpretation must flow with radius.  The archives do not
distinguish these alternatives.

## Statistical boundary

The calculation reuses the pinned radius-four 80k, radius-five 1.2M, and
radius-six 1.2M batch archives and their complete within-radius covariance.
Because the radii came from separate streams, the joint covariance is block
diagonal across radii.  No new Monte Carlo was run.  This result identifies no
continuum field, Jordan block, physical state count, or ordered-memory
mechanism.

