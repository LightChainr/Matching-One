# P334 birth-age iota and site-pair control freeze

## Question

The production slope at `2e99533` survives line fixed effects and the stored
local landing/H4 marks.  The same sparse tables also retain exact saturation
indices and both birth sites.  This no-sample diagnostic asks whether the age
association is merely an omitted Smith/embedding class.

The risk set, age and outcome remain unchanged:

```text
K1 <= k0 < K2,
age = (k0-K1)/N,
y = 1[K2=k0+1].
```

First and second orientations are fitted separately but remain one common-
batch covariance block per size.

## Frozen exact fixed-effect ladder

Every layer uses the same one-degree within-stratum linear-probability slope:

1. replay: primitive projective line `ell`;
2. index pair: `(ell,iota01,iota12)`;
3. index plus birth-site pair:
   `(ell,iota01,iota12,gcd((site12-site01) mod N,N))`.

The last coordinate is not an arbitrary bin.  All four production geometries
have cyclic Smith form `[1,N]` and HNF `[[N,h12],[0,1]]`; the stored site label
is the cyclic coordinate.  The gcd of the relative displacement with `N` is
invariant under global translation, sign, and a change of cyclic generator.
It is the minimal exact Smith order class of the two birth sites.  Absolute
site labels are not used because they depend on the arbitrary torus origin.

No interaction, nonlinear age basis, flexible learner, or stratum selection is
allowed.  For interpretability the controlled within-stratum age denominator
must retain at least 25% of the line-only denominator in every orientation.

## Covariance and attenuation

Delete one batch index from both orientations simultaneously and recompute all
three slopes.  The resulting `6 x 6` covariance supplies individual Student-t
scores, two-orientation joint Wald scores, and the paired differences from the
primary slope.  Report both magnitude retention

```text
abs(beta_control)/abs(beta_primary)
```

and retained within-stratum age information.  These are nested views of the
same raw block and cannot be counted as independent evidence.

The line-only point and covariance must replay `2e99533` before the new
coordinates are interpreted.

## Temporal boundary

`iota01/site01` are known at the first birth.  `iota12/site12` are known only at
the future second birth, after the k0 predictive state.  Therefore the full
pair stratum is a retrospective mechanism diagnostic: it can reveal that the
age association tracks two-birth embedding structure, but it cannot define a
causal mediator or a predictive state observable at k0.

If the integral saturation theorem forces both stored indices to one, the
iota layer must be reported as an exact no-op rather than an estimated null.
Even a surviving site-controlled slope still may proxy microscopic geometry
at k0 that this archive does not store.
