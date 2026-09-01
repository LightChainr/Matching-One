# P537: held-out L6 reveals an `L^-9/2` mixed-response fingerprint

## Production result

The frozen square-L6 (`N=36`) Bernoulli site-flip production completed
`100,000,000` off-root backgrounds in 100 counter-keyed batches.  It preserves
the exact L4/L5 contract: radius-one `ell_4`, the complete canonical pair
source, the `0->1` and `1->2` birth coordinates, matching-root conditioning,
and the root-coordinate mixed Hessian

```text
Psi4 = T_01*A_12 - T_12*A_01
S    = T_01 + T_12
C4   = 2*Psi4/S
G4   = [1/S] * partial_logit(C4).
```

The held-out result is

```text
p6 = 0.5923930092376508
G4 = 0.002297919737362416
L^4 G4 = 2.978103979621691 +/- 0.032636532053504
95% interval for G4 = [0.0022485620, 0.0023472775].
```

`S=-0.01906065786` remains safely separated from zero.  `Psi4` stays negative,
while `C4` and `G4` stay positive, so no sign reversal or denominator crossing
occurs between L4, L5, and L6.

## Frozen L^-4 decision

Before this production, the permissive L5-to-L6 `L^-4` continuation band was
frozen as

```text
G4(6) in [0.0022905, 0.0027995].
```

The L6 interval overlaps only the lower edge, so the declared decision is
`UNRESOLVED_OVERLAP`.  That wide-band decision is retained literally.

The sharper fixed-power comparison, however, changes the scientific picture.
Using the exact L5 value without refitting amplitudes gives

| fixed power | predicted `G4(6)` | observed-minus-prediction / SE |
|---:|---:|---:|
| `L^-4` | `0.0025450431` | `-9.81` |
| `L^-9/2` | `0.0023232959` | `-1.01` |
| `L^-5` | `0.0021208693` | `+7.03` |

Thus the exact `L^-4` point continuation and the neighboring integer
`L^-5` law are both strongly separated from the new data; the half-integer
`L^-9/2` law is not.

## Three-size phenomenon

The full L4/L5/L6 sequence is

| L | source | `G4` | `L^4 G4` | `L^(9/2) G4` |
|---:|---|---:|---:|---:|
| 4 | exact population | `0.01400841687` | `3.586154718` | `7.172309435` |
| 5 | exact population | `0.005277401459` | `3.298375912` | `7.375392754` |
| 6 | 100M production | `0.002297919737` | `2.978103980` | `7.294835151` |

The global L4-to-L6 effective power is `4.4582`; the L5-to-L6 power is
`4.5602`.  More importantly, `L^(9/2)G4` stays within about 2.8% across all
three sizes, while `L^4G4` decreases monotonically.

This exposes a new mechanism candidate:

```text
canonical pair normalization L^-4
times a root-conditioned landing mismatch L^-1/2
gives G4 ~ L^-9/2.
```

The half power has a natural percolation interpretation.  The alternating
four-arm channel has the expected exponent `5/4`, while the thermal scaling
direction has eigenvalue `y_t=3/4`; their difference is `1/2`.  A
root-conditioned landing response can therefore carry

```text
L^-4 * L^-(5/4-3/4) = L^-9/2.
```

This is a mechanism hypothesis, not a derivation or a rigorous square-site
exponent theorem.  It is more specific than the previous generic “local
N^-2 defect” interpretation and makes a new numerical prediction.

## Next frozen prediction

The mean of the three observed `L^(9/2)G4` amplitudes is `7.28084578`.
Before any L8 output is inspected, the half-integer mechanism predicts

```text
G4(8) = 0.000628459              [three-size amplitude]
G4(8) = 0.000629666              [direct L6 continuation].
```

For comparison, direct L6 continuation gives `0.000727076` for `L^-4` and
`0.000545307` for `L^-5`.  A single square-L8 production can distinguish the
three without changing the collar, source, minor, or geometry.

## Computation

The producer samples the off-root Bernoulli law at `p*=0.5923`, retains K for
exact local importance reweighting, and solves the matching root separately
in every delete-one replicate.  The global source mean uses one uniformly
counter-keyed endpoint per configuration as an unbiased estimator; the rare
landing events use the complete `O(N^2)` pair source exactly.  No descriptor,
radius, source, or minor was scanned.

The 16 local shards completed in about 66 seconds wall time.  The committed
100-batch sufficient-statistics table is enough to reproduce the score; the
100 million configurations themselves are not stored.

## Boundary

The L6 block is Monte Carlo, not a complete population, and its delete-one
normal interval is an uncertainty estimate rather than a rigorous enclosure.
The `9/2` observation was recognized from the three-size pattern and is
therefore exploratory at L6.  Its L8 prediction above is prospective.
