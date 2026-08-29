# P275 A_top field-identity Phase 1

## Frozen answer

None of the five frozen subspaces survives. The full 18-real covariance score
selects `none`:

| frozen model | chi-square / df | survival p |
|---|---:|---:|
| ordinary thermal `Q4 epsilon` | 222271475.8 / 16 | numerical zero |
| thermal `Q4` energy-Jordan | 19334819.9 / 12 | numerical zero |
| generic pure H4 completion | 203874360.9 / 12 | numerical zero |
| generic affine-log H4 completion | 3984901.7 / 6 | numerical zero |
| zero response | 520818635.3 / 18 | numerical zero |

This does not reopen H8/H12 and it does not reject the independently selected
global H4 or thermal Q4. It rejects the selector's scaling map from this
particular rank-birth line source to an `N^-13/8` correction amplitude.

## Nine finite-root responses

Every `p_N` was solved from `M_N(p)=E_p[A_top]=0` inside each delete-one
replicate before forming the connected response. `Y=N^(13/8) Gamma` is the
quantity seen by the frozen Q4/Jordan models.

| geometry | p_N | Re Y | Im Y | birth mass B |
|---|---:|---:|---:|---:|
| N50/i | 0.592615702394 | 197.52875 | 0.00133 | 7.60178 |
| N50/2i | 0.593162621815 | 250.31170 | -0.10077 | 6.09156 |
| N50/5i/2 | 0.593010690818 | 268.17776 | -0.14322 | 5.01589 |
| N130/i | 0.592737699040 | 914.58739 | -0.02777 | 10.84197 |
| N130/2i | 0.592688170172 | 1174.93927 | -0.24546 | 8.70161 |
| N130/5i/2 | 0.592830160517 | 1258.76020 | -0.25786 | 7.19970 |
| N170/i | 0.592721625884 | 1410.77150 | 0.02735 | 11.98037 |
| N170/2i | 0.592726207363 | 1814.54380 | -0.51547 | 9.62574 |
| N170/5i/2 | 0.592709301529 | 1946.97470 | -0.30331 | 7.95786 |

## Discovery reading: scale-zero response

The canonical unscaled response does not decay as `N^-13/8`. It approaches an
order-one, almost real modulus profile:

| modulus | Re Gamma at N50 | N130 | N170 | last effective decay exponent |
|---|---:|---:|---:|---:|
| i | 0.3426132 | 0.3357896 | 0.3349478 | 0.00936 |
| 2i | 0.4341651 | 0.4313775 | 0.4308122 | 0.00489 |
| 5i/2 | 0.4651537 | 0.4621522 | 0.4622542 | -0.00082 |

The N170 real standard errors are respectively `9.66e-5`, `6.27e-5`, and
`4.62e-5`. After exact physical-line phase transport the imaginary parts are
only `6.5e-6`, `-1.22e-4`, and `-7.2e-5`. Thus the failure is not a lost sign or
a winding-coordinate frame error.

At N170 the normalized modulus ratios are

```text
Gamma(2i)/Gamma(i)       = 1.28621
Gamma(5i/2)/Gamma(i)     = 1.38008
```

instead of the ordinary-Q4 targets `2.75` and `4.29344`. Multiplying a
scale-zero response by `N^(13/8)` creates the rapidly growing Y rows that
reject even the generic affine-log model.

The raw canonical covariance `C=Gamma*B` shows what the birth-mass
normalization changes:

| N | Re C(i) | Re C(2i) | Re C(5i/2) |
|---:|---:|---:|---:|
| 50 | 2.60447 | 2.64474 | 2.33316 |
| 130 | 3.64062 | 3.75368 | 3.32736 |
| 170 | 4.01280 | 4.14689 | 3.67855 |

At N170 its two non-square ratios to `i` are `1.03342` and `0.91671`.
Dividing by B rotates those to `1.28621` and `1.38008`, so B normalization
materially changes the modulus shape, but neither raw nor normalized shape is
the E4 shape. The failure cannot be repaired by undoing the normalization.

The economical mechanism is that `chi4(ell)` is not a local Q4 insertion.
`ell` is the global primitive ambient-H1 line at a rank birth, and division by
the birth mass turns the observable into a conditional projective-line
polarization. Such a scale-zero integrated/zero-mode Ward background is
allowed to remain order one. The raw covariance consequently tracks the
leading birth mass instead of an `x=21/4` irrelevant-field amplitude.

This is a productive negative result: H4 remains selected in the global
matching channel, but the line-resolved rank-birth source is not its local
field-identity probe. Commit `6899b119db5b16e9918db53abf5280d990eb6653`
already preserves, on the same marked-birth stream, the first/second local
landing-H4 `S/D` marks and their `q` products. The next preregistration should
use that stream to compare two fixed landing radii on the same birth path and
form a frozen UV-annihilator before the modulus score. This reuses common-field
covariance while removing the local/zero-mode component; it must not rerun or
rescale the global `chi4(ell)` observable and must not reopen H8/H12.

## Execution and integrity

- Frozen selector authorization: `3f1f4f81fb5bc6b563dab51436f4cca146e855dd`.
- Runner commit: `cb83673fb5f221616a47d53f564635c11e7d0680`.
- Identical ARM64 binary SHA256 on all hosts:
  `462cd70be52b18bba69d1337b551110bb0df1527a1bacaebda85215e13ffb356`.
- 20,000,000 permutations per geometry, 100 batches, counters
  `[9500000000,9520000000)`; 180,000,000 geometry-replicas total.
- Within every N, all three priority-field digest sequences are byte-identical.
- Remote and local SHA256 values match for every batch CSV and metadata file.
- N50 ran on ZyTrST, N130 on XPk2PZ, and N170 on HZsCM6.

The first N50 launch failed before entering the runner because that image lacks
`/usr/bin/time`; the 72-byte failure logs are retained. The retry used the same
frozen counter interval and runner, and runner metadata records its own elapsed
time. No target data existed from the failed launch.
