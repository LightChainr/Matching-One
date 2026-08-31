# N400 ordered-scalar target: still not identified, not a restored mechanism

The N100 source-defined scalar-warp ratios have now been read from an
independent N400 stream. They do **not** reject the scalar null, but their
uncertainty is too large to show that N100's tension has disappeared.
This outcome is separate from the density-warp analysis.

## Source, target and unchanged definition

N100 analysis commit `8ac86fb` selected the two ordered A height ratios:
central valley / first peak and second peak / first peak. It was committed
before this agent's N400 reveal, while N400 acquisition was already in
progress. This is a **source-frozen auxiliary readout before target reveal**,
not a claim of preregistration before target acquisition.

N400 was frozen at `894b3d8`, with every N100 period multiplied by two. Its
raw archive is committed at
[`3e01b49`](https://github.com/LightChainr/Matching-One/commit/3e01b495b5b637b0070705e37b4137a9a0ef0d8b),
8,000,000 permutations per shape pair, 400 common batches, seed
`20260831134001`, offset `267400000000`. N100 and N400 are independent;
within each scale, the three shapes form one common-stream block.

The target's empirical A derivative polynomials each have exactly three
interior roots, certified by exact integer-Bernstein subdivision. Their
ordinal correspondence is therefore the same two peaks and intervening
valley; no old p window is imposed, and no extra extrema are silently
dropped. All 400 LOO replicas retain the curvature types and ordered roots.

## Actual target readout

| A height ratio | N400 D | N400 U | U-D | same-batch SE |
|---|---:|---:|---:|---:|
| central valley / first peak | 0.0251833 | -0.0688078 | -0.0939911 | 0.4178808 |
| second peak / first peak | 0.8623188 | 0.8324037 | -0.0299152 | 0.1966682 |

The joint zero-null statistic is **chi2=0.0769793/2, p=0.96224**. A high
p-value does not establish the scalar mechanism. For comparison, N100's
differences were `(0.3061505,0.1607900)` with SEs `(0.1143783,0.0763016)`.
The independent N400-minus-N100 difference has `chi2=1.70508/2, p=0.42633`:
the cross-scale change is not resolved either.

The nominal sign flip of the N400 U valley is not a new topological
feature: the absolute central values are

```text
D_A:  0.00004948 +/- 0.00013578,
U_A:  0.00003705 +/- 0.00019999.
```

Both are unresolved. Empirical extra zero crossings around a tiny central
value would not certify population sign changes. The first peaks themselves
remain visible, but their amplitudes are only about 17–18% of N100's;
four times the samples do not preserve the precision of amplitude-normalized
fine-shape readouts. The full LOO covariance, amplitude uncertainty, and
finite-sample bias estimates are saved.

This does not mean all physical signal has vanished. It means these two
fine-shape contrasts, especially the target U valley, are weak at N400.
Other broad thermal-width readouts can have substantially better precision.

## Reproduce without new sampling

```sh
python3 scripts/p267_scalar_clock_transport.py \
  --source-directory /path/to/results/etop-n400-three-modulus \
  --source-commit 3e01b495b5b637b0070705e37b4137a9a0ef0d8b \
  --output results/p267-scalar-clock-n400
```

Every source blob is hash-checked against the declared commit. Numerical
root brackets follow exact ordered empirical root intervals, not N100's
locations. If the number of extrema changes, the script reports the new
pattern rather than forcing a three-landmark score.

## Scientific card

- **Mechanism space:** ordinary scalar transport remains unresolved at this
  new scale; N100's exploratory tension has not been independently reproduced,
  but its disappearance has not been established either.
- **Not proved:** no scalar recovery, no continuum fixed point, no new
  zero-crossing topology, and no relation to a density-map pass/fail.
- **Observer / source / geometry:** ordinary P4[A_top], identical ordered
  height ratios, homothetic N400 `2i,4i,1/2+i` period pairs, source 3e01b49.
- **Dependency:** one independent N400 block; full within-block covariance;
  zero additional Monte Carlo for this analysis.
- **Next useful science:** explain which broad profile component persists
  while the normalized valley contrast loses precision. A negative p-value
  ranking alone is not an explanation.
