# P205 quotient-character `A_top/E_top` ray test

## Answer

The completed P205 N25/N50/N125 equal-area quotient contrasts do **not**
show a detectable Smith/quotient-character rotation of the
`P4(A_top,E_top)` ray.  At each quotient pair's intrinsic matching center,
the common-ray fit is

```text
angle = -14.9089 degrees
E_top / A_top = -0.266246
chi2 = 2.21999 on 2 df
p = 0.329561
```

The independent common-`p_ref=0.59274605079` sensitivity is essentially
unchanged: `angle=-14.8673 degrees`, `chi2=2.27573/2`, `p=0.320503`.
All three pairwise ray-sharing penalties also survive (`p=0.2880, 0.2914,
0.3688`).

This gives a direct answer to the motivating question: changing the Smith
pair across the completed P205 prism does not explain the existing
production global-ray failure.  Instead, the prism supplies one stable
quotient-contrast direction that is itself incompatible with the old
eight-row global direction (`-30.6185 degrees`; separate-versus-common
`Delta chi2=38.7150/1`, `p=4.90e-10`).  At the prespecified 0.01 threshold,
the prism direction is incompatible with P43 (`Delta chi2=43.7994/1`,
`p=3.64e-11`) but not eliminated against P49, P50, or P57 individually.

## Reconstructed rows

| quotient pair | intrinsic p0 | P4(A_top) | P4(E_top) | SE(A), SE(E) |
|---|---:|---:|---:|---:|
| N25 | 0.5926739969 | 0.004318746 | -0.001130424 | 7.24e-5, 5.73e-5 |
| N50 | 0.5927671259 | 0.001304370 | -0.000420000 | 9.63e-5, 6.49e-5 |
| N125 | 0.5927445153 | 0.000230058 | -0.000195226 | 2.13e-4, 1.12e-4 |

The complete 2x2 covariance matrix for every row, source-blob SHA-256
hashes, ray profiles, optimizer certificates, and pairwise scores are in
`latest.json`.

## Scientific card

- **Mechanism space changed:** a quotient/Smith-dependent rotation versus a
  common quotient-contrast ray.
- **Result:** the common ray survives across all three completed quotient
  pairs.  Quotient character is therefore not the observed source of the
  earlier within-production ray heterogeneity.
- **Not proved:** this does not make the quotient ray universal, assign it to
  a field, or turn the quotient contrasts into individual geometry states.
- **Observer / sector / source / geometry:** `P4(A_top,E_top)` / Alexander
  rank plane / threshold-rank histograms / P205 equal-area quotient pairs.
- **Dependency:** N25, N50, and N125 use independent completed raw blocks;
  each row retains the full same-batch A/E covariance.  The production
  comparison is block diagonal and is not counted as another discovery.
- **Natural upweight:** acquire or reuse another equal-area quotient contrast
  selected to rotate the P205 ray relative to P43; do not rerun harmonic
  voting.

## Coordinate and claim boundary

The primary comparison uses an intrinsic center reconstructed separately for
each quotient pair and the declared first-order center influence covariance.
The common-fixed-`p_ref` block is reported only as an internal P205
sensitivity.  It is deliberately **not** combined with the production
fixed-center rows, because those condition on each archive's own plug-in
center and are not the same coordinate.

These rows are first-minus-second quotient contrasts normalized by exact
`Delta cos(4 theta)`, not observations of either quotient in isolation.

## Reproduction

```bash
python3 scripts/p205_quotient_etop_ray_test.py \
  --output results/etop-p205-quotient-ray-test/latest.json
python3 -m unittest discover -s tests \
  -p 'test_p205_quotient_etop_ray_test.py'
```

No new Monte Carlo samples were generated.  The scorer reads the locked P205
raw blobs from commit `fc14817bb8c0b2f6e7cbde41778e715dcb62bc64`.
