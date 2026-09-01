# P418: normalized production-archive reanalysis

This is a per-sample-unit correction of the existing radius4/radius5/radius6 production archives, not a new simulation.

All four joint raw/masked spectra are compatible at the original alpha=0.01. The historical radius-flow interpretation loses the support of this score. Compatibility does not identify a unique physical spectrum.

## Corrected shared-spectrum result

| channel | raw shared d2 (p) | masked shared d2 (p) | inherited masked sharing penalty (p; unreliable) |
|---|---:|---:|---:|
| plus_r1 | 70.878 (0.38247) | 70.878 (0.38247) | 17.9424 (0.7251) |
| plus_r2 | 75.692 (0.183267) | 75.692 (0.183267) | 33.13 (0.0677291) |
| minus_r1 | 76.7658 (0.175299) | 76.7658 (0.175299) | 32.8142 (0.0796813) |
| minus_r2 | 68.0672 (0.366534) | 68.0672 (0.366534) | 25.4723 (0.23506) |

At the unchanged alpha=0.01, corrected shared raw rejections: `[]`; corrected shared masked rejections: `[]`. Inherited sharing decisions are not scoreable because the radius5 fit is numerically unreliable.

## What changed from the historical result

| channel | historical masked shared d2 | corrected masked shared d2 | historical masked sharing penalty | corrected penalty |
|---|---:|---:|---:|---:|
| plus_r1 | 584.022 | 70.878 | 527.042 | 17.9424 |
| plus_r2 | 668.646 | 75.692 | 617.166 | 33.13 |
| minus_r1 | 1135.99 | 76.7658 | 1093.32 | 32.8142 |
| minus_r2 | 1153.18 | 68.0672 | 1111.51 | 25.4723 |

## Separate-radius compatibility

The radius5 rows and sharing penalties below are preserved outputs of the inherited solver, not reliable mechanism statistics: the saved-point diagnostic detects amplification of numerical null directions. Joint fits do not have this defect.

| channel | radius | raw d2 (p) | masked d2 (p) | resolved covariance modes |
|---|---|---:|---:|---:|
| plus_r1 | radius4 | 48.3248 (0.191235) | 48.3248 (0.191235) | 81 |
| plus_r1 | radius5 | 8.67998 (0.131474) | 4.61083 (0.350598) | 40 |
| plus_r1 | radius6 | 0 (1) | 0 (1) | 14 |
| plus_r2 | radius4 | 41.2627 (0.438247) | 41.2627 (0.438247) | 81 |
| plus_r2 | radius5 | 11.1425 (0.0278884) | 1.29924 (0.916335) | 40 |
| plus_r2 | radius6 | 0 (1) | 0 (1) | 14 |
| minus_r1 | radius4 | 41.2394 (0.478088) | 41.2394 (0.478088) | 81 |
| minus_r1 | radius5 | 3.88797 (0.661355) | 2.71218 (0.601594) | 40 |
| minus_r1 | radius6 | 0 (1) | 0 (1) | 14 |
| minus_r2 | radius4 | 39.2154 (0.494024) | 39.2154 (0.494024) | 81 |
| minus_r2 | radius5 | 8.20273 (0.131474) | 3.37941 (0.609562) | 40 |
| minus_r2 | radius6 | 0 (1) | 0 (1) | 14 |

## Method and finite scope

The CSV coordinates are sums over each batch. Every radius4 batch contains 200 samples, whereas every radius5/radius6 batch contains 3000. The original reader used these sums as observations without rescaling the common Fourier design. This wrapper divides each row by its own `samples` before calling the unchanged historical whitening and family fitter.

For block exposure n, replacing batch sums by per-sample means sends mean to mean/n and covariance to covariance/n². The whitened response is unchanged (up to the eigensystem basis), but the correctly whitened design acquires the factor n. An independent spectrum for each radius can absorb this factor; one common spectrum cannot absorb three inconsistent units. Thus the old common-spectrum rejection cannot itself be interpreted as physical radius flow.

The method keeps the original covariance-of-the-batch-mean formula, relative eigenvalue cutoff 1e-10, nonnegative 101-frequency Fourier cone, exact CRT mask, 250 Gaussian parametric-bootstrap draws, original channel seeds and alpha=0.01. The sharing penalty is calibrated under the common fitted center for both nested fits. `score.json` retains every original family-fit output plus the point spectra; `normalized-inputs.json` retains the normalized means and full within-channel/radius covariance matrices.

Different radius streams retain the original block-diagonal covariance convention. Hands/charges in one radius share configurations; their separate p values are not independent pieces of evidence and are not combined. The bootstrap draws are uncertainty calculations on archived means, not new Monte Carlo lattice configurations. Basis/sign choices in covariance eigendecomposition and NNLS nonuniqueness can change finite-bootstrap details without changing the statistical convention.

Non-rejection is compatibility of this finite spatial-spectrum model, not a unique recovered spectrum, physical state count, local field identification, or proof of a continuum mechanism. A saturated radius6 fit is not affirmative model identification. If corrected tension remains, it belongs to this normalized finite observation contract; no new mechanism is assigned automatically.

The radius4 block observes zero displacement, while radius5 and radius6 do not. Adding a constant to every spectral weight changes only zero displacement, so positivity alone cannot restrict an isolated nonzero-displacement shell beyond its signed Fourier span. This argument does not trivialize the common fit, whose total spectral mass is constrained by radius4's zero-lag row. Equal common raw/masked distances here are not a proof that their cones are always equal.

A small saved-point diagnostic (solver-note.json; no NNLS or bootstrap rerun) finds all eight joint fits consistent with their reconstructed residuals within 7e-13, with the resolved-rank-69 least-squares floor within 7e-11 and scaled KKT violations below 1e-14. The radius5 fits instead amplify floating-point null directions: their saved-weight residuals fall about 10–24 below the exact/resolved-rank-20 Fourier-span floor and scaled negative-gradient violations reach 0.066. The missing zero-lag row permits arbitrarily large uniform null mass in that isolated shell. Those radius5 distances and the derived sharing penalties are numerically unreliable, not physical changes caused by rescaling. The robust result is the disappearance of the large joint masked penalty; the inherited sharing-penalty p values are not used for mechanism inference.

The exact CRT/root-translation certificates, correctly normalized P250 Hankel/radius5/radius6 scorers and paired-anchor pilot's own statistics are not modified. Historical sum-unit outputs remain available; they are included for comparison, not silently overwritten.

## Provenance and reproduction

- Data, geometry and mask source: `8704eee790403e14e5ad75d3465ee1496eaa9c0e`.
- Radius-family fitter and historical comparison: `588ca452dedd47213a424d79fc119ad67f8f77df`.
- Input byte SHA256 values must match the historical score's frozen inputs. No replacement archive is accepted.
- Exact imported module commits and SHA256 values, Python/NumPy/SciPy/BLAS information, elapsed time and command are saved in `score.json` and `manifest.json`.

```bash
/Users/lc/python-envs/research-py311/bin/python scripts/analyze_p418_normalized_archive.py
/Users/lc/python-envs/research-py311/bin/python scripts/analyze_p418_normalized_archive.py --diagnose-existing
```

Elapsed wall time: 7.287 seconds.
