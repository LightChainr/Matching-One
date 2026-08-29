# Gaussian information-per-CPU optimizer v1

This first executable design layer ranks primitive Gaussian multiplier lineages
without reading any candidate target result.  Exact integer genealogy, Smith
invariants and H4/H8/H12 fractions are kept in the machine-readable table.

Expected simple-Gaussian separation uses

```text
KL(model i || model j) = (mu_i-mu_j)^2 / (2 variance).
```

The target variance and CPU cost are two-point power-law planning fits to the
committed N325/N425 pilot.  They are deliberately labeled machine-specific and
are inflated in the robust maximin view.  They must be replaced or extended
when general-period backend timings become the relevant design space.

Two rankings are reported:

1. H4 versus H12 information per CPU, which regression-checks the successful
   norm-5 geometry;
2. worst-pair information across H4 x=21/4, H4 x=17/4, H8, H12 and zero.

The second can be dominated by a nearly degenerate adversary and is therefore a
diagnostic, not an automatic production order.  A selected candidate still
requires a frozen typed observable, covariance plan and allocation before its
first target batch.

Reproduce with:

```bash
python scripts/design_next_gaussian_experiment.py \
  analysis/gaussian_experiment_design_manifest.yaml \
  --json-output results/gaussian-experiment-design/latest.json \
  --markdown-output results/gaussian-experiment-design/latest.md
```
