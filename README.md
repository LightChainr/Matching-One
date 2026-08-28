# Matching One

Matching One is a computational research project on square-lattice site percolation, its matching-lattice identity, and the finite-size structure behind the threshold.

The exact structural anchor is

\[
p_c^{\mathrm{site}}(\mathbb Z^2)
+ p_c^{\mathrm{site}}(\mathrm{NN+NNN}) = 1.
\]

## What we currently observe

The strongest current numerical evidence is an orientation-dependent finite-size sector on primitive Gaussian tori.

1. **Independent five-size confirmation.** With 100 million paired replicas per size, the frozen same-`N` orientation differences at `N=65,85,130,145,170` all have the sign predicted by `Delta cos(4 theta)`. The corresponding z-scores are `16.03, 11.23, 5.22, 5.27, 2.58`; the pooled scaled amplitude is

   ```text
   A4 = N^(13/8) DeltaM / DeltaCos4
      = 0.7885 +/- 0.0352.
   ```

2. **Two parameter-free Gaussian-doubling tests.** Multiplication by `1+i` doubles `N` and rotates the microscopic square lattice by `pi/4`. The frozen prediction

   ```text
   DeltaM(2N) / DeltaM(N) = -2^(-13/8) = -0.3242098887...
   ```

   gives fresh observed ratios

   ```text
   65 -> 130: -0.31382 +/- 0.0908
   85 -> 170: -0.34095 +/- 0.1118.
   ```

3. **A third prospective lineage also passes.** For the independently frozen `145 -> 290` lineage,

   ```text
   DeltaM_290 = -0.000160648 +/- 0.000040542
   frozen target = -0.0001376564 +/- 0.000024997.
   ```

4. **The residual moves the finite matching root in the expected way.** Threshold-rank reconstructions give

   ```text
   -DeltaRoot * mean(M') / DeltaM ~= 1
   ```

   across the tested sizes. A clean high-stat root-amplitude test gives `A_p=0.4203 +/- 0.0216` at `N=65` and `0.3949 +/- 0.0308` at `N=85`, against the frozen prediction `0.4510 +/- 0.0201`.

The numerical archive, production scripts, exact checks, failed models, and raw sufficient statistics are now on `main` under `results/server-20260828/`.

## Working interpretation

A compact working law is

\[
\Delta M_N
\approx A\,\Delta\cos(4\theta)\,N^{-13/8}
\]

with Gaussian-integer multiplication acting simultaneously on the radial scale and microscopic orientation.

The evidence is strongest for an **odd square-harmonic orientation sector with approximately `N^-13/8` radial behavior**. It is not yet enough to say that the leading harmonic is uniquely H4 rather than H12/H20/etc., or that the corresponding continuum field has been uniquely identified as the proposed `x=21/4` LCFT operator.

For the current integrated research view, see [`notes/SYNTHESIS-20260828.md`](notes/SYNTHESIS-20260828.md). The more formal claim ledger remains [`docs/STATUS.md`](docs/STATUS.md).

## Next three experiments

The project is currently optimized around three discriminators:

1. **Full-curve Gaussian doubling triptych** (#49/#50): test `DeltaM`, slope, and root ratios on the three exact lineages, including the parameter-free root target `-1/4`.
2. **Norm-5 H4 versus H12 test** (#57): use a Gaussian multiplier for which the two harmonic hypotheses predict different signs/magnitudes.
3. **Exact parity control** (#44, then #42/#48): separate generic square-lattice anisotropy from the matching-odd sector using self-matching/self-dual controls and derivative parity.

Prospective `N=185,265`, paired motif controls, finite-width annihilators, and kappa3 work continue in parallel when they do not displace those three tests.

## What this does not claim

Numerical estimates place the infinite square-site threshold near `0.59274605079`, but the project does not treat a rounded estimate as a definition and does not claim a known closed form. Published last digits are method-dependent and are tracked in `data/literature_threshold_sources.json`.

Likewise, the current orientation results do not yet prove:

- a unique asymptotic exponent;
- a unique H4 harmonic;
- an `x=21/4` LCFT operator identification;
- universality of the proposed derivative invariants;
- an exact expression for `p_c`.

Those are research questions, not merge blockers for exploratory code and data.

## Repository layout

```text
constants/      reference values and exact relations
data/           literature/source datasets and provenance
notes/          synthesis, theory, derivations, and negative results
scripts/        analysis, exact checks, and experiment tooling
experiments/    frozen protocols and compute queues
predictions/    preregistered numerical predictions
results/        raw/derived research result archives
tests/          smoke, regression, and exact-contract tests
docs/           status, roadmap, and project governance
src/            production C++ simulation engines
```

## Research workflow

This is an exploratory mathematics/computational-physics repository, not a production service. Useful code and result archives should enter `main` quickly once they are understandable and the relevant tests run. Stronger scientific language is controlled by evidence level rather than by keeping exploratory work off the main branch.

For local checks:

```bash
python3 -m unittest discover -s tests -p 'test_*.py' -v
```

See [`GOVERNANCE.md`](GOVERNANCE.md) for the lightweight evidence levels and [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) for how to preserve expensive run provenance.

## License

MIT. See [`LICENSE`](LICENSE).
