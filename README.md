# Matching One

Matching One is a computational research project on square-lattice site percolation, its matching-lattice identity, and the finite-size structure behind the threshold.

The exact structural anchor is

\[
p_c^{\mathrm{site}}(\mathbb Z^2)+p_c^{\mathrm{site}}(\mathrm{NN+NNN})=1.
\]

## What we currently observe

The strongest surviving numerical structure is a **matching-odd orientation sector on primitive Gaussian tori**.

- Independent 100M same-`N` confirmation at `N=65,85,130,145,170` gives the prescribed `Delta cos(4 theta)` sign at all five sizes, with z-scores `16.03, 11.23, 5.22, 5.27, 2.58` and pooled

  ```text
  A4 = N^(13/8) DeltaM / DeltaCos4 = 0.7885 +/- 0.0352.
  ```

- Three prospective `1+i` Gaussian lineages are compatible with the no-fit raw-contrast law

  ```text
  DeltaM(2N)/DeltaM(N) = -2^(-13/8).
  ```

- The first genuinely new-geometry full-curve test, `N=185,265` with 500M paired permutations per size, gives

  ```text
  DeltaM, x=21/4 H4: chi2 = 3.046 / 2
  DeltaM, zero:      chi2 = 29.409 / 2
  DeltaM, x=17/4:    chi2 = 30.246 / 2
  ```

  so the matching-odd `x=21/4` H4-like radial law survives prospectively and clearly beats both zero and the larger lower-dimensional adversary over these targets.

The same prospective run also **breaks the old simple two-sector picture**. The frozen matching-even `DeltaS ~ +N^-1` prediction reverses sign at both new sizes:

```text
N=185 DeltaS = -6.08154e-5 +/- 8.08957e-6
N=265 DeltaS = -7.02495e-5 +/- 9.38562e-6
frozen positive-law chi2 = 240.247 / 2
```

Likewise, the pure normalized derivative law `P4[S'] ~ N^-5/4` fails prospectively, while two predeclared corrections survive:

```text
pure N^-5/4:  chi2 = 52.716 / 2
rank-2/log:    chi2 =  1.204 / 2
analytic 1/N:  chi2 =  0.862 / 2
zero:          chi2 = 1278.555 / 2
```

The current picture is therefore **not** one clean two-field pure-power model. It is a robust matching-odd leading sector plus unresolved finite-size structure in the matching-even and derivative channels.

## Working interpretation

A compact empirical law for the surviving central odd sector is

\[
\Delta M_N \approx A\,\Delta\cos(4\theta)\,N^{-13/8}.
\]

This now has independent-seed, held-out, Gaussian-semigroup, root-closure, and prospective-new-geometry support. It still does not prove that H4 is unique rather than H12/H20, nor uniquely identify an `x=21/4` LCFT operator.

Full-curve data also show that the bare finite-size slope ratio is not exactly `2^(3/8)` at current precision; a small but resolved correction is required. For normalized `P4=DeltaX/DeltaCos4`, the Gaussian angular sign cancels, so pure H4 normalized transfer is positive `Q^(-alpha)` even when the raw contrast changes sign.

## Next experiments

The execution order is now intentionally short:

1. **Norm-5 H4 versus H12 — #57.** This is the highest-information dedicated discriminator. `N=325,425` are already supported by the production engine. Start with a small threshold-rank variance pilot, then choose production size from measured power rather than assuming billions of replicas.
2. **Third full-curve lineage — #50.** Score `145 -> 290`, including the already-frozen finite-size slope correction and induced root prediction.
3. **Even/derivative correction analysis — #48.** First use existing `N=65..265` data to explain the `DeltaS` sign reversal and distinguish corrected `S'` models. Do not schedule a dedicated new run until #50/#57 have been reused.

The exact N=1105 four-angle projector, axis-annihilator work, complex-zero maps, kappa3, and PSLQ are useful secondary tracks but should not displace these discriminators.

## Repository

The numerical archive, production source, exact checks, failed/null models, and raw sufficient statistics are on `main` under `results/`. The execution-facing synthesis is [`notes/SYNTHESIS-20260828.md`](notes/SYNTHESIS-20260828.md); the claim ledger is [`docs/STATUS.md`](docs/STATUS.md).

```text
constants/      reference values and exact relations
data/           literature datasets and provenance
notes/          synthesis, theory, derivations, negative results
scripts/        analysis and exact checks
experiments/    frozen protocols
predictions/    preregistered predictions
results/        raw and derived research archives
tests/          smoke, regression, exact-contract tests
docs/           status and roadmap
src/            production C++ engines
```

This is an exploratory mathematics/computational-physics repository. Useful research assets should enter `main` quickly; claim strength is controlled by evidence, not by hiding exploratory work on branches.

```bash
python3 -m unittest discover -s tests -p 'test_*.py' -v
```

No closed form for square-site `p_c` is claimed. Published numerical estimates remain method-specific and are tracked in `data/literature_threshold_sources.json`.

## License

MIT. See [`LICENSE`](LICENSE).
