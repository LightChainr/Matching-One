# Matching One

Matching One is a computational research project on square-lattice site percolation, its matching-lattice identity, and the finite-size structure behind the threshold.

The exact structural anchor is

\[
p_c^{\mathrm{site}}(\mathbb Z^2)+p_c^{\mathrm{site}}(\mathrm{NN+NNN})=1.
\]

## What we currently observe

The strongest numerical structure is an orientation-dependent sector on primitive Gaussian tori.

- Independent 100M same-`N` confirmation at `N=65,85,130,145,170` gives the prescribed `Delta cos(4 theta)` sign at all five sizes, with

  ```text
  A4 = N^(13/8) DeltaM / DeltaCos4 = 0.7885 +/- 0.0352
  chi2 = 1.53 / 4.
  ```

- Three prospective `1+i` Gaussian lineages are compatible with the no-fit raw-contrast law

  ```text
  DeltaM(2N)/DeltaM(N) = -2^(-13/8).
  ```

- The first genuinely new-geometry full-curve test, `N=185,265` with 500M paired permutations per size, gives

  ```text
  DeltaM, x=21/4 H4: chi2 = 3.046 / 2
  DeltaM, zero:       chi2 = 29.409 / 2
  DeltaM, x=17/4:     chi2 = 30.246 / 2
  ```

  so the matching-odd `x=21/4` H4-like radial law survives prospectively and clearly beats both zero and the tested x=17/4 adversary.

The same new geometries also provide a stronger test of the intrinsic-center P48 parity channels. Using only amplitudes frozen from `N=65,85,130`, with zero target refits:

```text
P4[S]   ~ N^-1:     chi2 =  1.139 / 2
P4[D]   ~ N^-13/8:  chi2 =  0.281 / 2
P4[D']  ~ N^-5/8:   chi2 =  0.088 / 2
P4[S']  ~ N^-5/4:   chi2 = 52.716 / 2
```

Thus `S`, `D`, and `D'` pure laws survive on prospective new geometries. `S'` is decisively nonzero but its pure law fails; two predeclared corrections remain viable:

```text
rank-2/log:    chi2 = 1.204 / 2
analytic 1/N:  chi2 = 0.862 / 2
zero:          chi2 = 1278.555 / 2
```

## Important protocol correction

The original Issue #43 registered matching-even `DeltaS` score used a P31 **`either/even`** source amplitude, while the threshold-rank target engine/scorer reconstructs rank-2 **`cross/even`**. The original positive registered prediction therefore fails (`chi2=240.247/2`) and remains preserved as a failed preregistration artifact.

But P31 already contains both channels and verifies the exact orientation-difference map

```text
DeltaS_cross = -DeltaS_either.
```

Applying only that pre-existing map, with zero target refits, gives a post-reveal protocol-repair diagnostic

```text
corrected cross/even chi2 = 0.570 / 2.
```

This is not a retroactive preregistered pass. It means the `240/2` failure is evidence of a source/target channel-contract bug, not evidence that the physical matching-even N^-1 sector is falsified. The audit and original score are both retained.

## Working interpretation

A compact empirical law for the central odd sector is

\[
\Delta M_N \approx A\,\Delta\cos(4\theta)\,N^{-13/8}.
\]

It now has independent-seed, held-out, Gaussian-semigroup, root-closure, and prospective-new-geometry support. It still does not prove that H4 is unique rather than H12/H20, nor uniquely identify an `x=21/4` LCFT operator.

The broader intrinsic-center parity pattern is also stronger than an early reading of the N=185/265 result suggested: `P4[S]`, `P4[D]`, and `P4[D']` all transfer successfully. The specific unresolved channel is `P4[S']`, where a finite-size correction is required.

Full-curve data additionally show that the bare finite-size slope ratio is not exactly `2^(3/8)` at current precision; a small but resolved correction is required. For normalized `P4=DeltaX/DeltaCos4`, the Gaussian angular sign cancels, so pure H4 normalized transfer is positive `Q^(-alpha)` even when the raw contrast changes sign.

## Next experiments

The execution order is intentionally short:

1. **Norm-5 H4 versus H12 — #57.** Highest-information dedicated discriminator. `N=325,425` are supported by the production engine. Start with a small threshold-rank variance pilot and choose production size from measured power.
2. **Third full-curve lineage — #50.** Score `145 -> 290`, including the frozen finite-size slope correction, induced root prediction, and derivative transfer.
3. **S-prime correction / coordinate mapping — #48.** Use existing data first to distinguish q=2 versus log/Jordan structure and relate fixed-coordinate cross/either observables to intrinsic-center P48 projectors. No dedicated new run yet.

The exact N=1105 projector, axis-annihilator work, complex-zero maps, kappa3, and PSLQ are secondary tracks and should not displace these discriminators.

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
