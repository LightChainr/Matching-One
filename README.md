# Matching One

Matching One is a computational research project on square-lattice site percolation, its matching-lattice identity, and the finite-size structure behind the threshold.

The exact structural anchor is

\[
p_c^{\mathrm{site}}(\mathbb Z^2)+p_c^{\mathrm{site}}(\mathrm{NN+NNN})=1.
\]

## What we currently observe

The strongest numerical structure is an orientation-dependent matching signal on primitive Gaussian tori.

- Independent 100M same-`N` confirmation at `N=65,85,130,145,170` reproduces the prescribed odd-square-harmonic sign and is compatible with

  ```text
  DeltaM ~ DeltaCos4 * N^(-13/8).
  ```

- Three prospective `1+i` Gaussian lineages are compatible with the no-fit raw-contrast transformation

  ```text
  DeltaM(2N)/DeltaM(N) = -2^(-13/8).
  ```

- The genuinely new N=185/265 full-curve target block, 500M paired permutations per size, gives

  ```text
  DeltaM, x=21/4 H4-like: chi2 = 3.046 / 2
  DeltaM, zero:            chi2 = 29.409 / 2
  DeltaM, x=17/4:          chi2 = 30.246 / 2
  ```

  so the matching-odd `x=21/4` H4-like radial law survives the new-geometry test and strongly outperforms these two frozen alternatives.

### Important channel-map erratum

The original #108 report described the matching-even N=185/265 result as a prospective sign reversal. That interpretation was a protocol error, not a new physical effect.

The frozen matching-even source amplitude was fitted from P31 `either/even`, while the threshold-rank target reconstructs rank-2 `cross/even`. Complementary torus topology gives

```text
DeltaS_cross = -DeltaS_either.
```

After applying only this exact source-to-target channel map, with no refit,

```text
corrected DeltaS chi2 = 0.5700315 / 2
marginal residual z   = +0.67, -0.12
```

The matching-even prospective result is therefore compatible with the frozen amplitude after channel conversion. The raw #108 files remain preserved; the erratum/scorer is on `main` via PR #134.

This is an important methodological lesson for the whole repository: every frozen prediction must name its wrapping channel, and every scorer must make any exact channel map explicit.

### Prospective P48 parity score

The same N=185/265 full curves are also genuinely new targets for the intrinsic-center P48 amplitudes frozen from N=65,85,130. With zero target refits,

```text
P4[S]   ~ N^-1:     chi2 =  1.139 / 2
P4[D]   ~ N^-13/8:  chi2 =  0.281 / 2
P4[D']  ~ N^-5/8:   chi2 =  0.088 / 2
P4[S']  ~ N^-5/4:   chi2 = 52.716 / 2
```

So `S`, `D`, and `D'` pure laws all survive on prospective new geometries. `S'` is the unique clear pure-law failure among these four channels.

The predeclared corrections for `S'` remain viable:

```text
rank-2/log:    chi2 =  1.204 / 2
analytic 1/N:  chi2 =  0.862 / 2
zero:          chi2 = 1278.555 / 2
```

The correction mechanism is not yet uniquely identified.

## Working interpretation

A compact empirical law for the central matching-odd sector remains

\[
\Delta M_N \approx A\,\Delta\cos(4\theta)\,N^{-13/8}.
\]

It now has independent-seed, held-out, Gaussian-semigroup, root-closure, and prospective-new-geometry support. It does **not** yet prove that H4 is unique rather than H12/H20, nor uniquely identify an `x=21/4` LCFT operator.

The matching-even central sector is compatible with its frozen `N^-1` amplitude once cross/either semantics are aligned, and the independent intrinsic-center P48 score confirms pure `S`, `D`, and `D'` transfer on the new geometries. This strengthens the empirical parity picture, but does not prove a local matching/OPE automorphism.

Full-curve data also show that the bare finite-size center-slope ratio is not exactly `2^(3/8)` at current precision; a small but resolved finite-size correction is required. For normalized `P4=DeltaX/DeltaCos4`, the Gaussian angular factor is already divided out, so pure-H4 normalized transfer is positive `Q^(-alpha)` even when the raw contrast changes sign.

## Next experiments and analysis

The execution order is intentionally short:

1. **Norm-5 H4 versus H12 — #57.** Highest-information new-compute discriminator. N=325/425 are supported by the threshold-rank engine. Use a frozen variance/power pilot and preserve raw-versus-normalized channel conventions.
2. **Third full-curve lineage — #50.** Score `145 -> 290`, including the already-frozen finite-size slope correction and induced root prediction.
3. **Use existing full curves harder before simulating more.** Priority analyses include the S-prime/channel-coordinate program (#48), prequential evidence ledger (#95), pivotal/Russo bridge (#100), intrinsic quantile-center spectroscopy (#101), multi-u thermal response (#119), and joint operator-mixing treatment (#125).
4. **Choose later expensive geometry by information gain — #102.** Do not default to a larger N.

Parallel theory/control routes are tracked in the post-P43 frontier map: FK/Potts torus sectors, four-arm anisotropy, torus-modulus spectroscopy, exactly-critical isoradial controls, Euler/Betti identities, universal amplitude ratios, exact self-matching Beta-family tests, and finite-polynomial Galois certificates.

The exact N=1105 four-angle projector, axis-annihilator work, complex-zero maps, kappa3, rigorous-bound feasibility, and PSLQ are secondary/gated tracks.

## Repository

The numerical archive, production source, exact checks, failed/null models, raw sufficient statistics, and protocol errata live on `main`. The execution-facing synthesis is [`notes/SYNTHESIS-20260828.md`](notes/SYNTHESIS-20260828.md); the claim ledger is [`docs/STATUS.md`](docs/STATUS.md).

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

This is an exploratory mathematics/computational-physics repository. Useful research assets should enter `main` quickly; claim strength is controlled by evidence and chronology, not by branch location.

```bash
python3 -m unittest discover -s tests -p 'test_*.py' -v
```

No closed form for square-site `p_c` is claimed. Published numerical estimates remain method-specific and are tracked in `data/literature_threshold_sources.json`.

## License

MIT. See [`LICENSE`](LICENSE).
