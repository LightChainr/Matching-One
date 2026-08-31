# Matching One

Matching One is a computational research project on square-lattice site percolation, its matching-lattice identity, and the finite-size structure behind the threshold.

The exact structural anchor is

\[
p_c^{\mathrm{site}}(\mathbb Z^2)+p_c^{\mathrm{site}}(\mathrm{NN+NNN})=1.
\]

## Current scientific boundary

The strongest numerical result is a reproduced matching-odd orientation signal on square-site Gaussian tori. Independent primary blocks disfavor global zero and remain compatible with the frozen H4 / `N^-13/8` finite-size prediction. The norm-5 child block discriminates H4 from tested H12/H8 aliases but, by itself, remains compatible with zero.

Several simpler stories are already retired: the N145->290 full curve is not one scalar multiplier; pure `P4[S'] ~ N^-5/4` is insufficient; one scalar width does not close the higher thermal jet. q=2/Jordan is not uniquely identified, and primitive square-bond `x≈4` spin-4 behavior remains a distinct sector from the square-site thermal-Q4 `x=21/4` candidate.

The project does **not** claim a closed form for square-site `p_c` or a uniquely identified continuum operator.

## Current decision

**P0 production is empty.** The requested independent tests are complete; their failed parameterizations are stopped, while the broader #154/#334/#337 questions remain open at P1.

- **#154:** the [165M-permutation result](https://github.com/LightChainr/Matching-One/blob/f4999e29612da16a3650f24d124fb59137f053d7/experiments/p154-prospective-transmission-20260831/REPORT.md) rejects the frozen strong entry/completion predictions and satisfies the weak net-response stop rule. Do not prioritize this lag1 policy as the main H4 explanation or top up its sample.
- **#334:** the [1M-prefix source-normal intervention](https://github.com/LightChainr/Matching-One/blob/d0a9daf1132779205f119e9b4470f4eea9cb89c1/notes/p334-independent-normal-intervention-result.md) rejects complete two-score conditional-mean closure; the separate [600k-prefix contact intervention](https://github.com/LightChainr/Matching-One/blob/14b2c98ed3a252a2fe79ce5e124d9484b23a264f/experiments/p334-prospective-intervention-20260831/REPORT.md) rejects both frozen residual forecasts. No same-block feature rescue or new half-amplitude law follows.
- **#337:** [finite two-coupling profile closure is excluded](https://github.com/LightChainr/Matching-One/blob/f5aa94c3c1d619da2272f2409623ed52c876463d/experiments/p337-two-coupling-closure-20260831/results/REPORT.md). The [2690f665 uniform bounds](https://github.com/LightChainr/Matching-One/blob/2690f665bc8029cb2370d3f1efcef5eb2853705c/notes/closed-source-uniform-projection-tail.md) prove opposite original-U signs for Sstar/Sdrop at fixed N25 for every real m>=64; the [Poisson joint limit](https://github.com/LightChainr/Matching-One/blob/2690f665bc8029cb2370d3f1efcef5eb2853705c/notes/closed-source-poisson-double-scaling.md) suppresses U when N/m² stays bounded on growing-systole quotients. The remaining fixed-m question concerns actual phase/sector control; the local-pair spatial question is stated below. Neither opens an automatic sampling window. These branch-delivered assets are not represented as merged.

**Further completed finite calculations:** [homogeneous N50](https://github.com/LightChainr/Matching-One/blob/ef3b2c68f824e29421747c805ea7a505aca41908/experiments/p337-homogeneous-n50-20260831/RESULT.md) now gives original U=1.0615603877 and V_S=+0.0543457827, with full 2^50 counts per geometry and independent derivative review. The finite zero-transmission null is excluded; positive-sign continuation survives without mechanism confirmation. This contract is finished, with no N100/t/epsilon extension. The separate [specified Q1 closed trace](https://github.com/LightChainr/Matching-One/blob/bea717e826df5a22518774b1725ae7bcbe2cb801/notes/p337-q1-closed-trace-transmission-result.md) also has completed finite landing and V=-0.0019048362. The [canonical local pair](https://github.com/LightChainr/Matching-One/blob/2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb/notes/regular-pair-interaction-result.md) has also completed finite original-U transmission: the direct Q1 response is identically zero, while the fixed mixed response is -0.04503611398. The [new spatial support result](https://github.com/LightChainr/Matching-One/blob/baa5d33b2f87b2868aa0cb9d3f6518c93dbf3bff/experiments/p337-regular-spatial-support-20260901/RESULT.md) excludes first-Q transmission through at most one shared component and bounds the canonical annealed response by (43/16) times the two-component event probability. All 4140 fixed-kernel entries agree between two independent exact algorithms. Remaining work concerns the weighted spatial probabilities/scaling and the same original-U response, not another finite single-mark score. No new production follows.

The exact birth-current organization is

```text
d_p A_ell = j_in,ell - j_out,ell
J_grad    = j_in - j_out
J_act     = j_in + j_out
```

A strong local `J_act` response is not automatically a mechanism for global `U`. The validation rule is deliberately hard: if a frozen model fails, it is removed or demoted; the same block is not used to add another rescue feature.

## Research navigation

There are three current navigation documents, in this order:

1. [`docs/STATUS.md`](docs/STATUS.md) — the sole authoritative current claim ledger;
2. [`docs/ROADMAP.md`](docs/ROADMAP.md) — current information-gain priorities and stop rules;
3. [`docs/RESEARCH-MAP.md`](docs/RESEARCH-MAP.md) — durable scientific-track map.

Supporting machine/provenance indexes:

- [`analysis/research_ledger.yaml`](analysis/research_ledger.yaml) — current machine-readable research/decision state;
- [`analysis/artifact_registry.yaml`](analysis/artifact_registry.yaml) — canonical and unmerged-asset pointers;
- [`analysis/evidence_ledger_manifest.yaml`](analysis/evidence_ledger_manifest.yaml) — raw-data groups, primary/sensitivity roles and freeze chronology;
- [`docs/CLEANUP-20260831.md`](docs/CLEANUP-20260831.md) — the 2026-08-31 issue/PR cleanup and evidence-dedup audit.

Branch-local `NEXT-TARGETS`, Draft PR bodies, old synthesis notes and closed Issue handoffs remain provenance, not competing current-status documents.

## Execution constraints

Claim-bearing work keeps three hard constraints:

1. do not rewrite frozen predictions or committed result history;
2. do not silently compare incompatible observable semantics;
3. do not count correlated views of one raw random block as independent primary evidence.

The repository test command is

```bash
python3 -m unittest discover -s tests -p 'test_*.py' -v
```

See [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) for provenance conventions and [`GOVERNANCE.md`](GOVERNANCE.md) for the broader exploratory execution policy.

## License

MIT. See `LICENSE`.
