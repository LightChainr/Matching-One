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

Only two issues are live for the 2026-08-31 decision:

- **#334:** execute the already frozen independent projective birth-current intervention;
- **#154:** freeze and then test a prospective source -> ingress/egress -> `J_grad` -> original pooled-root/slope-normalized `U` transmission.

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
