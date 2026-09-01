# P337 thermal/pivotal preflight descriptive score

Status: **completed descriptive preflight only**.

All 64 frozen counters, 32 original pairs per configuration and 5,242,880 pair/site callbacks passed the structural checks. This report contains no significance test, ensemble centering, moving-root reconstruction, full J2, field ratio or independent evidence claim.

## Validation

- configuration rows: `64`
- shell rows: `19464`
- pair/site callbacks: `5242880`
- every configuration/pair partitions exactly N sites across complete carrier masks
- all q/E midpoint residual sums and maxima are zero
- pair/site absolute values were checked before descriptive aggregation

## Signed and absolute primitives by L

| L | sites | Δg signed/abs | q observable signed/abs | q kernel signed/abs | E observable signed/abs | E kernel signed/abs |
|---:|---:|---:|---:|---:|---:|---:|
| 32 | 1048576 | -1.625 / 2.875 | +0.125 / 0.125 | -0.25 / 1.25 | +0.125 / 0.125 | +0 / 1.25 |
| 64 | 4194304 | +1.75 / 1.75 | +0 / 0 | -1.75 / 1.75 | +0 / 0 | +1.75 / 1.75 |

The JSON preserves the same fields by dyadic shell, complete carrier mask, nonexclusive carrier bit, and endpoint/square-NN/external relation. Carrier-bit views overlap by construction and must not be added as independent or exhaustive votes.

## Interpretation boundary

These are uncentered finite callback sums over a deterministic subset replay. They do not estimate the centered observable-pivot or kernel-pivot expectation, and they do not apply the baseline q/E jets, common root, denominator or slope terms. No absence or sign in this 64-counter report authorizes more counters or a new seed.

Frozen scorer CLI:

```bash
python3 scripts/analyze_p337_thermal_pivotal_preflight.py \
  --config results/p337-thermal-pivotal-preflight/raw/preflight.config.csv \
  --shell results/p337-thermal-pivotal-preflight/raw/preflight.shell.csv \
  --metadata results/p337-thermal-pivotal-preflight/raw/preflight.metadata.json \
  --output-json results/p337-thermal-pivotal-preflight/latest.json \
  --output-md results/p337-thermal-pivotal-preflight/REPORT.md
```
