# Independent N85 reveal of the natural A-current coordinate

Status: independent fresh N85 scale point, generated only after the complete
normalization, geometry, sample size, seed and three model targets were pushed
in preregistration commit `5e601dc`.

## Why this coordinate

The N65 activity/net decomposition showed that the charged A orientation
contrast lives in the derivative current, while the plateau amplitude and
common activity do not carry it separately.  The amplitude-free natural
coordinate was therefore frozen as

```text
K_A = p(1-p) Jminus_A/W_A = d_eta log W_A,
eta = log[p/(1-p)].
```

This ratio is a logarithmic response, so an overall source normalization of
`W_A` cancels.  No N85 result is used to choose the denominator, `p`, exponent,
orientation order, or source.

## Archive search and fixed design

No existing N85 or larger file contains the required
`(tau1,primitive ell,tau2)` sufficient statistics.  The older C01 and P34 N85
metadata have neither `projective_births` nor a `births_csv`.  They cannot
reconstruct the charged plateau and its projective source/sink current.

The smallest next same-semantic pair was frozen as

```text
N=85: first=9+2i, second=7+6i,
200000 samples/shape, 20 aligned batches,
seed=202608337, replica offset=0.
```

N65 variance projected `SE(Delta_K)=0.00644`; 200k exceeds the N65-only
three-sigma requirement for separating scale-neutral and project-H4 numeric
targets.  The observed N85 SE is `0.00793`, 1.231 times the projection.

The block ran on Huawei `DevEnvC_ZyTrST`
(`f415a4bcbd9a438b85f5f29e4a507ea4`, AArch64, 16 vCPU) using the exact
`1714141` projective-birth source plus freeze `5e601dc`.  Engine time was
`0.3697 s`; self-test passed and stderr is empty.

## Independent result

At the unchanged `p_ref`:

```text
N85 first  (9+2i): W_A=0.337113008108,
                    Jminus_A=-0.014785600461,
                    K_A=-0.010587597656;

N85 second (7+6i): W_A=0.335492737739,
                    Jminus_A=+0.020853818557,
                    K_A=+0.015005015228.
```

Therefore

```text
Delta_K_A(N85)=0.025592612884 +/- 0.007926872163,
z=3.229.
```

The natural orientation response remains independently nonzero at N85, but it
has attenuated from the N65 value `0.069480659660 +/- 0.020366733052`.
The descriptive ratio is `0.368 +/- 0.157`.  This ratio was not used to fit or
select an exponent.

## Frozen model comparison

The three N85 numeric targets were committed before production:

| target | frozen value | N85 residual | measurement quadratic | predictive quadratic |
|---|---:|---:|---:|---:|
| zero | 0 | +0.0255926 | 10.424 | 10.424 |
| source-fitted scale-neutral | 0.0694807 | -0.0438880 | 30.654 | 4.033 |
| source-fitted project H4 | 0.0449306 | -0.0193380 | 5.951 | 1.583 |

The measurement-only column treats each frozen number as exact.  The predictive
column adds N65 fit variance for the two source-fitted transfers, as declared
in the preregistration.

Three statements follow without a post-N85 fit:

1. Zero is disfavored by the independent N85 block (`3.229 sigma`).
2. Scale-neutral transfer is much too large; even after N65 fit uncertainty it
   has predictive quadratic `4.033`.
3. The observed attenuation is stronger than the nominal `N^-13/8` H4 target,
   but H4 remains the closest and predictively compatible frozen transfer
   (`1.583`).

The correct research update is therefore not “H4 scaling confirmed”.  It is:
the charged natural current persists at a second scale, rules out no-decay as
the leading description of this pair, and narrows the live region to H4-or-
faster attenuation.  A third independently frozen scale, not a fitted exponent
through these two points, is the next clean discriminator.

## Exact and dependency gates

- production metadata matches every preregistered geometry/RNG/sample field;
- `dW_A/dp=J_A,birth-J_A,exit` closes to `7.9e-15`;
- raw source and metadata SHA256 hashes are stored in the score and checksum
  ledger;
- N85 uses a fresh seed/counter block and is independent of the N65 archive;
- N65 enters only through the three precommitted transfer equations and their
  fit variance.

## Boundary and reproduction

This is one independent scale discrimination.  It does not estimate an
exponent, select a new normalization, or establish a continuum field.

```bash
python3 scripts/score_natural_current_scale.py \
  --preregistration analysis/p337_natural_current_scale_preregistration.json \
  --births results/server-20260830/P337-natural-current-scale-N85/raw/n85_200k.births.csv \
  --metadata results/server-20260830/P337-natural-current-scale-N85/raw/n85_200k.metadata.json \
  --json results/server-20260830/P337-natural-current-scale-N85/score.json \
  --markdown results/server-20260830/P337-natural-current-scale-N85/REPORT.md

python3 -m unittest discover -s tests \
  -p 'test_score_natural_current_scale.py'
```
