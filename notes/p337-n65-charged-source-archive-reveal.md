# Existing-archive reveal of the frozen F3 charged sources

Status: a no-new-simulation reveal from the aligned N65 projective-birth archive
at commit `1714141`.  The source definitions were frozen in `539b629`:

```text
q_A=T_01-T_10,  q_D=T_12-T_11.
```

No line combination, phase, reference `p`, random block, or covariance grouping
was selected after inspecting this reveal.

## What the archive identifies

For every archived `(tau1,ell,tau2)` cell, reduce the primitive line `ell`
modulo 3.  At fixed `p`, its rank-one plateau contributes to

```text
W_C(p)=E[q_C^2],
J_C,birth(p)=E[flux(tau1) q_C^2],
J_C,exit(p)=E[flux(tau2) q_C^2].
```

Thus the complete frozen six-vector is identifiable without a new field:

```text
(W_A,J_A,birth,J_A,exit,W_D,J_D,birth,J_D,exit).
```

All 20 aligned batches for both N65 orientations are retained in one covariance
matrix.  The two shapes use the same counter stream, so their difference is
formed within batch.

## Charged activation is clear

At `p_ref=0.592746050790`:

| orientation | sector | `W_C` | SE | z | birth | exit |
|---|---|---:|---:|---:|---:|---:|
| `8+i` / first | A/B1 | 0.336204267896 | 0.00204844 | 164.13 | 3.42110897046 | 3.47899626478 |
| `8+i` / first | D/B2 | 0.0370407065872 | 0.00103950 | 35.63 | 0.640526567485 | 0.627038212802 |
| `7+4i` / second | A/B1 | 0.338344544086 | 0.00189693 | 178.36 | 3.48303234796 | 3.44390381009 |
| `7+4i` / second | D/B2 | 0.0369750143646 | 0.000644550 | 57.37 | 0.628804925877 | 0.625908683408 |

The explicit A and D source responses are therefore numerically easy to
resolve even in the old 20k engineering block.  This does not contradict the
unweighted A/D null: the measured quantity is the even susceptibility
`E[q_C^2]`, not the forbidden one-point `E[q_C]`.

The F3 phase is fixed, not fitted:

```text
O_C(omega)=(omega-omega^2)W_C/2.
```

Because all four `W_C` estimates are positive, both charged one-points lie on
the same `+i` ray.  There is no empirical phase choice.

## Null and current controls

The exact response matrix explains which controls should and should not vanish:

```text
d<H>/ds_A = E[H A] proportional to E[A] = 0 in the ensemble,
d<H>/ds_D = E[H D] proportional to E[D] = 0 in the ensemble,
d<A>/ds_D = d<D>/ds_A = E[A D] = 0 statewise.
```

The A-D cross entries are exactly zero in every archive row.  The H cross
entries reduce to the sampled unweighted A/D coordinates and fluctuate because
the archive did not antithetically pair each path with its quarter-turn.  Their
two-dimensional null quadratics are `4.510/2 df` for first and `3.207/2 df` for
second, compatible with their role as noise controls rather than charged
responses.

For both sources, direct differentiation of the plateau and the source/sink
construction agree to at most `4.1e-15`:

```text
dW_C/dp=J_C,birth-J_C,exit.
```

## Same-N orientation response

Source activation and orientation modulation are different questions.  The
within-batch second-minus-first scores are:

```text
A (W,birth,exit): 12.153 / 3 df,
D (W,birth,exit):  1.509 / 3 df,
joint six-vector: 15.530 / 6 df.
```

The A current triplet carries a visible orientation response, although its
plateau `W_A` difference alone is only `0.900 sigma`; the information is in the
joint source/exit timing.  D orientation modulation is not resolved.  These
are same-block descriptive scores from an engineering archive, not a new
independent discovery block.

## Representation transport

The archive also closes the algebraic representation checks:

- quarter-turn sends both `q_A` and `q_D` to their negatives, so each has
  projective C4 charge 2;
- reflection leaves `q_A` even (B1) and makes `q_D` odd (B2);
- direct permutation of all four F3 line bins by `T` agrees with the frozen
  H/A/D shear matrix for both the state vector and its response matrix, with
  maximum residual `4.2e-17`.

This last item is an internal relabeling certificate.  The `8+i` and `7+4i`
archives are two microscopic Gaussian quotients, not an identity/T-shear
source pair, so their difference is not advertised as an independent shear
experiment.

## Boundary and reproduction

The reveal establishes finite-volume charged-source identifiability and shows
that the existing archive carries a nontrivial A current timing contrast.  It
does not establish large-N survival, a continuum operator, or a new shear
experiment.

```bash
python3 scripts/reveal_n65_charged_source_archive.py \
  --births /path/to/n65_20k.births.csv \
  --metadata /path/to/n65_20k.metadata.json \
  --births-label results/local-20260830/P334-projective-birth-N65-smoke/n65_20k.births.csv@1714141 \
  --metadata-label results/local-20260830/P334-projective-birth-N65-smoke/n65_20k.metadata.json@1714141 \
  --json results/p337-n65-charged-source-reveal/latest.json \
  --markdown results/p337-n65-charged-source-reveal/latest.md

python3 -m unittest discover -s tests \
  -p 'test_n65_charged_source_archive_reveal.py'
```
