# P334: shared next-label birth gates and the complete even clock response

Status: completed. Formula freeze `47d6eb41` and exact-schema reader freeze
`874b1025` preceded target counts. Producer `a3249a59` completed the forty
gzip batches in immutable source `e32a85939279b8574278024d647b56d2d1485247`.
The following fixed definitions were scored once on that complete source.
This analysis started no sampling; it uses the producer's1,280,000 fresh
conditional tails on the original40,000 paired prefixes.

## One within-prefix spatial question

In original checkpoint cells01 and10, name the lower-rank orientation R0 and
the higher-rank orientation R1. For a common next insertion label u define

```
g(u) = (R0 gains its first ambient-H1 birth,
        R1 gains its second ambient-H1 birth).
```

The first entry is `next_rank>=1`, including a0->2 simultaneous birth. The
second is `next_rank==2`. Thus01/10 use the same named gate roles even though
their old first/second orientation order is reversed.

The question is whether a common position promotes the two births together
or puts them in competition **within the same fixed prefix**, and how this
spatial contrast aligns with its remaining complete A/E response.

## Unbiased half-difference estimators

Let U,V be independent common next labels drawn from that prefix's remaining
labels. Their positions need not differ. Conditional on each label, draw two
independent suffixes; the two original orientations share each suffix order.
For any two-vector X, define m_U as the average of its two suffix outcomes.

```
G_q = (g(U)-g(V))(g(U)-g(V))^T / 2,
C_q = (g(U)-g(V))(m_U-m_V)^T / 2.
```

At a fixed prefix, these have expectations `Cov_U(g,g)` and
`Cov_U(g,E[X|U])`, respectively. Suffix noise drops out of the second formula
because the gate is fixed by the label. No squared noisy conditional mean is
used here. A negative off-diagonal G entry indicates competition; a positive
one indicates co-promotion relative to the conditional gate rates.

For clarity, the saved same-label joint gate rate and cross covariance obey
the expectation identity

```
E G_01 = E_U[g0(U)g1(U)] - E_U[g0(U)] E_U[g1(U)].
```

The second term is a product **inside each prefix**. Substituting the product
of population-averaged rates would insert between-prefix heterogeneity and
answer a different question.

## Two nonredundant clock coordinates

Use complete integrated A/E:

```
A=1-(K1+K2)/(N+1),
E=1-(K2-K1)/(N+1),
X=(A_first-A_second, E_first-E_second)/delta_cos4.
```

This retains the old physical orientation sign and H4 normalization; it does
not identify the gate with a harmonic. K1,K2 andW need not be added as redundant
columns. The second X coordinate directly addresses the full even response
already identified with the lifetime difference.

Keep01 and10 separate. Their sum describes gate covariance against the old
global orientation contrast. For a **raw unscaled** role-aligned low-minus-high
observer, the same saved vectors give `delta_cos4*(01-minus10)`; this removes
the old normalization and is not a new H4 coefficient. Neither combination is
chosen after looking at signs. The delta values remain those in `bb79fd47`:
N325=-0.7634556213017751 and N425=-0.8928996539792388.

## Batch handoff and interpretation

Average the eight quartets within every original prefix. For each cell, save
the13 fixed fields in `scripts/p334_next_label_gate_coupling.py`, with zero
contribution from all other prefixes. Divide by each **original1000-prefix
batch**, not the number in01/10. This retains full20k denominators and cell
prevalence. The two cells'26-vector is handed to the common covariance
coordinator on the original20 batch IDs; no independence between cells,
readouts, providers or old/new uses of the same prefixes is asserted.

The covariance is conditional on the existing common-random-number label
coupling and the selected geometry pair. It is not an intrinsic or universal
covariance of space. These cross moments do not determine a reliable fraction
of full Doob information without the needed conditional gate probabilities
and response covariance; no such ratio will be substituted.

The completed-source command is:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
/Users/lc/python-envs/research-py311/bin/python scripts/p334_next_label_gate_coupling.py \
  --source-commit e32a85939279b8574278024d647b56d2d1485247
```

## New spatial result: local alignment and transpose cancellation coexist

Every number below retains the full20k denominator. For a cell c, a covariance
entry estimates `E_prefix[1_c Cov_U(...|prefix)]`; it is not a population-rate
covariance or a statement about the sign at every individual prefix.

| Joint cell / average conditional gate cross-covariance | N325 | N425 |
| --- | ---: | ---: |
| 01 prevalence | 2406/20000=12.030% | 2550/20000=12.750% |
| 10 prevalence | 2457/20000=12.285% | 2402/20000=12.010% |
| 01 gate cross-covariance | 0.00020625 +/-0.00007355 | 0.000053125 +/-0.00005503 |
| 10 gate cross-covariance | 0.00017500 +/-0.00005241 | 0.000006250 +/-0.00004202 |
| 01+10 gate cross-covariance | **0.00038125 +/-0.00009521** | **0.000059375 +/-0.00007820** |

The N325 selected population shows positive average conditional co-promotion
at about4 shared-batch standard errors, with positive point estimates in both
transpose cells. N425's average sign is unresolved. This is evidence about the
specified shared-label coupling, not proof of a universal geometric alignment
or an inferred cross-size exponent.

The clearer clock-response structure appears after applying the **previously
specified** low-minus-high, unscaled projection. Entries are gate-row ×(A,E):

| N / fixed gate role | Cov with low-minus-high integral A | Cov with low-minus-high integral E |
| --- | ---: | ---: |
| 325 R0-first | +0.000692863 +/-0.00001209 | -0.000149688 +/-0.00000931 |
| 325 R1-second | -0.000386043 +/-0.00001226 | -0.000403269 +/-0.00001040 |
| 425 R0-first | +0.000543358 +/-0.00001357 | -0.000110303 +/-0.00000586 |
| 425 R1-second | -0.000319876 +/-0.00000998 | -0.000317499 +/-0.00000937 |

**Both gate roles align with a decrease of the low-minus-high E difference,
while their A directions oppose.** The E magnitude attached to the higher-rank
second-birth gate is larger. This gives a concrete two-birth response geometry
using two nonredundant complete clock readouts, not a fitted latent model.

That local alignment does not automatically explain the global E contrast.
Restoring the original first-minus-second/H4 convention sums opposite-signed
01 and10 responses. The remaining gate×E entries are

```
N325: R0 -1.71577e-5 +/-1.03317e-5; R1 +9.83129e-6 +/-1.53318e-5
N425: R0 +3.79970e-6 +/-1.07966e-5; R1 -1.17893e-6 +/-1.17608e-5.
```

None resolves the global gate×E alignment at two standard errors. In contrast,
the N425 R0×A residual is `-5.03327e-5 +/-1.18803e-5`, approximately4.24 standard
errors. The observed transpose asymmetry is sharper in this A coordinate;
selecting a new mechanism after these signs is not part of this calculation.

## R0 can jump directly to R2

The prescribed R0-first gate included all `next_rank>=1` from the outset. The
existing next-rank field also gives the following descriptive counts of direct
0->2 jumps. A next-label draw is counted once, not once per suffix; independent
U andV draws remain two draws even if they happen to choose the same label.

| N / cell | R0 next-rank2 / same-cell label draws | Frequency |
| --- | ---: | ---: |
| 325 /01 | 26/38496 | 0.0675395% |
| 325 /10 | 26/39312 | 0.0661376% |
| 425 /01 | 16/40800 | 0.0392157% |
| 425 /10 | 14/38432 | 0.0364280% |

These retain the simultaneous-birth support without adding a candidate model
or another covariance block. They are raw repeated-next-label frequencies on
the fixed prefix population, not independent trials across all suffix rows.

## Shared handoff

`results/p334-next-label-gate-coupling/batch_vectors.json` contains the fixed
26-vector for each original batch, its matching delete-one vectors, original
cell counts, raw0->2 counts, both fixed orientation projections and all40
input hashes. The covariance coordinator received these vectors directly,
before a scientific conclusion was independently aggregated elsewhere. The
displayed standard errors come from those same20 paired batches; no covariance
inverse or extra evidence block is introduced.
