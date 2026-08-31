# P334: shared next-label birth gates and the complete even clock response

Status: fixed thin analysis **before** new quartet counts are read. The exact
long-CSV schema is supplied by producer `a3249a59`; the producer will separately
supply the immutable completed source identity. The reader takes only that
commit and requires all forty gzip batches. It does not consume partial
production or start sampling.

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
  --source-commit PRODUCER_CONFIRMED_COMPLETE_SHA
```
