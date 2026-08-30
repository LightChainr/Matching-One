# P250 Z5 charged three-point stream

This note turns the exact charged-fusion prediction into one auditable
measurement.  It is an operator-algebra experiment, not another vote on the
ordinary global `A_top` channel.

## Observable and exact map

The parent is the Gaussian torus `(8+i)` of area 65.  The two children are
the same-parent hands `(8+i)(2+i)` and `(8+i)(2-i)`, both of area 325.  A
counter-derived uniform parent translation moves the fixed three-anchor
pattern

```
(0,0), (1,0), (0,1).
```

At every anchor and in each hand, the runner evaluates the existing local
primal-minus-matching landing-pivotal H4 row at all five deck fibers.  With
the transported generator fixed by the exact cover map,

```
O_r = (1/5) sum_k zeta5^(-r k) local_H4(k).
```

One common Bernoulli field therefore yields `C113` and `C122` in both hands,
their independently accumulated conjugates `C244` and `C334`, and the
nonneutral controls `C111` and `C112`.  The response saves the complete 8x8
real covariance of the four primary complex channels and the complete 24x24
joint covariance including controls.

The exact gate checks all 65 parent translations, all 1,950 child root labels,
parent projection, anchor noncollapse, and DFT conjugacy.  The random
translation changes location only; it does not refit the anchor shape or a
phase.

## Frozen score

The score order comes unchanged from
`predictions/z5_charged_ope_fusion_20260829.json`:

1. the two-real-dimensional cross-product closure
   `C113_plus*C122_minus-C113_minus*C122_plus=0`, with delete-one covariance;
2. the joint eight-real GLS for H4, H8 and H12 using the exact fixed
   `q_(3s)=((2+i)/(2-i))^(3s)` phases and only two complex amplitudes;
3. the joint nonneutral zero control and configurationwise DFT conjugacy.

No channel receives its own fitted hand phase.  The raw cubic phase is not
called universal, and this measurement alone does not construct the
normalization-free `I_rst`, whose neutral two-point factors are not in this
minimal stream.

## Engineering smoke and budget

The 20,000-replica smoke used the already frozen score solely to verify the
schema and estimate variance.  It took 194.69 seconds on eight local workers.
All exact gates and conjugacy controls passed.  The largest primary variance
of the mean was `4.568e-12`, corresponding to a component standard error of
`2.14e-6` at 20k.  Scaling the covariance without using the observed channel
means gives a worst-component target standard error of about `3.02e-7` at one
million replicas.

Accordingly the production acquisition is frozen at one million replicas,
100 batches and 16 workers.  The production seed differs from the smoke seed.
The anchor, DFT convention, p, radius, hypotheses, score order, and covariance
contract are unchanged.  The 20k point estimates are model-development data
only and are not combined with production.

### Execution-only amendment

Before production produced any response, the one-million-replica counter
domain was partitioned into three disjoint ranges of 340k, 330k and 330k.
Every range retains 10k batch units and uses the same runner commit, seed, p,
radius, anchor, DFT, and frozen scores.  The merge oracle verifies exact
coverage `[0,1000000)`, then recomputes all point estimates and covariance from
the 100 retained batch sufficient-statistic rows.  This is wall-clock
parallelism only; it creates no new scientific degree of freedom.
