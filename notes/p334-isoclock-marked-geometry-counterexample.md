# P334: equal complete birth clocks do not identify the trigger geometry

The five-vertex double-star factor found in the actual P334 residual pair core
has the same independence polynomial as a four-cycle plus one isolated inert
site. Their **entire unmarked pair-birth clocks are identical**, but their
final insertion marks differ under every relabeling. The most economical
scalar discriminator is the probability that two independent continuations
of the same fixed state finish at the same site: **5/18 versus1/4**.

This is a constructed exact mechanism counterexample motivated by a real
P334 factor. It is not a new N425 checkpoint or a claimed realization of the
alternative graph on the specific square torus.

## The fixed two graphs

```
Double star G:             Cycle plus isolate H:

    2                     0---1       4 (inert)
    |                     |   |
3---0---1---4              3---2
```

G has edges01,02,03,14. H has edges01,12,23,30 and an isolated site4. Sites are
inserted once in a uniformly random order. Birth is the first insertion that
completes any graph edge. The isolate remains one of the five insertion labels.
These are pair-trigger graphs; their cycles are not ambient torus winding
cycles. Both are bipartite and define ordinary two-terminal connection events
by adjoining s to one bipartition and t to the other.

The independence polynomials obey the paper-level identity

\[
I_G(z)=1+5z+6z^2+2z^3
      =(1+z)(1+4z+2z^2)=I_H(z).
\]

Because `P(T>k)=[z^k]I/C(5,k)`, both have

\[
P(T=2)=\frac25,\qquad P(T=3)=\frac25,\qquad
P(T=4)=\frac15,\qquad E[T]=\frac{14}{5}.
\]

Every unmarked canonical occupancy response obtained from the same five-site
safe coefficients also agrees for every p. Additional accuracy in this
unmarked clock cannot decide tree versus cycle-plus-inert trigger geometry.

## Exact time-resolved final-site mark

V is the actual site inserted at the first birth. The following entries count
sequences among the120 equally likely insertion orders. They are obtained by
examining only32 subsets per graph: a safe preceding (k−1)-subset that becomes
unsafe on adding v contributes `(k−1)!(5−k)!` sequences to `(T=k,V=v)`.

| Graph / vertex role | Count T2 | Count T3 | Count T4 | π_v=P(V=v) | Inert-knockout mean increase |
| --- | ---: | ---: | ---: | ---: | ---: |
| G:0, degree3 hub | 18 | 20 | 12 | 5/12 | 6/5 |
| G:1, degree2 hub | 12 | 12 | 6 | 1/4 | 7/10 |
| G:2, leaf at hub0 | 6 | 4 | 0 | 1/12 | 1/5 |
| G:3, leaf at hub0 | 6 | 4 | 0 | 1/12 | 1/5 |
| G:4, leaf at hub1 | 6 | 8 | 6 | 1/6 | 1/2 |
| H:any one cycle vertex | 12 | 12 | 6 | 1/4 | 7/10 |
| H:4, inert isolate | 0 | 0 | 0 | 0 | 0 |

The unmarked column totals are48,48,24 in either graph. But the π multisets
are `{5/12,1/4,1/6,1/12,1/12}` and `{1/4,1/4,1/4,1/4,0}`. No relabeling can
remove this difference. In particular:

| Relabel-invariant readout | Double star | C4 plus isolate |
| --- | ---: | ---: |
| Positive pivotal support | 5 | 4 |
| Σ_v π_v² | **5/18** | **1/4** |
| Normalized mean-knockout concentration | 113/392 | 1/4 |

For the last row, let Δ_v be the mean-clock increase when v is permanently
blocked but kept as an inert insertion label. The previously derived identity
`Δ_v=Σ_k k P(T=k,V=v)` gives Σ_v Δ_v=E[T]=14/5. The reported concentration is
`Σ_v (Δ_v/E[T])²`, without any additional blockade simulation.

## Cheapest mark that breaks the degeneracy

With two independent fresh continuation orders from the **same fixed prefix**,
retain their final-site IDs V1,V2. Then

\[
P(V_1=V_2\mid\text{prefix})=\sum_v\pi_v^2.
\]

The single binary mark `1[V1=V2]` separates these mechanisms by exactly1/36;
it needs neither a full sitewise probability estimate nor a graph classifier.
It is invariant under any common relabeling of the fixed state's sites.
Independent continuation is essential; this is not a collision statistic for
unrelated prefixes or coupled/same-order continuations.

The temporal structure is informative too: conditional on T2, both mark
collision probabilities are1/4. At T3 they are5/18 versus1/4; at T4 they are
3/8 versus1/4. Thus this particular scalar mark does not distinguish the
graphs at the earliest possible birth, but does distinguish their later
branching geometry. The full joint(T,V) table supplies all these readouts.

The conclusion is deliberately specific: a complete unmarked clock can close
the event's cardinality distribution without identifying the underlying
branch/cycle trigger graph. A final-site mark resolves this exact degeneracy;
it is not claimed to identify arbitrary graphs uniquely.

## Artifact

`scripts/p334_isoclock_marked_counterexample.py` scores only these prescribed
two graphs,32 subsets each, with exact integer/fraction arithmetic in the
existing local research environment. Results are saved in
`results/p334-isoclock-marked-counterexample/exact_marked_clock.json`.
No graph-family scan, newN425 sample, MonteCarlo stream, networkDP, or remote
job was used. The P334 source connection is the factorization exposed in
`notes/p334-middle-bridge-physical-interventions.md` at0143632d.
