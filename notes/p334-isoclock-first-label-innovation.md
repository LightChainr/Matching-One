# P334: equal full clocks can have different first-label information ranks

The prescribed double-star and C4-plus-isolate mechanisms have identical
complete unmarked birth clocks and identical uniform-blockade mean response.
Nevertheless, their next-label conditional clock responses have **different
exact covariance ranks, 2 versus 1**. The next label explains 17.6587% versus
12.6984% of total clock variance, although the direct-versus-safe floor is
zero in both examples.

This is a new readout of the two fixed graphs from `250c5899`, not a graph
search or a new real N425 prefix. The double star was motivated by a real
residual factor; the alternate C4-plus-isolate remains a constructed
pair-trigger mechanism, not a claimed embedding on the same square torus.

## Fixed graphs and conditioning time

Keep the original labels and edges:

- Double star: 01,02,03,14.
- C4 plus isolate: 01,12,23,30, with label4 inert.

T is the first insertion completing an edge. The new mark is V1, the **first
inserted label**, uniform on all five labels. It is not the label inserted
at the random stopping time T. We enumerated all 5!=120 orders per graph,
exactly 24 orders for each specified first label. This is finite exact
enumeration, not Monte Carlo or a repeated noise-semigroup computation.

For a general current prefix and first choice v, define the child wait
T_child(v) after that insertion. A directly absorbing choice enters a
cemetery child with wait zero. Then

\[
m_v=E[T\mid V_1=v]=1+E[T_{child}(v)],\qquad
\operatorname{Var}(E[T\mid V_1])=
\frac1d\sum_v(m_v-E T)^2.
\]

Both prescribed graphs start with no direct singleton absorber, h=0. The
binary direct/safe mark is constant and its bound B=0, but the identity of
the safe label can still expose nonzero clock information.

## Exact conditional means

Both unconditional clocks have P(T=2,3,4)=(2/5,2/5,1/5), E T=14/5 and
Var(T)=14/25. The conditional means are:

| Graph / first label | Role | E[T given first label] | E[remaining child wait] |
|---|---|---:|---:|
| Double star 0 | degree3 hub | 9/4 | 5/4 |
| Double star 1 | degree2 hub | 8/3 | 5/3 |
| Double star 2 | leaf at hub0 | 37/12 | 25/12 |
| Double star 3 | leaf at hub0 | 37/12 | 25/12 |
| Double star 4 | leaf at hub1 | 35/12 | 23/12 |
| C4 labels0,1,2,3 | cycle vertex | 8/3 | 5/3 |
| C4 label4 | inert isolate | 10/3 | 7/3 |

Consequently,

\[
\begin{array}{c|cc}
 &\operatorname{Var}(E[T\mid V_1])&
 \operatorname{Var}(E[T\mid V_1])/\operatorname{Var}(T)\\
\text{double star}&89/900&89/504=17.658730\%\\
\text{C4+isolate}&16/225&8/63=12.698413\%
\end{array}
\]

Their explained-fraction difference is exactly 25/504, or 4.960317
percentage points. Mean remaining clock variance after the complete first
label is 83/180 versus 22/45. These are exact population moments over the
finite uniform orders, not estimates with sampling standard errors.

## The full survival-vector Gram has different rank

Let s_v(k)=P(T>k|V1=v), k=0,...,5. Its mean vector is the common
S=(1,1,3/5,1/5,0,0). Define

\[
\Gamma_{k\ell}=\frac15\sum_v[s_v(k)-S(k)][s_v(\ell)-S(\ell)].
\]

Only cuts k=2,3 vary. The exact active covariance blocks are

\[
\Gamma_G=\begin{pmatrix}1/25&13/600\\13/600&7/450\end{pmatrix},
\qquad
\Gamma_H=\begin{pmatrix}1/25&1/75\\1/75&1/225\end{pmatrix}.
\]

Double-star determinant is 11/72000>0, with two positive eigenvalues
(50+-sqrt(2005))/1800, approximately .05265401464 and .002901540914.
C4 has determinant zero and one nonzero eigenvalue 2/45. The full six-cut
Grams have respectively four and five additional zero eigenvalues.

The earliest nonconstant cut has the **same variance 1/25** in both
mechanisms. The difference is in how first-label effects at cuts2 and3
couple: every C4 conditional curve moves in one common direction, whereas
the double-star has two linearly independent temporal response directions.
This is a rank of conditional observable responses, not a count of
continuum fields or an asserted minimal physical state dimension.

The scalar innovation is recovered by summing all Gram entries:
ones^T Gamma ones=Var(E[T|V1]). More generally this Gram gives the next-label
information for every fixed linear functional of the clock survival curve.

At the combinatorial level, if I_k(v) counts safe k-subsets containing v,

\[
s_v(k)=I_k(v)/\binom{4}{k-1},\quad k\ge1,\qquad
\sum_v I_k(v)=k I_k.
\]

The unmarked safe polynomial determines the average of these local
profiles, not their covariance across first labels. The present pair shows
that this missing second-order structure can change rank, not just amplitude.

## First label is not the final winner

The earlier collision statistic concerns V_final=the insertion at time T:
two independent continuations share that final label with probability
sum_v P(V_final=v)^2. Its already published values 5/18 versus1/4 were not
recomputed here. That statistic is not a next-step Doob variance.

The isolate makes the distinction sharp: C4 label4 **never** wins the
terminal birth, but is the first label with probability1/5 and has
conditional mean10/3 instead of8/3. The entire rank-one first-label
survival variation in C4 is the distinction between that inert first label
and the four cycle labels. A site irrelevant as a final trigger can carry
substantial early timing information.

The calculation therefore gives an exact example of information hidden
behind a zero direct-event floor and a fully known clock. It does not yet
identify a unique graph from this Gram, add independent evidence to the
same two constructed mechanisms, or imply intrinsic temporal memory in
the sampled torus process.

Executable: `scripts/p334_isoclock_first_label_innovation.py`.
Complete per-label counts, conditional curves, rational Grams and spectra:
`results/p334-isoclock-first-label-innovation/exact_first_label.json`.

Subsequent interpretation (`31c17d48`): the **full time-resolved** final-site
law, together with the unmarked clock, determines these first-label profiles
by an invertible triangular temporal transform. Rank2 versus1 is therefore
also the response-rank difference of the centered full final-site tables.
It is a new exact interpretation of the same marked laws, not independent
evidence. Marginal final-winner probabilities alone do not determine it.
