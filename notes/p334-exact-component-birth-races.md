# P334: which parallel channel wins the complete birth race?

## Outcome

The twelve frozen prefixes now have exact **winning-component versus birth
time** distributions. Their small component safety polynomials are saved
permanently, so future channel readouts need no further rank oracle or DP.
This increment exports coefficients that the preceding whole-clock run had
computed but not retained; it is not described as a zero-DP archival lookup.
The subsequent boundary algebra and direct/collective split require no DP.

For the representative clock-crossing pair, counter83 has an extra immediate
gate but a weaker collective completion route than counter1006. The complete
race allocation makes that distinction concrete:

| Prefix / channel | Core sites | Original H2 | Wins overall | Share among T>40 births | E[T given channel wins] |
| --- | ---: | ---: | ---: | ---: | ---: |
| 83 / 0 | 71 | 2 | 45.9074% | 59.8469% | 20.0549 |
| 83 / 1 | 5 | 2 | 25.7058% | 18.1718% | 15.8927 |
| 83 / 2 | 11 | 2 | 28.3869% | 21.9814% | 16.5054 |
| 1006 / 0 | 107 | 3 | 79.0869% | 88.6684% | 17.7346 |
| 1006 / 1 | 2 | 2 | 20.9131% | 11.3316% | 13.2890 |

Counters here omit the common43042500000 prefix. Channel indices are ordered
by the minimum original vertex label of the corresponding component in the
saved map. The conditional times are for the **global race given that channel
wins**, not isolated-channel waiting times.

In counter83 the larger channel is increasingly important in the long tail;
the two small channels dominate a larger part of early completion. In1006,
the107-site channel dominates both the whole event and the long tail and
finishes earlier even conditional on winning. This gives an explicit channel
interpretation of the previously observed crossing at steps10/11.

## Boundary polynomial and exact first-hit allocation

Let factor i contain n_i random sites, with safe-subset polynomial F_i(z).
The coefficient of z^k in

\[
B_i(z)=n_iF_i(z)-(1+z)F_i'(z)
\]

is

\[
(n_i-k)f_{i,k}-(k+1)f_{i,k+1}.
\]

It counts pairs `(safe k-set U, unselected site v)` for which inserting v
connects that component's two ports. Therefore, with r free sites and d=173,

\[
W_i(z)=B_i(z)\prod_{j\ne i}F_j(z)(1+z)^r,
\qquad
\Pr(T=k+1,\ I_{\rm winner}=i)=
\frac{[z^k]W_i(z)}{d\binom{d-1}{k}}.
\]

The factors have disjoint variable sets, so the final insertion belongs to
exactly one competing channel. The saved integer boundaries obey

\[
\sum_i W_i(z)=dF(z)-(1+z)F'(z),
\]

and the winning probabilities sum to exactly one. This reconstructs the
whole physical clock at every step; it is not an independence assumption
about fixed-k survival probabilities.

All local F_i, local B_i, full W_i coefficient arrays, original site lists,
port addresses, rational winning probabilities, conditional times and tail
shares are under `results/p334-component-birth-race/`.

## Direct original gates versus collective completion

Let h_i be the component's **original checkpoint** singleton-trigger count.
The boundary splits exactly as

\[
B_i=h_iF_i+\left[(n_i-h_i)F_i-(1+z)F_i'\right].
\]

The first term counts final insertions at an originally unsafe singleton;
the second counts collective completion at all other sites. This is not a
claim about a site's microscopic arm type, nor a redefinition using the
evolving instantaneous H2 after other insertions.

For the global system, a component's original-direct contribution is simply
`h_i F(z)`. Hence the already exported coefficients suffice for this further
split. The result for the crossing pair is:

| Prefix / channel | Original-direct wins | Collective wins | Collective share among all T>40 births |
| --- | ---: | ---: | ---: |
| 83 / 0,71 sites | 22.5959% | 23.3114% | 46.6933% |
| 83 / 1,5 sites | 22.5959% | 3.1098% | 5.0181% |
| 83 / 2,11 sites | 22.5959% | 5.7909% | 8.8278% |
| **83 total** | **67.7878%** | **32.2122%** | **60.5392%** |
| 1006 / 0,107 sites | 31.3697% | 47.7172% | 71.6711% |
| 1006 / 1,2 sites | 20.9131% | 0 | 0 |
| **1006 total** | **52.2828%** | **47.7172%** | **71.6711%** |

Thus all collective completion in1006 resides in its107-site channel. In83,
the71-site channel is the main but not exclusive collective route. Its
overall collective share is smaller, despite three competing channels and
more direct gates. This agrees with the preceding exact survival ratio:
at step40,83's direct-avoidance disadvantage0.761905 is overcome by its
collective-survival advantage2.301917, leaving total survival ratio1.753842.

## Beyond the displayed pair

The readout covers all twelve selected prefixes and their47 two-port factors.
For example, counter904 has14 competing channels: its largest winning share
is only17.9230% (23-site factor), rather than a single dominant bottleneck.
Counter48's67-site factor grows from40.6026% of all births to46.5512% of the
T>40 tail; counter1013's49-site factor grows from27.4265% to36.2576%.
The four one-factor prefixes have winning probability one as required.

These are exact allocations conditional on the same frozen real prefixes,
not a new population sample or a new model fit. No sitewise forcing sweep,
server, new MC, or higher-order trigger enumeration is used. Small factors
were exported once and reused for the additional algebraic split.

```sh
/Users/lc/python-envs/research-py311/bin/python scripts/p334_component_birth_race.py
```

Parent `bd95f2a048d5780568b689bd42e0a684daf74315`; original selection remains
manifestb9cbe13e. The allocation labels which route actually completes a
sampled-order birth; it is not, by itself, a causal decomposition of the
between-prefix mean difference or a continuum operator identification.
