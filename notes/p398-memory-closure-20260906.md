# #588 Phases B and C — the compact object exists, and #580 did not find it

`scripts/p398_memory_closure.py`, `results/p398-memory-closure/latest.json`,
`tests/test_p398_memory_closure.py` (7 tests). Widths 4–8, declared intervention,
`eta` in `{-1/4, 0, +1/4}`, everything frozen at `eta = 0`.

**This block revises #580's headline.** #580 reported that a frozen rank-6
observable-Krylov span transports but does not represent five held-out readouts, and
read that as a fact about the dictionary. Phase C says it was a fact about the *span*.
A finite-horizon balanced realization of **order 4**, built from baseline data only and
frozen at `eta = 0`, predicts all eight readouts' contrast responses at `eta = ±1/4`
to **5.9%** at width 8, where the rank-6 Krylov span sits at 41.8% on the held-out
block and never passes at any rank up to 12.

---

## Phase B — the growth profile under dictionary refinement

Five numbers per dictionary, kept apart because they are five questions. At the frozen
rank 6, widths 4 → 8:

```text
dictionary   outs  seed   r_linear                  r_transport      r_memory   balanced
                   span                             (Krylov)         (numeric)  order
D0            3      4    10, 26,  72, 218, >=150   4,  6, 8,  8,  8   4, 9,12,13,14   3,3,3,3,3
D1            6      7    10, 26,  76, 232, >=150   6,  8, 8, 12, 12   4, 8,14,16,17   3,3,3,3,3
D2            8      9    10, 42=, 76, 415, >=150   --, 12, --, --, --  4, 8,14,16,17   3,3,4,4,4
```

`=` means the Krylov space filled the state space; `>=` means the modular-rank budget
was reached at width 8 and the true value is larger; `--` means no rank in `(3,4,6,8,12)`
reached the declared threshold.

Refining D0 → D2 (three readouts to eight):

- **`r_linear` roughly doubles** and at width 7 reaches 415 of 429 states — an ordinary
  linear realization of the full dictionary is essentially the whole chain;
- **`r_transport` degrades from 8 to "never"**;
- **`r_memory` moves from 14 to 17**;
- **the balanced order moves from 3 to 4.**

Everything grows except the last one. That is the profile, and it is the answer to
"what is the state dimension": four different numbers that move differently under the
same refinement.

## Phase C1 — more state, and a defect of the Krylov ladder

Held-out balanced error, declared intervention, width 8:

```text
rank  +dof   eta=0    eta=+1/4   eta=-1/4   declared(eta=0)
   6     0   0.4176     0.4439     0.3887        0.1310
   7     1   0.4847     0.5054     0.4578        0.0563
   8     2   0.2350     0.2745     0.2179        0.0404
   9     3   0.2517     0.3439     0.2259        0.0373
  10     4   0.3136     0.3833     0.2985        0.0233
  12     6   0.1844     0.3030     0.1954        0.0066
```

**The ladder is not monotone.** Rank 7 is worse than rank 6 on the held-out block, at
widths 5, 6, 7 and 8; so is rank 9 against rank 8, and rank 10 against rank 9. The
declared column *is* monotone, which locates the problem: a rank cutting inside a
Krylov level adds half a level, and the Galerkin closure then uses that half-direction
to fit the declared readouts at the expense of everything else.

#580 noted that a prefix of a block-Krylov construction is still a declared span. It is
— but it is not a *nested improvement*, and "the smallest rank that transports" is
therefore a weaker notion than it sounds.

## Phase C2 — more memory, at the same rank 6

Fitted by Ho-Kalman on the `eta = 0` kernel and then frozen — poles and residues — and
applied unchanged at `eta = ±1/4`. Width 6, `eta = 0`, balanced error:

```text
variant                       declared   held-out   kernel fit
markov_only                    0.1074     0.2960        --
poles_1                        0.0948     0.2671      0.2779
poles_2                        0.0894     0.2475      0.1053
poles_3                        0.0786     0.2454      0.0556
poles_4                        0.0785     0.2449      0.0069
exact_kernel                   0.0785     0.2449        --
exact_kernel_with_forcing      0.0785     0.1884        --
```

Three findings.

**The memory closure saturates at three poles.** `poles_3` already reaches the
exact-kernel value to three digits, at every width. Four poles match it exactly. So the
low-order memory description is real: the kernel of this projection is a three-pole
object for every purpose this task has.

**It buys much less in the response metric than in the trajectory metric.** Phase A
measured the state-space trajectory `Phi^T exp(tG) f` and found the declared closure
error was 96% memory. Here, in the source-sampled contrast response, exact memory
removes 27% of the declared error and 17% of the held-out error. Both are correct: the
contrast response subtracts the common part that every model gets for free, so what
remains is the part memory is worst at. Reporting only the first number would have
overstated what a memory closure is worth.

**The unresolved initial condition is the larger half of the held-out gap**, again:
`0.2960 -> 0.2449` from memory, `-> 0.1884` when the exact forcing term is added too.
Consistent with Phase A's ~40/43 split.

## Phase C3 — the balanced adversary, and what it settles

Finite-horizon Gramians on the zero-sum contrast subspace. The stationary mode is
removed by construction and not by regularization: the sources enter as contrasts
against the declared reference, so they sum to zero, and a zero-sum functional
annihilates the constant function exactly. Nothing marginally stable is fed to a
Lyapunov solve, the horizon is the declared lag grid, each readout is divided by its own
baseline contrast signal *before* the Gramians are built (a balancing transform is
exactly the object a rescaled readout would silently rotate), and the transform is
frozen at `eta = 0` and transported.

Width 8, `error at eta = 0 / worst at eta = ±1/4`:

```text
dictionary   hankel sigma (normalized)          order 3        order 4        order 6
D0    1.000 0.583 0.136 0.040 0.006 0.003   0.0398/0.0577  0.0097/0.0135  0.0018/0.0054
D1    1.000 0.391 0.296 0.089 0.009 0.008   0.0835/0.0875  0.0197/0.0273  0.0054/0.0184
D2    1.000 0.418 0.325 0.120 0.021 0.014   0.1175/0.1340  0.0470/0.0591  0.0092/0.0241
```

The comparison that matters, on the full eight-readout dictionary at width 8:

```text
Krylov rank 6  + nothing          held-out 0.418
Krylov rank 6  + exact memory     (width 6) 0.245, and 0.188 with exact forcing
Krylov rank 12 + nothing          held-out 0.303
balanced order 4                  all eight readouts, 0.059
```

**A four-dimensional balanced state beats rank-6-plus-exact-memory by a factor of four
and rank-12 by a factor of five, on more readouts, with fewer coordinates.** And the
Hankel spectrum barely moves with width: D0 is `1, 0.45, 0.09, 0.014` at width 4 and
`1, 0.58, 0.14, 0.040` at width 8, across a factor of 102 in state count.

So the answer to the owner's three-way question is unambiguous:

```text
balanced state stays compact and transports
    -> #580's Krylov choice was not special; a genuine compact
       input/output module exists.
```

### What the compact object actually is

It is **not** a state of the system, and it is not a state of the dictionary either. It
is a state of the triple `(sources, readouts, horizon)`. The Hankel singular values are
computed from `mu_c exp(tG) f` for the declared contrasts and readouts over the declared
horizon; change any of the three and the spectrum changes. What is remarkable is only
that the number stays at 3–4 while `r_linear` climbs toward the full state space and
`r_positive` reaches 750.

### The caveat that has to travel with this

The balanced construction uses strictly more baseline information than Krylov: it is
built from the exact baseline Gramians on the declared sources, readouts and horizon,
i.e. it is fitted to the exact metric it is then scored on, at `eta = 0`. Krylov uses
only the readouts. Its advantage at `eta = 0` is therefore expected and is not the
finding.

The finding is that a transform frozen at `eta = 0` **transports**: on D2 at width 8 the
error goes from 0.047 at baseline to 0.059 at `eta = ±1/4`, an excess of `+0.012`,
which is the same order as the `+0.006` #580 reported for its own span on its own
readouts — but here on all eight readouts and at order 4 instead of rank 6.

Stability was checked rather than assumed: a Petrov–Galerkin reduction of a generator
need not be stable, and an unstable reduced model can score well on a finite grid before
diverging. `||exp(2T A_r)||` at twice the scoring horizon is below `10^-2` for every
retained order at every width — the reduced models decay, as the contrast subspace they
replace does.

Where the task does not expose enough coordinates, the construction says so instead of
inventing them: order 12 on D0 comes back `gramian_rank_deficient` at widths 4–8,
because three readouts and three sources over this horizon do not span twelve
input/output directions.

## What this changes

**#580's held-out failure was a property of the span, not of the dictionary.** The five
readouts `max_block, linked_pairs, halves_linked, covering_depth, boundary_span` are
representable and transportable — by an order-4 model. #580's conclusion that "the
low-dimensional state is a state of its readouts, not of the system" survives, and in
fact strengthens: it is a state of its readouts *and its sources and its horizon*. But
the accompanying suggestion that the held-out readouts needed either more state or a
memory description was answered by neither: they needed a **better-chosen** state of
smaller dimension.

**Phase A's candidate description is dominated.** "A small state plus three or four
poles" is a valid description and the poles are real — three of them reproduce the
kernel exactly. It is simply beaten by a four-dimensional Markov state with no memory
at all.

**The three-rank table gains a fourth column that behaves differently from all of them.**

```text
width   states   r_linear(D0)   r_positive   r_transport(Krylov)   balanced order (D2)
  4       14         10             10             4                     3
  5       42         26             26             6                     3
  6      132         72             76             8                     4
  7      429        218            232             8                     4
  8     1430      >=150            750             8                     4
```

## Boundary

One declared projection, source set, readout set and horizon on one exactly known
finite process. The balanced order is a property of that quadruple. None of this is a
percolation statement, and none of it transports to square-site Matching One without a
declared intervention with common observables before and after and a supplied map
between microscopic state spaces.
