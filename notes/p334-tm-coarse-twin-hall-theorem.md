# P334: the corrected reservoir has an exact second Hall compression

## Theorem

Let `Q` be the finite translation group of an arbitrary HNF quotient and let
`N=|Q|`.  Translation acts freely on every ordered D or F face because a
stabilizer fixes its first marked site.  After the existing simultaneous-
translation normalization, every corrected-reservoir source has a unique
description

```text
(replica, [D face], [F face], relative phase)
  in Z2^2 x (D/Q) x (F/Q) x Q.
```

Fix the first three coordinates and vary the relative phase `h`.  The
corrected target builder translates the flat face over every `g in Q`.
Replacing `h` by `h0` only reindexes that loop because

```text
{g h : g in Q} = Q = {g h0 : g in Q}.
```

The occupied-to-vacant carrier exchange, transverse mark release, motif test
and target normalization are translation equivariant.  Therefore all `N`
relative-phase sources have exactly the same normalized target neighbourhood.

Call the phase-free source classes `C`, and write `R(c)` for this common
neighbourhood.  The normalized orbit graph has a saturating matching if and
only if the capacitated graph

```text
source -> c       capacity N
c -> t            infinite, t in R(c)
t -> sink         capacity 1
```

has flow `N|C|`.  By max-flow/min-cut, this is equivalent to

```text
|union_(c in A) R(c)| >= N |A|       for every A subset C.
```

The exact remaining Hall obstruction is consequently

```text
Delta_res = max_(A subset C) [N|A|-|union R(A)|].
```

An integral coarse flow expands to distinct relative-phase targets and then,
by the previously proved translation-orbit theorem, lifts to a collision-free
matching of the raw configuration graph.  This is a general theorem for every
finite HNF; no cyclic Smith assumption is used.

## Compression and frozen gates

The source reduction is now by `N^2` relative to the raw graph:

| gate | raw sources | translation-orbit sources | demand-N classes |
|---|---:|---:|---:|
| each N6 row | 1,152 | 192 | 32 |
| N8 Smith-(2,4), k=4 | 46,080 | 5,760 | 720 |

All `N6` twin neighbourhoods were also compared directly: each of the 32
classes contains six identical rows.  The N8 equality follows from the exact
group-reindexing proof and is protected by the class-size gate.

## The target channel changes

Because `MM` and `YN` targets are disjoint, the same network can be scored on
each channel separately.

| gate | combined | MM only | YN only | pure saturating channel |
|---|---:|---:|---:|---|
| N6 row 0 | 192/192 | 152/192 | 192/192 | YN |
| N6 row 1 | 192/192 | 150/192 | 192/192 | YN |
| N6 row 2 | 192/192 | 152/192 | 192/192 | YN |
| N6 row 3 | 192/192 | 150/192 | 192/192 | YN |
| N8 Smith-(2,4) | 5,760/5,760 | 5,760/5,760 | 4,800/5,760 | MM |

The failed pure-channel flows return exact minimum cuts:

- the two N6 MM cut types have `(classes,demand,neighbours,deficiency)` equal
  to `(16,96,56,40)` and `(24,144,102,42)`;
- the N8 YN cut is `(576,4608,3648,960)`, evenly split as 144 coarse classes
  from each source replica.

Thus neither one universal YN injection nor one universal MM injection can be
the arbitrary-HNF proof.  The successful channel already changes between the
minimal and first nontrivial Smith gates.

## Constructive next form

A general injection may still be layered.  Assign each coarse class an
integer `b(c)` and route

```text
b(c)       units through MM,
N-b(c)     units through YN.
```

Pure YN and pure MM are the endpoints `b=0` and `b=N`.  A row on which the
combined graph saturates but neither endpoint does is the first genuinely
mixed-channel gate.  Its flow supplies the exact candidate allocation that a
topological proof must explain.

Future scans should therefore report:

1. exact size-`N` twin classes;
2. combined, MM-only and YN-only coarse flow;
3. pure versus genuinely mixed classification;
4. on failure, the residual minimum-cut class ranges, demand, target union,
   deficiency and deterministic descriptor hashes.

## Boundary

The second compression and minimum-cut formula are proved.  Corrected-
reservoir saturation for every HNF is not.  Digital Alexander complement has
not yet supplied a uniform bound on overlaps between `R(c)` for distinct D/F
orbit pairs, nor a layer/carrier rule selecting `b(c)`.

The implementation is `scripts/p334_tm_coarse_reservoir_hall.py`; it reads the
unchanged corrected-reservoir observable from `4bb7517`.
