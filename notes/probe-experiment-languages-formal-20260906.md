# Experiment languages as formal objects: separating what a language can see from what the process is

**C4 note (#599 Program A/G, north-star formalisation).**  The probe's working
hypothesis is that "complexity is indexed by a declared experiment language".
This note makes that claim precise enough to test: it defines *static* and
*test* languages, the class count `kappa`, the response rank `r`, and the
relation between language structure and these quantities, then attaches the
exact facts established in C2 (affine obstruction, root-count = depth) and new
P398 dictionary-signature numbers.  It does not claim a new universal theorem
of realisation theory; it claims a *language-first* notation in which the
repository's separate results become entries of one table.

---

## 1. Two kinds of experiment language

A **static (readout) language** `D` is a finite set of functions
`f: X -> R` on the state space.  Its *signature class count* is

```
kappa_D = #{ distinct value vectors (f(x))_{f in D} over x in X }.
```

This is the static Myhill–Nerode quotient of the dictionary.  A **test
language** `L` is a finite set of experiments `E_j`, each with a success
probability `p_j(x)` on `X`.  Its class count and its response rank are

```
kappa_L = #{ distinct rows (p_j(x))_j },
r_L     = rank of M = ( p_j(x_i) )_{i,j}.
```

Always `kappa_L >= r_L`?  No: `kappa_L = rank of the row set in the affine
sense only if distinct rows are linearly independent; in general kappa_L can
exceed r_L (C2) or equal it.  This is the whole point of the language-first
notation: class counting and rank counting answer different questions, and the
language object keeps them separate instead of collapsing them into "the state
dimension".

**Relation to earlier notions.**  The repository's #401 `Z0/Z1/Z2` hierarchy
is a test-language filtration; #588's D0/D1/D2 is a static-language
filtration; #582's quantile laws are a one-parameter family of *test
languages* indexed by scale word; #598's even/odd dictionary is a static
language with a declared symmetry.  All of them are special cases of "a
language" in the sense above, and none of them is the state space.

## 2. Structural facts (exact, from the probe)

**(F1) Monotonicity.**  `L subset L'` implies `kappa_L <= kappa_{L'}` and
`r_L <= r_{L'}`; a symmetry-invariant language can never separate states in
one `K`-orbit (Theorem 1/2, round 2): `kappa_{L^{K}} <= |X/K|`.

**(F2) Affine obstruction (#549, C2).**  For the parallel-gadget family with
class coordinate `a`, every single-root fork test has `p(a)` affine, so any
finite single-root language has `r_L <= 2` while `kappa_L = k+1`.  Language
class counting and linear independence are provably different coordinates.

**(F3) Root-count = depth = algebraic degree (C2).**  In the exact abstract
d-root model, degree(P)=d; a language of depth `d` (d independent roots) has
`r_L <= d+1` and attains it when the types are distinguishable.  Hence:
**the minimal distinguishing depth of the k+1 classes is k**, and clone width
at fixed depth does not raise the rank.

**(F4) Static language vs dynamics (P398, new numbers).**  Signature class
counts of the declared readout dictionaries on the *static* state space:

| w | n | C2-orbits | D0 | D1 | D2 |
|---|---|---|---|---|---|
| 4 | 14 | 10 | 8 | 9 | 10 |
| 5 | 42 | 26 | 12 | 19 | 25 |
| 6 | 132 | 76 | 18 | 40 | 49 |
| 7 | 429 | 232 | 24 | 81 | 109 |
| 8 | 1430 | 750 | 32 | 156 | 209 |

(data: `results/probe-dictionary-signature-classes/latest.json`).  Three
readings that would be invisible without the language-first notation:

1. **The full declared dictionary (D2) never resolves the C2 orbits at
   w >= 5** (209 < 750 at w=8) and at w=5–8 resolves far fewer classes than
   the state space.  The "8 readouts" are a genuinely coarse language.
2. **The D0 initial blocks (8,12,18,24,32) are exactly the Gate-1 "initial"
   counts** (#598), and D2 = 10 at w=4 equals the orbit count: at w=4 the
   declared dictionary happens to separate every C2 orbit (kappa_D2 = 10 =
   |X/R|), at larger w it does not.  So the "dictionary-compatible quotient"
   of #597/#598 is a statement about D0 admissibility under the dynamics, not
   about D2 resolving the orbits.
3. **The hierarchy of numbers 1430 (micro) > 750 (orbit/lumping) > 209
   (D2-static) > 10..32 (D0-static) > 3–4 (balanced task order)** is now each
   term attached to a different formal object: state space, dynamics
   symmetry/lumping, static dictionary, task compression.  None of them is
   "the state dimension"; all of them are well-defined functionals of a
   declared language.

## 3. Minimal distinguishing language

Define `d*(Sigma)` for a family `Sigma = {x_0,...,x_m}` of distinguishable
states as the smallest depth `d` such that some depth-`d` test language has
`kappa_L = m+1` and `r_L = m+1` (class *and* rank separation).  Facts F2–F3
give: for the #549 family `d* = m = k` (single-root languages never reach it;
depth-k languages do).  For P398 with the declared D2 dictionary, `d*` is not
known; the cheapest determination is to run the D2 signature filtration under
the dynamics (a Krylov/lumping computation at fixed width), which would also
decide whether the static gap `209 < 750` closes once time evolution is
allowed to mix states.

## 4. Boundaries and next uses

* F2/F3 are proven for the published #549 protocol algebra and the abstract
  gadget model; F4 is exact arithmetic on P398.  Nothing here is a
  percolation-threshold statement.
* "Minimal distinguishing depth" is defined for test languages of the
  fork/root type; the general definition for arbitrary compositional languages
  (#550 depth-2 composition) is left open and is the next formal target
  (it should subsume #401's interaction-order hierarchy).
* The language-first notation is a *tool*, not a claim that all languages are
  equivalent; the round-2 symmetry language (`K`-invariant tests = orbit
  quotient; character-marked tests = isotypic fibres) slots into F1 unchanged.

**Open ranked items from this note:** (1) P398 dynamical `kappa` for D2 and
the minimal distinguishing depth; (2) a general "compositional depth vs
algebraic degree" theorem for tree languages; (3) whether #550's witness is
captured by degree-2 root languages on the real N16 network.
