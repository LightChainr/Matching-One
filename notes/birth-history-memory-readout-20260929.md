# A concrete process question: does rank one forget its birth history?

2026-09-29. #778 continuation in Draft #838. This is a new analysis of
existing joint-birth archives, not new sampling or independent validation.

**Outcome.** A single preselected three-time comparison finds a small positive
history association in the large-size insertion data. It does not settle the
continuum Markov question. The existing exact L3 archives provide a stronger,
finite conclusion: the square rank process is Markov under both insertion
and iid-label clocks, while the triangular rank process is not. Closed-form
two-time label kernels certify the distinction without repeating enumeration.

This separates three objects which should not be conflated: an exact finite
counterexample, an exploratory large-size signal, and an unidentified limit.

## 1. One null with an actual consequence

For the two births T1<=T2, define

```text
R(t) = 1{T1<=t} + 1{T2<=t},
H(s,t) = P(T1<=s, T2>t),          s<=t.
```

If current rank alone is a time-inhomogeneous Markov state, then for every
a<b<c with positive rank-one risk,

```text
H(a,c)*H(b,b) = H(a,b)*H(b,c).                           (1)
```

The [kernel characterization](birth-kernel-markov-characterization-20260929.md)
proves the converse when (1) holds for all triples. Direct double jumps are
allowed; they never enter the rank-one risk set. This concerns the ordinary
Markov property in the natural rank history, not the full configuration,
an arbitrary hidden-state process or a strong Markov assertion.

Here is a directly interpretable violation. Condition on R(b)=1 and compare
future survival to c between paths whose first birth occurred by a and paths
whose first birth occurred later:

```text
Delta = P(T2>c | T1<=a, T2>b)
      - P(T2>c | a<T1<=b<T2).                           (2)
```

If both groups have positive mass, Delta=0 is equivalent to (1). A nonzero
value rules out every current-rank-only Markov model at this time triple,
including models with arbitrarily time-dependent rates. It does not rule
out models that retain connectivity, direction or another hidden state.
Conditioning on rank one is part of this probabilistic null, not a causal
intervention on birth time.

Static rank curves and adjacent fluxes admit a Markov completion. That
completion can reproduce those inputs while violating the actual longer
history. Equation (1) tests a constraint that fitting the static curve
cannot certify.

## 2. Existing large-size block: modest evidence, no cutoff hunt

The calculation uses the 112 production batches in
[`cloud-data-5k/run.json`](../analysis/birth-gap-20260929/cloud-data-5k/run.json):
14 independent batches x 5000 filtrations for each of two lattices and
L=64,128,256,512. No benchmark rows, local-pilot pooling or new samples.

Before this scoring, choose exactly one triple in each cell: the .25, .50,
.75 **lower quantiles of the equal mixture of J1,J2**. Time here is occupied
insertion count, not iid-label time. L512 is the primary available-size
readout; smaller sizes are context. These choices are fixed in the script,
but the block was already examined for gap statistics: this is retrospective,
not a prospective model-elimination certificate.

| L512 insertion process | Early-group survival | Late-group survival | Delta, percentage points |
|---|---:|---:|---:|
| Square | 0.44914 | 0.43491 | +1.422 +/- 0.511 |
| Triangular | 0.43248 | 0.42546 | +0.702 +/- 0.515 |

Errors are 14-batch delete-one standard errors, **not confidence intervals**.
Every deletion recomputes the three empirical cutoffs and both risk groups.
The JSON retains raw four-cell counts, all deleted estimates, cutoffs,
source batch identities and the full six-metric covariance. A standard error
does not account for the wider adaptive research programme.

The four size-specific differences are all positive in each lattice, but
we do not fit a limiting constant or exponent, combine eight signs into a
new significance claim, or treat related gap readouts as extra evidence.
At L512 the square difference is suggestive and the triangular one less
resolved. No continuum non-Markov theorem follows from this table.

The source geometries are primitive-coordinate square and triangular tori
with different moduli. Differences between the rows are not a same-modulus
universality test. The insertion/label clock comparison proved earlier is
asymptotic, not a finite-N license to replace this observable by a label-time
one.

## 3. Exact finite result from the saved complete tables

The existing L3 archives each contain all 9! permutations. Their paired
birth histogram determines the entire rank-only insertion path. Evaluating
all 120 triples 0<=a<b<c<=9 with exact integer counts gives:

- Square: all Markov kernel identities hold. By the characterization,
  its **L3 insertion rank process is Markov**.
- Triangular: two identities fail. At (a,b,c)=(4,5,6), the early/late by
  survive/exit table is `[[23328,101088],[15552,93312]]`. Survival probabilities
  are 3/16 and 1/7, so Delta=5/112 and conditional covariance=1/90.
  The other violating triple, (3,4,5), has Delta=2/25 and covariance=1/90.

This is an exact finite positive control alongside an exact finite negative
control, **conditional on the stored complete engine tables**. No 9!
enumeration was repeated. The all-triples check is a finite identity
calculation, not a search for a significant production cutoff. It also shows
why nontrivial microscopic connectivity does not by itself prove a rank
projection is non-Markov at every size.

## 4. Continuous-label kernels, computed exactly rather than substituting Np

The uniform permutation is independent of the ordered iid uniform labels.
For p<=q, the numbers of sites with label at most p, between p and q, and
above q have Multinomial(N; p,q-p,1-q) law. Thus the exact clock bridge is

```text
H_label(p,q) = sum_(0<=k<=ell<=N)
  [N!/(k!(ell-k)!(N-ell)!)]
  p^k (q-p)^(ell-k) (1-q)^(N-ell) H_J(k,ell).          (3)
```

This is not H_J evaluated at expected counts. Applying (3) to the two
complete L3 tables and expanding over the rationals gives:

```text
H_square(p,q) = 6 p^3 (1-q^2)^3,

H_triangular(p,q) = 9 p^3 (1-q)^3 [A(p)+B(q)],
A(x) = 2x^3-6x^2+3x,
B(x) = -2x^3+3x+1,                  0<=p<=q<=1.        (4)
```

The square kernel separates as a function of p times a function of q.
It therefore satisfies (1) identically, not just at the rational check
points: **the square L3 rank process is Markov in iid-label time too**.
On its rank-one interval its survival transition is

```text
P(R(t)=1 | R(s)=1) = [(1-t^2)/(1-s^2)]^3,
q_12(t) = 6t/(1-t^2),                         0<s<t<1. (5)
```

This is the rank-one exit rate, not the 0->1/0->2 rate or a large-L claim.

For the triangular kernel the determinant in (1) is exactly

```text
81 a^3 b^3 (1-b)^3 (1-c)^3
   [A(b)-A(a)] [B(c)-B(b)].                            (6)
```

It is not identically zero. For example at label times (1/3,1/2,2/3),
Delta=-1156/716639. This certifies finite non-Markovianity under the actual
iid-label clock. The sign need not be positive at every triple; the
determinant's explicit factors show how it changes. Do not attribute the
sign difference from the insertion example solely to the clock: the tested
time triples are different too.

The coefficient matrices of the two polynomials have ranks one and two.
These are kernel separation ranks, **not** assertions about a minimum
positive latent-state realization or the topology's full predictive memory.
Neither L3 identity is extrapolated to large size or to a continuum limit.

## 5. What this changes next

The useful question is no longer whether one can invent a history-sensitive
readout. We now have an exact null identity, finite controls of both signs
of the Markov verdict, and a small large-size contrast on the actual archive.

The next discriminating theoretical object is the continuum rank-one kernel
H(s,t), not more equivalent gap descriptors. A rank-only mechanism must
predict its multiplicativity; a proposed geometric mechanism should predict
the defect's dependence on entry history or a specified homology mark.
The current birth archive has no such direction marks, so it cannot test
conditional sufficiency of a direction state by retrospective invention.

Adding first-birth time as a state variable is not by itself a mechanism
discovery. When the current rank is one, the whole observed rank past is
already encoded by that one birth time. A time-dependent rank-plus-age
description is therefore a complete encoding of this observation history;
it does not explain how microscopic geometry produces the conditional
completion law, nor imply a time-homogeneous semi-Markov law.

No larger-size production or finer-cutoff search is launched here. If a new
block is used later to assess the small large-size contrast, its clock,
fixed time triple and microscopic mechanism prediction should be stated
before scoring. That is attention guidance, not a research lock.

## Reproduce and evidence boundaries

From the repository root:

```sh
python3 analysis/birth-markov-kernel-20260929/analyze.py
python3 analysis/birth-markov-kernel-20260929/exact_clock.py
```

The first script uses the standard library; the second uses SymPy exact
polynomial arithmetic. Executed with managed research Python 3.11 on arm64;
runtime/version details are in the outputs. The analysis skill's validation
pass checks source identities, full batch counts, weighted risk denominators,
exact determinant algebra and complete-batch covariance; it does not certify
a continuum conclusion or turn this reused block into independent evidence.

- [Readable archive result](../analysis/birth-markov-kernel-20260929/RESULT.md)
- [Raw counts, all delete-one estimates and covariance](../analysis/birth-markov-kernel-20260929/result.json)
- [Exact label-kernel polynomial certificate](../analysis/birth-markov-kernel-20260929/exact-clock.json)
- [General kernel characterization and Markov completion](birth-kernel-markov-characterization-20260929.md)

These are author analysis/proof outputs, not independent certification of
the older engine or a novelty claim. Original archives and previous results
are unchanged. No cloud job, original-U contract change or PR merge.
