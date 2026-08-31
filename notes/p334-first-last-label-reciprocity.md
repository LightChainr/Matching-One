# Exact first/last label reciprocity for a monotone birth clock

The saved time-resolved **last** birth-site table already determines every
conditional survival curve after observing the **first** insertion label.
This gives the new first-step Doob readout on a solved prefix without another
network solve or a new ensemble of child continuations. The two complete
site-by-time readouts are equivalent information, not independent evidence.

## Fixed-cardinality identity

Fix a safe prefix with d remaining labels and a hereditary family of safe subsets.
Insert the labels once in a uniform random order. Let T be the first unsafe
insertion and V its label. If absorption is impossible, retain T=d+1 as a
cemetery time and do not invent a winner. Define

```
I_k = number of safe k-subsets;
c_v(k) = number of safe k-subsets containing v;
b_v(k-1) = number of safe (k-1)-subsets not containing v
           that become unsafe when v is added.
```

The safe (k-1)-subsets avoiding v are counted by `I_(k-1)-c_v(k-1)`.
After adding v, exactly c_v(k) of them are still safe, by heredity and the
bijection that removes v. Therefore, for 1<=k<=d,

```
c_v(0)=0,
c_v(k)=I_(k-1)-c_v(k-1)-b_v(k-1).
```

Each triggering set has `(k-1)! (d-k)!` compatible insertion orders, giving

```
b_v(k-1)=d*choose(d-1,k-1)*Pr(T=k,V=v).
```

Conditional on the first label U1=v, the first k labels form a uniform k-set
containing v, so

```
Pr(T>k | U1=v)=c_v(k)/choose(d-1,k-1),  1<=k<=d.
```

The probabilities at horizon zero are one. This recurrence recovers the
entire conditional survival vector using only the old safe counts and
time-resolved final-site table. It covers initial direct gates, inert labels,
and later collective births; there is no h>0 assumption or approximation.

Conversely, the first-label survival vectors recover c_v(k). Their uniform
average is the unconditioned survival, so they also recover I_k. The same
recurrence then recovers every b_v(k-1), hence the time-resolved winner table.
If there are nonabsorbing full sets, their probability is already retained
in I_d. No winner is assigned to them.

## First-step information and its response rank

Write `s_v(k)=Pr(T>k|U1=v)`, `s(k)=Pr(T>k)`. Then

```
G_(k,l) = (1/d) sum_v [s_v(k)-s(k)]*[s_v(l)-s(l)]
```

is the covariance Gram matrix of the next-label conditional survival vector.
For any terminal function expressed as a linear combination of survival
indicators, its first-label Doob variance is the corresponding quadratic
form in G. In particular,

```
m_v = E[T|U1=v] = 1 + sum_(k=1)^d s_v(k),
Var(E[T|U1]) = sum_(k,l=1)^d G_(k,l).
```

For a monotone absorption event, once v is inserted this also equals one plus
the remaining child-clock mean, with the absorbed-child clock set to zero.
This identity must not replace the full child law by its mean when computing
nonlinear canonical readouts.

Center the c and b tables across v at every horizon. Their relation is

```
centered_c(k) = -centered_c(k-1)-centered_b(k-1).
```

This is an invertible triangular temporal transform. The nonzero binomial
normalizations are invertible as well. Thus the full first-label response
rank equals the rank of the centered, time-resolved final-site table. An
observed rank difference is a genuine difference in spatial/temporal response
structure, but first-label and final-label descriptions do not double the
number of independent discoveries.

Only the **full time-resolved** winner table has this equivalence. The marginal
winner probabilities pi_v or their collision sum alone discard the necessary
timing information. Likewise one first-label mean per vertex does not recover
the whole final-site law. No minimal continuum-state dimension is inferred
from this finite response rank.

## Immediate scientific use

The prescribed five-label double-star and C4-plus-inert examples have equal
unmarked clocks and no initial direct gate. A binary direct/safe first-step
readout is therefore identically uninformative in both, while their full
first-label response profiles can have different rank. The ongoing exact
example computes that difference; it is not another unmarked-clock comparison.

The actual P334 prefixes already have solved safe coefficients and
time-resolved final-site marks in `1c06230b`. When their coefficient convention
is read correctly, the recurrence above provides a thin physical readout of
their first-label information. It calls for no new DP, path replay or random
sample. Both readouts inherit the same selected-prefix dependency group.
