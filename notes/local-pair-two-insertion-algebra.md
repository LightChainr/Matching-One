# Two local pair insertions and an entry-regular Q1 completion

The four-line double closure of the C4-averaged pair kernel has a genuine
Q1 pole, despite the finite single-insertion Bell4 closures. Adding the
specified local singlet with coefficient one cancels the entrywise pole
and gives a consistent finite-network interaction with an exactly invisible
Q1 endpoint. Its Q derivative is nonzero, including already at **two**
insertions. This is a new declared microscopic completion, not an assertion
that a counterterm leaves the old pure-pair first-order score unchanged.

Base `7681eedd`; the single-insertion definitions are those of
[the local four-port pair kernel](closed-source-local-four-port-pair-kernel.md).
This note uses exact equality-pattern algebra only. It does not rerun
the N25 population, the previous U score, or any coupling point.

## 1. Kernels and the four-line double closure

Use ports a=N,b=E,c=S,d=W. Let `A2=iP_[Q-2,2]i^dagger` be the single-cut
kernel, and let `R` be the quarter-turn reshuffling of its four indices.
Define

```
K2=(A2+R A2)/2,
A0=(1-delta_ab)(1-delta_cd)/[Q(Q-1)],
K0=(A0+R A0)/2.                                             (1)
```

These K2,K0 are the C4 averages; R is a tensor reshuffling, not a
conjugation on the old two-port transfer space. For real integer Q>=4,
the complete four-line double closure is the Frobenius pairing

```
<K,L> = sum_(a,b,c,d) K(a,b,c,d)L(a,b,c,d).
```

Its rational continuation need not remain a positive norm outside
the integer representation range.

Write `d2=Q(Q-3)/2`. The unaveraged A2 is an orthogonal projector, so
`<A2,A2>=d2`. Its single-insertion Bell4 closures vanish except for
`NS|EW` and `NW|ES`, each equal to d2. In the partition-diagram expansion
of R A2, the coefficients of those two diagrams are respectively
`1/2` and `1/[(Q-1)(Q-2)]`. Therefore

```
<A2,R A2> = d2 [1/2+1/((Q-1)(Q-2))],

||K2||^2 = (d2+<A2,R A2>)/2
 = Q(Q-3)(3Q^2-9Q+8)/[8(Q-1)(Q-2)].                        (2)
```

In particular,

```
||K2||^2 at Q4 = 5/3,
Res_(Q=1) ||K2||^2 = 1/2.                                 (3)
```

This proves that the finite one-insertion continuation alone did not
define a pole-free arbitrary-strength local interaction. Two kernels
joined by four ordinary identity wires are already an obstruction.
An application may restrict its allowed spatial wirings, but cannot
infer unrestricted consistency from the old Bell4 table alone.

## 2. The singlet cancels the entrywise pole

The frozen completion is

```
Kreg=K2+K0 = C4-average of i(I-P1)i^dagger.                   (4)
```

Its unaveraged rational entries are explicitly

```
Areg(a,b;c,d) = (1/2)(1-delta_ab)(1-delta_cd) {
  delta_ac delta_bd + delta_ad delta_bc
  -(delta_ac+delta_ad+delta_bc+delta_bd)/(Q-2)
  +4/[Q(Q-2)] }.                                            (5)
```

There is no Q1 pole in any coefficient of (5). In contrast, writing
`U=(1-delta_ab)(1-delta_cd)` and using C4 averages,

```
K2 = -average(U)/(Q-1)+O(1),
K0 = +average(U)/(Q-1)+O(1).                                (6)
```

For a scalar coefficient c(Q) regular at Q1, entrywise regularity of
`K2+c(Q)K0` consequently requires **c(1)=1**. It does not require c(Q)=1.
The choice (4) is the declared constant-coefficient completion; it is
not a uniqueness theorem. Writing

```
c(Q)=1+(Q-1)a(Q)
```

leaves the additional regular local counterterm
`a(Q) average(U)/Q`. It has real consequences for the Q activation
specified below. No a(Q) is fitted or selected here.

The complete two-insertion cancellations can also be displayed directly:

```
||K0||^2 = (2Q^2-4Q+3)/[2Q(Q-1)],
<K2,K0> = (Q-3)/[4(Q-1)].                                  (7)
```

For the first identity, the overlap `<A0,R A0>` counts proper colourings
of a four-cycle, whose polynomial is `(Q-1)^4+(Q-1)`. For the second,
same-cut singlet/pair orthogonality gives zero, and the rotated singlet's
two-crossing diagram contributes `d2/[Q(Q-1)]` before averaging.
Combining (2) and (7) gives

```
||Kreg||^2
 = (Q-1)(3Q^3-12Q^2+20Q-24)/[8Q(Q-2)],

||Kreg||^2 at Q1 = 0,
partial_Q ||Kreg||^2 at Q1 = 13/8.                          (8)
```

The residues of the two squared pieces are each +1/2; the cross term
`2<K2,K0>` has residue -1. The cancellation is exact, not a subtraction
of numerically large evaluations. The nonzero 13/8 also demonstrates
that the Q derivative of a finite-strength interaction can contain
terms quadratic in its insertion amplitude.

## 3. Every one-insertion Bell4 closure and its Q activation

For an outside connectivity partition pi, let b=|pi|. Define the
relative closed coefficient

```
B_pi(Q) = Q^(-b) sum_(a,b,c,d) Kreg(a,b,c,d) D_pi(a,b,c,d).   (9)
```

The letters a,b inside the sum are colour indices; b=|pi| only denotes
the block count in the prefactor. Distinct outside blocks may carry
the same colour. They are not forced unequal by D_pi.

All 15 cases are as follows. At Q1 **every B_pi is zero**.

| Outside partition pi | Exact B_pi(Q) | partial_Q B_pi at Q1 |
|---|---|---:|
| N\|E\|S\|W | `(Q-1)/Q^3` | 1 |
| NE\|S\|W | `(Q-1)/(2Q^3)` | 1/2 |
| NS\|E\|W | `(Q-1)/Q^3` | 1 |
| NW\|E\|S | `(Q-1)/(2Q^3)` | 1/2 |
| ES\|N\|W | `(Q-1)/(2Q^3)` | 1/2 |
| EW\|N\|S | `(Q-1)/Q^3` | 1 |
| SW\|N\|E | `(Q-1)/(2Q^3)` | 1/2 |
| NE\|SW | `(Q-1)(Q-2)/(4Q^2)` | -1/4 |
| NS\|EW | `(Q-1)(Q-2)/(2Q^2)` | -1/2 |
| NW\|ES | `(Q-1)(Q-2)/(4Q^2)` | -1/4 |
| NES\|W | 0 | 0 |
| NEW\|S | 0 | 0 |
| NSW\|E | 0 | 0 |
| ESW\|N | 0 | 0 |
| NESW | 0 | 0 |

These follow by adding the singlet contractions to the old K2 table.
For example, the all-free K0 contraction is Q(Q-1); a single opposite
pair gives Q-1; and two opposite pairs give 1. Hence the opposite
2+2 closure of Kreg is `d2+1=(Q-1)(Q-2)/2`, rather than d2 alone.
Quarter-turn averaging halves the corresponding adjacent cases.

For an insertion at the named vacant site x, the Q derivative is the
explicit occupation statistic

```
1_(x vacant) partial_Q B_(pi_x(A))(1).                       (10)
```

The derivative of the background occupation weight does not contribute
to this *one-insertion first Q derivative*, because B_pi(1)=0 for each
configuration. The same holds for the derivative of Q^(-b) in (9).
At Q1, partial_Q and partial_logQ coincide. Any subsequent normalization,
thermal derivative or pooled-root transmission must still be performed
on the specified original observable.

The completion ambiguity is quantitatively visible: replacing c=1 by
`c=1+(Q-1)a(Q)` leaves the singleton and 2+1+1 activation rows unchanged,
but changes the opposite 2+2 slope to `-1/2+a(1)` and each adjacent
2+2 slope to `-1/4+a(1)/2`. All 3+1/all4 rows stay zero. This is why
freezing c identically one is a model definition, not an irrelevant
choice of finite counterterm.

## 4. Finite-network zero theorem, with its scope

Consider a finite closed network made from ordinary equality tensors,
constant vacant tensors and finitely many Kreg tensors. Allow scalar
coefficients analytic at Q1. Evaluate its colour sums at integer Q>=4
and continue the resulting rational partition-diagram expression.
Then every term with **at least one Kreg insertion vanishes at Q1**.

Proof: (5) expands in finitely many equality diagrams with coefficients
regular at Q1. Every closed equality diagram evaluates to Q raised to
its number of free colour components. Therefore continuation to Q1 is
identical to evaluating these regular tensors on the single-colour set.
On that set every delta is one and `(1-delta_ab)(1-delta_cd)=0`, so
each unaveraged Areg, and hence Kreg, is zero. A term containing one
or more such tensors is consequently zero. There is no exchange of
limits with an infinite network or an unresolved pole.

The statement remains true after a finite sum over original occupations
A, with their fixed q(A), E(A), and analytic scalar factors such as
`Q^(-r(A)/2)`. It requires that r refer to that stipulated original
occupation, not a new rank assigned independently to each virtual
diagram join. It also holds with the regular counterterms of Section 2,
which retain the same unequal-pair factor.

For example, a finite lattice with local factors

```
vacant: 1+epsilon Kreg,       occupied: v delta_all4
```

has, at Q1, exactly its epsilon=0 partition and original q/E numerators,
to all powers of epsilon in this finite product. Provided the original
pooled root is simple, its root and U at Q1 are likewise independent of
epsilon. This is an algebraic endpoint statement; no positive transfer
operator or continuum limit is needed.

The theorem does **not** cover singular external projectors, factors
such as `1/(Q-1)`, uncontinued literal permutation seams, or taking an
infinite-volume limit before Q1. Nor does it say the Q derivative of
each nonempty-insertion term is zero: (8) is an explicit two-insertion
counterexample. The rational one-colour argument applies to the value,
not to differentiation of the number of colours.

At fixed integer Q>=4, sufficiently small epsilon gives a genuine
nonnegative local tensor perturbation because the vacant tensor is
strictly positive and Kreg has finitely many finite entries. At
noninteger Q near one, the claim here is the declared analytic finite
occupation completion, not a literal noninteger set of colour states.

## 5. Scientific consequence and exact companion

The old pure K2 one-insertion Q1 statistic was finite and nonzero.
Adding K0 changes it: the new full insertion is identically zero at Q1
and has the activation table above. This supplies a consistent local
finite-network completion, but does not retrospectively preserve or
reinterpret the previous pure-K2 U score as the same observable.

The accompanying small script
[`local_pair_two_insertion_algebra.py`](../scripts/local_pair_two_insertion_algebra.py)
computes the rational identities and table by the 15 exact equality
patterns. A pattern with k distinct colour blocks has multiplicity
`Q(Q-1)...(Q-k+1)`; no finite-Q interpolation is used. One run with the
managed Python research environment returned (2), (7), (8) and every
row of the table exactly. It prints the results and creates no archive,
sampling output, or fitted model.

The new mechanism is a local interaction invisible at Q1 but visible
to a colour derivative, with explicitly nontrivial higher insertion
orders. It combines the pair and singlet sectors under a frozen finite
counterterm. No Jordan, CFT-field identity, scaling exponent, or new
original-U numerical response is inferred.
