# Same-parent conjugate norm-5 coalescence tomography

## Research boundary

This note responds to a specific weakness in the post-P57 mechanism story. The
observed norm-5 conjugation-odd thermal-jet diagnostic compares two different
parent lineages, so multiplier conjugation is entangled with parent size and
geometry. It is useful as a discovery clue, but it is not by itself a clean
measurement of a conjugate Gaussian quadrature.

The exact Gaussian arithmetic supplies a direct control with no radial model:
complete the missing conjugate norm-5 branch from each P57 parent. For the two
parent pairs already used by P57, both missing branches coalesce to one new
noncyclic child at the same N.

## Exact coalescence

For N=65,

```text
(8+i)(2-i)   = 17-6i
(7+4i)(2-i) = 18+i

(8+i)(2+i)   = 15+10i
(7+4i)(2+i) = 10+15i.
```

The last two Gaussian numbers are D4-equivalent, so the conjugate branch has a
single child class `(15,10)`, norm 325 and Smith group `(5,65)`.

For N=85,

```text
(9+2i)(2+i) = 16+13i
(7+6i)(2+i) = 8+19i

(9+2i)(2-i) = 20-5i
(7+6i)(2-i) = 20+5i.
```

Again the conjugate branch coalesces, now to `(20,5)`, norm 425 and Smith group
`(5,85)`.

This is useful for two independent reasons. First, it removes the parent-lineage
confounder from conjugate-multiplier spectroscopy. Second, the coalesced node is
noncyclic while the observed P57 child nodes are primitive/cyclic, making the
same block a matched arithmetic-class control.

## Same-N H4 law

At fixed N and fixed p, a reflection-even scalar response with one H4 anisotropy
has the form

```text
M(theta) = S_N + A_N cos(4 theta).
```

The scalar `S_N` and amplitude `A_N` are unrestricted. Three orientations
therefore give one exact linear relation.

For N=325 the exact cos4 values imply

```text
M_(15,10) = 11/5 M_(17,6) - 6/5 M_(18,1),
```

or

```text
5 M_(15,10) - 11 M_(17,6) + 6 M_(18,1) = 0.
```

For N=425,

```text
M_(20,5) = -13/20 M_(16,13) + 33/20 M_(19,8),
```

or

```text
20 M_(20,5) + 13 M_(16,13) - 33 M_(19,8) = 0.
```

These relations are stronger than another cross-size H4 fit. They contain no
`13/8`, no threshold value, no root conversion, and no nonuniversal amplitude.

The same identities may be read as true same-parent conjugate quadrature
relations. The parent spin-4 sine ratios are exactly

```text
sin(4 theta_(8,1)) / sin(4 theta_(7,4)) = 6/11,
sin(4 theta_(9,2)) / sin(4 theta_(7,6)) = 33/13.
```

The corresponding child conjugate differences have the same ratios under a
single H4 response.

## Angular adversaries

A three-angle same-N block can test H4 against other single harmonics without
using a radial exponent. The frozen interpolation weights are in
`predictions/norm5_conjugate_coalescence_20260829.yaml` and are independently
recomputed by `scripts/verify_norm5_conjugate_coalescence.py`.

The H8/H12 weights are deliberately very different from H4, especially at
N=325. The control is therefore informative even though P57 already disfavored
the old H8/H12 transfer aliases: it asks whether the same angular law survives
when the third point changes the quotient arithmetic class.

## Common-field three-geometry block

The clean implementation uses new disjoint data rather than treating historical
P57 A/B values as a fixed source.

At each N run two general-period pairs with the same seed and counter interval,
keeping C in the same slot:

```text
N325: (15,10) vs (17,6)
      (15,10) vs (18,1)

N425: (20,5)  vs (16,13)
      (20,5)  vs (19,8).
```

The rank priority permutation is a function of N and the counter stream, not of
the second geometry. Therefore the repeated C histogram should be byte-identical.
That statement is an implementation gate, not an assumption: verify it on a
tiny run. If it fails, retain both C streams and propagate their measured
covariance explicitly.

With the identity gate passed, the two pair runs form one synchronized
three-geometry block. Compute each frozen linear residual inside delete-one
batches. This can yield much smaller variance than combining a new C mean with
the old P57 A/B block.

## Why this should precede a strong interpretation of norm 4

The planned norm-4 targets N260/N340 are also noncyclic. If their q2/Jordan
closure fails, the first interpretation question will be whether the failure is
a radial/Jordan mechanism or a quotient-arithmetic effect. The coalescence
triangle asks that question at fixed N, where radial scaling is absent.

A pass therefore strengthens the use of norm 4 as an even-generator test. A
failure specifically at the coalesced C node would instead make quotient
arithmetic or nonlocal topology the leading explanation before another LCFT
operator is introduced.

## Compute gate

Do not automatically displace the frozen balanced norm-4 production. First run
a 1M/10M variance pilot, using target means only for variance estimation. If the
H4 residual and fixed-adversary separations are underpowered per CPU, keep the
construction as a later arithmetic control and proceed with norm 4.

The point is not to create another campaign. It is to cheaply test the exact
confounder that currently limits interpretation of both the P57 conjugation
residual and a future noncyclic norm-4 result.
