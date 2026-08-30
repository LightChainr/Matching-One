# F2/F3 flat-twist representation and the charged F3 doublet

Status: exact finite representation theory plus a reinterpretation of the
exploratory N65 contrast from `a7cb19a`.  There is no new Monte Carlo block and
no continuum-field assignment.

## Convention and exact action

Write primitive homology lines as columns and order them by

```text
F2: [x=(1,0), d=(1,1), y=(0,1)]
F3: [x=(1,0), d+=(1,1), d-=(1,-1), y=(0,1)].
```

The realizable generators and a reflection are

```text
S = [0 -1; 1 0],   T = [1 1; 0 1],   m = [1 0; 0 -1].
```

They act on homology lines by `ell -> g ell`.  A twist covector transforms
contragrediently, so its kernel line follows the same projective permutation.
The machine certificate verifies `S^2=(ST)^3=1`, `T^q=1` on
`P^1(F_q)`, all rational projector identities, and the group orders below.

## D4 does not finish the classification

For F2, reduction kills the reflection and the projective D4 image has order
two.  Its permutation representation is

```text
2 A1 + B1

A1 scalar:       (1, 1, 1)
A1 H4 alias:     (1,-2, 1)
B1 axis odd:     (1, 0,-1).
```

For F3 the projective D4 image is the Klein four group and

```text
2 A1 + B1 + B2

A1 scalar:       (1, 1, 1, 1)
A1 H4 alias:     (1,-1,-1, 1)
B1 axis odd:     (1, 0, 0,-1)
B2 diagonal odd: (0, 1,-1, 0).
```

Here `B1` is reflection-even and `B2` is reflection-odd for the declared
x-axis reflection.  There is no `A2` block.  More importantly, the balanced
axes-minus-diagonals H4 alias is only the second vector in a multiplicity-two
`A1` isotypic component.  D4 character theory alone cannot distinguish it
from the scalar; the zero-sum/augmentation condition supplies that extra
choice.

## Full modular image: one irreducible block

Adding `T` changes the useful classification:

```text
P1(F2): PSL(2,F2)=S3,  permutation representation = 1 + standard 2D;
P1(F3): PSL(2,F3)=A4,  permutation representation = 1 + standard 3D.
```

The exact character inner product of each nontrivial standard character is
one.  Thus the F2 H4 alias and axis-odd coordinates are two axes of a single
S3 irrep.  At F3, H4, axis-odd, and diagonal-odd are three axes of a single
A4 irrep; none is a standalone modular field.

In the orthonormal F3 basis `[H,A,D]`, the shear is

```text
T = [  0       1/sqrt(2)  -1/sqrt(2)]
    [1/sqrt(2)    1/2         1/2   ]
    [1/sqrt(2)   -1/2        -1/2   ].
```

Restricting this A4 triplet to `<T>=C3` gives a neutral line and a real
charged plane:

```text
u = H/sqrt(3) + sqrt(2/3) A,
v = sqrt(2/3) H - A/sqrt(3),
w = D,

T|_[v,w] = [-1/2 -sqrt(3)/2; sqrt(3)/2 -1/2].
```

Equivalently `v +/- i w` carry the two nontrivial C3 characters.  This is the
minimal charged/projective interpretation: diagonal-odd is one real
coordinate of a charged doublet, not a scalar anomaly.

## Placement of the existing N65 block

The same-modulus orientation contrast from `a7cb19a`, in the D4 basis, is

```text
[H,A,D] = [0.0011029842, 0.0051299851, 0.0020775613].
```

After the exact cyclic change of basis,

```text
[u,v,w] = [0.0048254235, -0.0020612155, 0.0020775613].
```

With the full shared-batch covariance, the charged pair `(v,w)` gives
`5.9050 / 2 df`; the neutral coordinate gives `z=1.602`.  The old marginal
`D` value remains approximately two sigma, but its correct unit of inference
is the correlated two-dimensional charged plane.  This is a representation
placement, not a discovery claim.

## Frozen no-fit prediction

If a future source is related by the declared projective shear and all source
normalizations are transported covariantly, the complete measured vector must
obey the exact matrix above.  The current vector therefore predicts

```text
T[H,A,D] = [0.0021583896, 0.0043837008,-0.0028238456]
```

with no fitted amplitude.  The strongest minimal next experiment is not a
second one-coordinate `D` score: it is a paired identity/shear source run that
freezes the full three-vector and tests this transport with one joint
covariance.  In a charged basis the sharper prediction is phase rotation by
`exp(+/- 2 pi i/3)` while `u` stays fixed.  Failure of that rotation after
declared source transport falsifies the projective-source interpretation;
success would connect the projective birth archive to the charged-sector
program without another harmonic vote.

## Reproduce

```bash
python3 scripts/flat_twist_representation.py \
  --score results/local-20260830/P334-projective-birth-N65-smoke/flat_twist_score.json \
  --json results/local-20260830/P334-projective-birth-N65-smoke/flat_twist_representation.json \
  --markdown results/local-20260830/P334-projective-birth-N65-smoke/flat_twist_representation.md

python3 -m unittest discover -s tests -p 'test_flat_twist_representation.py'
```
