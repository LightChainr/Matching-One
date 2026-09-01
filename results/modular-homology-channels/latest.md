# Modular homology-channel oracle

Exact integer-homology classification; no Monte Carlo or fitted continuum field is used.

| channel | classification | reason |
|---|---|---|
| `rank` | modular scalar | rank(MW)=rank(W) because M in SL(2,Z) is invertible over Q |
| `either` | modular scalar | either is exactly rank>0 |
| `cross` | modular scalar | cross is exactly rank=2 |
| `direction_0` | basis-dependent | depends on the selected homology generator |
| `direction_1` | basis-dependent | depends on the selected homology generator |
| `both` | basis-dependent | a rank-1 spiral can use both generators and shear to one generator |

The scalar classification lifts to primal, matching, even, and odd combinations of `cross` or `either`.
Complement commutes with geometric relabelling, and linear combinations of scalar channels remain scalar.

## Exact counterexample

The rank-1 spiral basis `(1,1)` has `both=true`.  The determinant-one shear
`[[1,-1],[0,1]]` maps it to `(0,1)`, so `direction_0` and `both` become false while
`rank/either/cross` stay unchanged.  The transpose shear similarly changes `direction_1`.

## Scalar-channel elliptic filter

| spin | tau=i | hexagonal rho | square and hexagonal |
|---:|:---:|:---:|:---:|
| H0 | yes | yes | yes |
| H4 | yes | no | no |
| H8 | yes | no | no |
| H12 | yes | yes | yes |
| H16 | yes | no | no |
| H20 | yes | no | no |
| H24 | yes | yes | yes |

Thus a homogeneous first-order response in scalar `cross/either` channels kills H4/H8 at the hexagonal point and permits H12.
This conclusion does not apply to a primitive winding character, `both`, or another vector-valued channel.

## Boundary

- This proves a topological channel-label statement under SL(2,Z) basis change, not an FK-to-CFT matrix element.
- The elliptic zero additionally assumes a homogeneous first-order spin response; logarithmic or mixed-spin terms must be typed separately.
- Primitive winding characters, direction channels, and the repository `both` channel are vector-valued or basis-dependent and do not inherit the scalar zero.
- No overlap or selection rule between thermal Q4 and V_(2,2) is inferred, and Issue 114 remains open.
