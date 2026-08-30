# P398 canonical Q-adic/Jantzen control at widths three and four

The unmarked width-three and width-four Gram degenerations have Loewy/Jantzen length two only: all positive local invariant-factor valuations equal one and J2 vanishes. Width three's tested scalar plus C3-charge-one responses span the entire grade-one layer. At width four the tested scalar, C4-charge-one and C4-charge-two responses cover only 6 of 13 grade-one dimensions, leaving exact deficits 2, 2 and 3 in the trivial, charge-one and charge-two sectors. The missing datum is therefore multiplicity/rooted-connectivity within existing sectors, not a deeper Q-adic order or another terminal C4 irrep.

| width | dim V | valuations | dim gr0 | dim gr1 | dim J2 | det leading radical form | tested rank | uncovered |
|---:|---:|---|---:|---:|---:|---:|---:|---:|
| 3 | 5 | `0^1,1^4` | 1 | 4 | 0 | -1 | 4 | 0 |
| 4 | 14 | `0^1,1^13` | 1 | 13 | 0 | -1 | 6 | 7 |

## Associated-grade coverage

### Width 3

| sector | gr1 dimension | tested projection rank | uncovered |
|---|---:|---:|---:|
| trivial | 2 | 2 | 0 |
| charge1_rational | 2 | 2 | 0 |

### Width 4

| sector | gr1 dimension | tested projection rank | uncovered |
|---|---:|---:|---:|
| trivial | 5 | 3 | 2 |
| charge1_rational | 4 | 2 | 2 |
| charge2 | 4 | 1 | 3 |

## Interpretation

- The full polynomial Gram family is available through degree `width`, but the exact unimodular leading radical form proves that first order already identifies every nonzero layer.
- At width three the tested responses cover the complete grade-one representation. At width four every C4 sector still has uncovered multiplicity, so changing only the terminal character cannot close the gap.
- The automatic dual-number `t` action is square-zero before specialization and zero after `Q=1`; it is a base-parameter nilpotent, not a fixed-Q marked or physical extension.
- Nonzero associated-grade projection is necessary descriptive information but not a closure certificate: the existing width-four fixed-Q tests fail despite nonzero projection in every tested sector.

## Boundary

- Exact local-Q Gram algebra for the repository's declared unmarked connectivity modules at widths three and four only.
- The associated-graded projections classify already-tested response covectors; they do not construct a fixed-Q extension.
- J2=0 rules out a deeper base-Q filtration in these modules, not nonsemisimple rooted, direct-sum or physical transfer modules.
- No continuum LCFT field, periodic-TL cell-module dictionary or physical Jordan block is identified.
