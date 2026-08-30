# P250 result: the minimal magnetic-translation state is not present

The existing fresh N505 bivariate data reject all four frozen canonical Z5
Weyl realizations.  The result is not a positive observation of projective
translation; it sharply moves the common-state search away from the minimal
clock/shift model.

## Exact result

The unit x/y lifts commute around every elementary plaquette in the actual
fivefold cover.  Both plus and minus have curvature count
`{0:101,1:0,2:0,3:0,4:0}`.  Thus the microscopic cover has center charge zero.
The nonzero `m` models tested below were, by construction, only candidate
effective representations after observable projection.

The canonical matrices themselves pass all relations to maximum residual
`3.42e-15`, including `R X R^-1=Z_m`, `R Z_m R^-1=X^-1`, `R^4=I`, the Weyl
commutator, and `D^5=I`.

## Frozen result

| model | mixed gate | degree-four heldout |
|---|---:|---:|
| Weyl `(m_plus,m_minus)=(1,4)` | `1390.04/24`, `p=6.68e-279` | `342.24/40`, `p=1.21e-49` |
| Weyl `(2,3)` | `1428.71/24`, `p=3.61e-287` | `312.81/40`, `p=5.44e-44` |
| Weyl `(3,2)` | `1939.97/24`, `p<1e-300` | `340.89/40`, `p=2.20e-49` |
| Weyl `(4,1)` | `1269.17/24`, `p=4.33e-253` | `329.16/40`, `p=4.01e-47` |
| free commuting rank 4 | used for training | `107.86/40`, `p=3.81e-8` |
| free commuting rank 5 | used for training | `74.90/40`, `p=6.82e-4` |

The free commuting controls were given every first-quadrant point through
degree three, including mixed points, while Weyl received axes only.  Even
with that favorable advantage rank five misses the frozen `p>=0.01` gate.
Nevertheless it is enormously closer than any Weyl model.  The data therefore
do not support noncommutation: the current information gradient points first
to a larger or non-normal commuting state, or to a periodic-image/context
kernel, rather than to canonical magnetic translation.

The rank-five rejection is multivariate rather than one visibly bad cell.  Its
largest marginal degree-four residual is only `1.56` standard errors; the
joint failure is distributed across several low-variance covariance modes.
That makes the next useful object a covariance-aware model-free bivariate
Hankel rank lower bound, not another hand-picked transfer spectrum.

## Scientific card

- **Mechanism space changed:** the minimal irreducible Z5 clock/shift state is
  removed, and a free diagonal state of dimension five is still insufficient
  under the frozen realization.
- **Not proved:** nonexistence of all rank-five projective states, a unique
  rank-six state, Jordan structure, or path/context memory.
- **Observer/sector/source/geometry:** neutral projective-leg charged pair;
  charges 1/2, plus/minus Gaussian children; fresh N505 radius-four diamond.
- **Dependency group:** reuses the 80k/400-batch block with seed
  `25050510120261130`; no new random data.
- **Next discriminator:** the degree-two monomial 6x6 multivariate Hankel
  matrix and its full-covariance low-rank distance, including Jordan-compatible
  commuting realizations and an explicit flat-extension consistency test.
