# Scientific-frontier delta after the PR267 snapshot

This note is the compact handoff for scientific work that landed after the
immutable PR267 content snapshot `3080d8bf`.  Every item below is pinned to a
commit.  A citation records provenance; it does not promote a branch result to
`main`.

## One-screen mechanism update

| Line | New result | What it removes | Live discriminator |
|---|---|---|---|
| P250 | A model-free degree-two Hankel score keeps separate plus/minus two-charge rank-five charts compatible while rejecting a raw shared rank-five chart.  A later null-direction comparison leaves both Alexander reflection-conjugation and identity maps compatible, but does not jointly enforce the two annihilation residuals. | Common raw state dimension at most five; the claim that one fitted Weyl spectrum explains the split.  It does not yet remove the identity map or identify an Alexander bridge. | Jointly score `H_plus q=0` and `H_minus A(q)=0` for every frozen map using the same 400-batch covariance; only then consider a radius-five shell. |
| P333 | The universal first-jet connectivity radical now has an exact Gram-compatible Jordan positive control, `K=D+J-DJ-JD=(D-J)^2`, on three marks. | Join-only dynamics, a single deterministic detach/join history, and the idea that an epsilon pivot alone selects a log pair. | Embed the signed connected-history subtraction in a physical Q-dependent transfer expansion and test the dimension-velocity collision separately. |
| P334 | Translation regularity collapses every proper Hall cut; aggregate TM is exactly `4DF<=M^2+4Y(T-D)`, and the finite typed rewrite has one unclosed global `D x F` critical pair. | Larger Hall atlases, displacement-local, quotient-order and Fourier-SOS proof routes. | Realize `D x F -> M x M` or `Y x nonD` globally across displacement classes; pursue BA independently. |
| P337 | The preregistered N680 child lands closest to the two-mode recurrence and excludes scale-neutral, but remains unresolved against free-single.  Four-generation recurrence, fixed identity dressing and same-base Jordan all survive; N1360 is source-covariance limited across them. | Fixed nominal single H4 and no decay on this lineage; not the live correction mechanisms. | Decompose all four archived generations into first-birth, completion and direct-rank-two recurrence components; if still aliased, add a geometry/modulus covector rather than another size. |

These moves change the project's center of gravity.  It is no longer useful to
say merely that a compact rank-two or rank-three state “survives.”  Survival is
word-, observer- and geometry-dependent.  The axis P250 recurrence closes at
rank three, while a statistically independent bivariate block excludes a
common diagonal rank at most three and then excludes the two frozen minimal
rank-five families.  The next state description must say which words are
observed, whether path order is retained, which source/readout projection is
being represented and whether a claimed sector map was tested against the
full annihilation residual or only an extracted singular-vector direction.

## P250: no minimal Weyl state

The fresh bivariate block is pinned at
`4c1f5d8cf98b9c39557acff0146332899777505b`.  It contains 80,000 samples in
400 aligned batches and is independent of the earlier N505 axis stream.  The
signed C4 closure passes (`66.858/64`, `p=.379`), but common commuting diagonal
rank three fails both its mixed gate (`66.981/24`, `p=6.21e-6`) and the
degree-four holdout (`141.732/40`, `p=2.69e-13`).

The post-result hypothesis and scorer were frozen at `6065312`; the reveal is
pinned at `ded78bb94e7051b6cecbb4f4328d68f5b4cf7655`.  The exact spatial cover
has zero plaquette curvature at all 101 parent sites in both hands.  Every
nonzero-flux Weyl model is therefore only an effective projected-state model,
not literal deck holonomy.

| Frozen model | Mixed score | Degree-four holdout |
|---|---:|---:|
| Weyl `(1,4)` | `1390.04/24`, `p=6.68e-279` | `342.24/40`, `p=1.21e-49` |
| Weyl `(2,3)` | `1428.71/24`, `p=3.61e-287` | `312.81/40`, `p=5.44e-44` |
| Weyl `(3,2)` | `1939.97/24`, `p<1e-300` | `340.89/40`, `p=2.20e-49` |
| Weyl `(4,1)` | `1269.17/24`, `p=4.33e-253` | `329.16/40`, `p=4.01e-47` |
| free commuting rank four | training includes mixed degree <=3 | `107.86/40`, `p=3.81e-8` |
| free commuting rank five | training includes mixed degree <=3 | `74.90/40`, `p=6.82e-4` |

The free rank-five result is much closer than the Weyl models, so the data do
not positively select noncommutation.  They redirect attention toward higher
dimension, non-normal/Jordan transfer, finite periodic images or a
context-dependent projection.  The archive does not record `TxTy` versus
`TyTx` path order; the model-fixed center `D` is not an independently observed
coordinate.  Those two facts must remain visible in every interpretation.

The invariant-first Hankel result at
`a770ac9f71460564dd2090210e1e0e64a1cd4979` removes the exponential-root
assumption.  For the two-charge blocks, rank at most four fails within each
hand (`p=2.80e-11` plus, `.00235` minus), while rank at most five remains
compatible at the frozen `.01` level (`p=.0543/.0655`).  The raw shared block
rejects rank at most five (`p=2.84e-5`).  This means “compatible truncated
five-dimensional chart in each hand,” not exact rank five.

The subsequent direction-only bridge result
`a46ed6343b4b6874259a0f24defeedf16800a8d0` keeps all four frozen Alexander
reflection-plus-conjugation maps, but also keeps identity comparators.  It
compares separately extracted smallest-singular-vector lines and does not
jointly require `H_plus q=0` and `H_minus A(q)=0`.  It therefore cannot by
itself reinterpret the shared-rank rejection as a coordinate mismatch or
identify an Alexander intertwiner.  The clean zero-production repair is a
candidate-constrained joint-null score with the complete paired covariance.

## P333: a minimal weighted-history Jordan mechanism

At `6c60b0e06360b035ed0cc54df230fe57fd6a7549`, the exact three-mark radical
contains

```text
K = D + J - D J - J D = (D-J)^2,
rank K = 1,  K^2 = 0,  H K = K^T H.
```

With `v=(1,-1,-1,0)` and `w=(0,0,0,1)`, `Kw=v`, `Kv=0`,
`<v,v>_H=0` and `<v,w>_H=1`.  This is the first finite positive control in
this line that has both a Jordan chain and first-jet Gram compatibility.  The
signed subtraction is not a positive stochastic transfer and is not an LCFT
field identification.  The next problem is no longer “find a generator”; it
is “derive this connected-history subtraction from a physical Q-dependent
transfer or prove that the physical cone cannot contain it.”

## P334: the Hall family has collapsed

At `c8405b59cc289c0b607f7f78a3bb9f513be4cef6`, torus translation symmetry
makes the demand and supply pair graphs regular and proves

```text
g(U) = |U|/N * g(V) + (D_cut(U)+S_cut(U))/2.
```

The all-site aggregate TM inequality therefore implies every proper Hall cut;
2,470,440 bounded cuts are exact regression checks, not additional theorem
assumptions.  In all 688 positive-demand rows the all-site cut is the unique
worst ratio.  Only aggregate TM and the independent BA concordance inequality
remain to be derived from arbitrary digital Alexander quotients.

The later branch tip `7ef99ae48a354c3a2199075746feef8081e637d0`
reduces aggregate TM to one four-face inequality,
`4DF<=M^2+4Y(T-D)`.  Ordinary Rayleigh and synergy-only coverage each fail on
bounded rows, and the inequality fails displacement by displacement, after
quotient-order grouping and under a Fourier-SOS attempt.  A terminating
count-level rewrite closes all 984 atlas rows but is not a configuration
injection.  The sole general TM problem is the global, collision-free critical
pair `D x F -> M x M` or `Y x nonD`; BA remains a separate marked-pair
concordance problem.

## P337: H4 sign is stable, radial transfer is curved

- `0db21b7c21fe7912b28804582411dc5966dd957d`: N170 gives
  `Delta K_A=-.0177169+/-.0015556` (`z=-11.389`).  The residual to the frozen
  H4 magnitude is `-3.035` predictive SE and remains in the H4 direction;
  the projective scalar is null (`z=-.783`).
- `e819f5ecd0074ee2250c1852457fd8dc758099ea`: N340 gives
  `A_H4=-.0048573+/-.0012489`; the sign flips back (`z=3.889`), the
  scale-neutral target is excluded by `3.946` predictive SE, and the
  projective scalar remains null (`z=.506`).
- `4024a7c381d7c269880734bbe61355aa8a4232c2`: fixing the leading transfer to
  `2^-13/8` yields `lambda1=.212+/-.236`.  This zero-residual three-point
  interpolation is a useful coordinate, not a resolved second eigenvalue.
- `02080a4c574914439e080299ad28ab66f5172826`: the preregistered N680 child
  gives `A_H=-.00216756+/-.00055693`.  It lands at `-.295` predictive SE from
  two-mode, `+.747` from free-single, `-2.024` from fixed single and `+7.159`
  from scale-neutral; its projective scalar is null.  The recorded historical
  server is outside the permitted environment and is not contacted here; only
  the committed compressed archive is used as provenance.
- `7263fbc49d641f9b20af8f19a5aa46280cda6e2c`: the four-generation fit gives
  `lambda1=.27068+/-.12801`, recurrence `q=.0773/1` (`p=.781`) and free-single
  `q=1.979/2` (`p=.372`).  Their descriptive AIC differs by only `.098`; the
  secondary mode is boundary-marginal, not a discovered eigenvalue.
- `efd6d31`: the theory-fixed adversary uses no new samples.  Fixed identity
  dressing with eigenvalues `2^-13/8` and `2^-21/8` passes (`.803/2`,
  `p=.669`) and leads descriptive AIC; a same-base rank-three Jordan
  polynomial also passes (`.084/1`, `p=.772`).  Source covariance makes a
  lone N1360 child unable to separate these mechanisms universally.

## Default attention after this delta

1. Repair the P250 bridge comparison with one candidate-constrained joint-null
   score on the existing 400 batches.  Add a radius-five shell only for maps
   that survive; ordered words remain later.
2. Transport P333's exact signed connected-history nilpotent into a physical
   Q-dependent transfer expansion; keep positivity and dimension velocity as
   separate checks.
3. Realize P334's unique global `D x F` critical pair across displacement
   classes and pursue BA as an independent concordance injection.
4. Use the committed N85/N170/N340/N680 birth archives to locate the P337
   correction in first activation, second completion or `DIRECT_RANK2`.
   Retain N130/N170 as a cross-lineage control and missing path fields as
   `not_scoreable`; if aggregate aliases survive, add a geometry/modulus
   covector.  Do not contact the historical N680 server or assume N1360.
5. Use proof-carrying model elimination where it cheaply converts a declared
   low-rank class into a verifiable no-go.  An inconclusive relaxation does not
   demote the mechanism or block exploratory work.

This order allocates attention only.  All lines remain open to exact work,
reuse, pilots and alternative conjectures.
