# Scientific-frontier delta after the PR267 snapshot

This note is the compact handoff for scientific work that landed after the
immutable PR267 content snapshot `3080d8bf`.  Every item below is pinned to a
commit.  A citation records provenance; it does not promote a branch result to
`main`.

## One-screen mechanism update

| Line | New result | What it removes | Live discriminator |
|---|---|---|---|
| P250 | Fresh N505 bivariate data preserve the exact signed C4 envelope but reject a common commuting diagonal state through rank three.  A frozen reuse then rejects all four canonical five-dimensional Z5 Weyl models and the favorable free commuting rank-four and rank-five controls. | A two-dimensional three-state continuation of the axis recurrence; the minimal clock/shift explanation; the frozen shared-eigenpair rank-five explanation. | A covariance-aware model-free bivariate Hankel rank lower bound with a flat-extension check, followed by an ordered-word/path observation only if rank alone cannot decide. |
| P333 | The universal first-jet connectivity radical now has an exact Gram-compatible Jordan positive control, `K=D+J-DJ-JD=(D-J)^2`, on three marks. | Join-only dynamics, a single deterministic detach/join history, and the idea that an epsilon pivot alone selects a log pair. | Embed the signed connected-history subtraction in a physical Q-dependent transfer expansion and test the dimension-velocity collision separately. |
| P334 | Translation regularity gives `g(U)=|U|/N g(V)+(D_cut(U)+S_cut(U))/2`; aggregate TM implies every proper Hall cut. | Larger proper-cut atlases, line-coset classifications, and an independent Hall-family conjecture. | Prove aggregate TM and the separate BA concordance inequality from digital Alexander homology geometry. |
| P337 | N170 flips with exact H4 sign but exceeds the frozen magnitude; N340 flips back and returns toward nominal H4 while the projective scalar remains null.  A three-generation two-mode coordinate is plausible but unresolved. | Scale-neutral charged amplitude; another scalar/projective common mode as the N170 residual. | The frozen N680 same-lineage forecast, rerouted to the only permitted TV2N0X environment before any production, or a zero-production N130/N170 common-subset current comparison. |

These moves change the project's center of gravity.  It is no longer useful to
say merely that a compact rank-two or rank-three state “survives.”  Survival is
word-, observer- and geometry-dependent.  The axis P250 recurrence closes at
rank three, while a statistically independent bivariate block excludes a
common diagonal rank at most three and then excludes the two frozen minimal
rank-five families.  The next state description must say which words are
observed, whether path order is retained, and which source/readout projection
is being represented.

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
- `ba4ca6f87c7fa97fd6a2d878de83f03e4309b209`: the N680 forecast is frozen,
  but its old infrastructure binding is not executable under the current
  server constraint.  Scientific forecasts may be copied to a new manifest;
  the old server must not be contacted.

## Default attention after this delta

1. Build the model-free P250 bivariate Hankel lower bound before inventing a
   sixth hand-picked spectrum.  Add ordered words only if rank/flat-extension
   information cannot distinguish the survivors.
2. Transport P333's exact signed connected-history nilpotent into a physical
   Q-dependent transfer expansion; keep positivity and dimension velocity as
   separate checks.
3. Prove P334 aggregate TM and BA from topology.  Do not spend on larger Hall
   cut scans.
4. Decide whether the N680 forecast is worth a TV2N0X reroute.  In parallel,
   reuse the N130/N170 common archive subset and label the missing completion
   winding/lift/transporter fields as not recoverable.
5. Use proof-carrying model elimination where it cheaply converts a declared
   low-rank class into a verifiable no-go.  An inconclusive relaxation does not
   demote the mechanism or block exploratory work.

This order allocates attention only.  All lines remain open to exact work,
reuse, pilots and alternative conjectures.
