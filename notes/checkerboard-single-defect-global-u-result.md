# One saturation defect rules out a scalar continuation of the closed source

## The decision

The exact closed-source endpoint does **not** extend even infinitesimally
as a source-independent overall gain of the original global U. The
predeclared necessary identity

`R = U U_st - U_s U_t = 0`

is violated by the complete N50 single-defect calculation:

`R = +27.766563581230237`.

The exact rational enclosure of R/A50^2 is strictly positive. This is a
finite model elimination, not an extra significance vote or a fitted
correction to the previous endpoint result.
The zero decision does not rely on choosing the exponent13/8: replacing
the fixed nonzero area/projector normalization scales R by its square.

## Fixed physical experiment

Use the same source throughout:

`S*=C+F+Bvac,  weight proportional to P_(p,s) exp(t S*)`.

The A/B probabilities are `pA=s+(1-s)p,pB=p`. Thermal derivatives hold
s,t fixed, and each U(s,t) uses its own pooled q=0 root and slope ratio
`U=A50 P4(E_p)/mean(q_p)`. The parent pair is `(5,5),(1,7)`, with
delta_cos4=-1152/625. At s=1 it maps to the N25 `(5,0),(4,3)` pair.

The new population removes the origin A site, keeps the remaining24 A
sites occupied and sums every one of the2^25 B configurations in each
geometry. Translation supplies the other24 defect positions. The intact
endpoint is taken from the existing exact child coefficients by complement;
it was not re-enumerated or sampled.

The source and the physical saturation chart are fixed before calculation.
The [freeze9024fdbf](checkerboard-endpoint-defect-decision-freeze.md) is
earlier than the new enumeration. The producer-only contract6c65157f
initially called s the A occupation probability; the integer conditional
counts do not depend on that label. Before scoring,6b5e66bc clarified it
to the physical chart already specified in9024fdbf. The original receipts
remain unchanged. No observable, graph, coefficient or decision was changed.

## Measured coefficients and the eliminated prediction

| Quantity at s=1,t=0 | Exact-coefficient numerical value |
|---|---:|
| Parent pooled root | 0.40733446067177326 |
| U | 2.715728877348466 |
| U_s | 3.708240929282322 |
| U_t | 0.38914717849771724 |
| U_st under source-independent gain: U_s U_t/U | 0.5313680267777353 |
| Actual U_st | 10.755718407564073 |
| Actual minus gain prediction | 10.224350380786337 |

The unknown geometric gain slope cancels in the decision R; the displayed
prediction merely makes that algebraic elimination interpretable. There
is no gain-fit uncertainty because these are complete finite sums, and
no free coefficient is adjusted after the result.

Equivalently,

`partial_s(U_t/U)=+3.764864182896271`.

Therefore moving **away** from saturation, epsilon=1-s, decreases the
relative source response at this endpoint. The source is changing its
coupling to the global observer, beyond the change in U's overall amplitude.
The exact endpoint closure and same-t transport still hold at s=1.

The overview's separately stated mixed thermal-only null is resolved too:

`Xi=U_t,epsilon=-U_st=-10.755718407564073 != 0`.

These two decisions come from the same finite calculation. They are not
independent confirmations, and the scalar-gain rejection is the stronger
predeclared amplitude comparison here.

In the original source units, the computed double first jet is

```text
U(1-epsilon,t) = 2.715728877348466
                -3.708240929282322 epsilon
                +0.38914717849771724 t
                -10.755718407564073 epsilon*t
                +O(epsilon^2)+O(t^2).
```

This is a local expansion derived from complete coefficients, not a fitted
surface. In particular it does not assert a finite-epsilon sign change.

## What the defect measures

The actual source change is a cycle-loss and ambient-rank-loss insertion.
Let k be the lost occupied graph-cycle dimension, ell the lost ambient
rank, and k_null=k-ell. Since one site is removed,

```text
Delta S* = 3-2k+ell = 3-2k_null-ell,
q_defect = q_intact-ell,
E_defect-E_intact = ell^2-2ell*q_intact.
```

These follow directly from S*=2beta1-3K-q+2N and E=q^2. The graph and
ambient-image inclusions give k_null>=0. Thus the elementary insertion
has a specified topological meaning; no fitted contact feature is added.
The [single-defect proof](checkerboard-single-defect-source.md) narrows the
mechanism further at this saturated endpoint: **ell is at most1, and only
an alternating child face can change rank**. All other neighbor patterns
have local lifted bypasses or cannot carry a cycle. The alternating case
switches the local diagonal connection; exterior winding decides whether
the global rank changes.

At zero source only these rank-changing switches directly alter q/E.
At nonzero source, rank-preserving holes can also enter through the exact
weight-redistribution term `Cov(exp(t Delta S),O_intact)`. The present
mixed result includes both terms, without claiming a separately measured
share for either. This is a two-term insertion identity, not a proposal
to add another empirical descriptor catalogue.
The defect must be evaluated on the parent graph: replacing it by an
ordinary independent child site would remove the four-terminal surgery
that this question is testing.

For unnormalized parent moments, the entire required s jet is

`H_O=H0_O+(1-s)*25*(1-p)*(Hd_O-H0_O)+O((1-s)^2)`.

Divide by each geometry's own H_1 **before** pooling. The mixed source
normalization, the p derivative of1-p, the moving root and the changing
thermal-slope denominator all contribute. The [complete jet derivation](checkerboard-endpoint-gain-discriminant.md)
shows that only first source moments and one defect are necessary; no
second defect or S*^2 moment is needed for this decision.

The fixed S* cannot be exchanged for a density-shifted action while
retaining this saturation chart. Such an exchange preserves the hard
endpoint U(t), but mixes saturation and temperature in its neighborhood
and can change R. This is why the same source units are retained throughout.

## Provenance and scientific boundary

- **Lifecycle:** freeze9024fdbf; producer contract6c65157f/code393ea7c4;
  committed exact counts13a2c197; root scorer695cadba/ebc3b790;
  one score at6b5e66bc, resultf5c4a74a. The outcome-blind theory companions
  are1df9b35d (mixed derivative) and4981e625 (cycle/rank defect).
  All branch-delivered, not main.
- **Cost:** two enumerations1.020/1.090seconds; compilation and both runs
  together2.497seconds; one rational score0.300seconds. Local Mac only,
  no Monte Carlo, cloud job, baseline replay or test campaign.
- **Observer/sector/source/geometry:** original pooled-root/slope U,
  ordinary homology q/E with normalized cos4 contrast, bulk S*, N50
  parents of the declared N25 pair. The two quotient Smith classes differ.
- **Dependency:** one exact defect block plus its pinned exact N25 endpoint;
  all derivatives, residuals and predictions share these coefficients.
- **Files:** [complete result](../results/p337-endpoint-defect/score/REPORT.md),
  rational enclosures and every ratio/root term in `score/score.json`,
  integer profiles and producer receipt in the parent result directory.
- **Boundary:** this is the one-sided endpoint first jet. It does not give
  a finite-epsilon zero crossing, the full interior curve, the homogeneous
  large-size law, a continuum H4 field or an RG semigroup. Linear
  extrapolation of the displayed large mixed slope would not establish
  any of those. P154/P334/F4 stop decisions remain unchanged.

The new requirement on any proposed interior theory is explicit: it must
carry a source-dependent defect correction to the original U, with this
fixed finite coefficient, instead of extending the endpoint by a scalar
geometric gain. No extra source, descriptor or automatic production round
is introduced to make a failed model survive.
