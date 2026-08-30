# Q-lift covariance of the matching projector

Status: exact finite-volume algebra and exact `L=2,3` square-bond FK oracle.
This closes the first ambiguity identified in Issue #333; it does not yet
construct a canonical LCFT connection.

## 1. The endpoint observable does not determine its Q tangent

Write the restricted torus FK state sums as `W_0D`, `W_1D`, and `W_2D`.
Two natural continuations of the `Q=1` matching projector are

\[
 H_Q=W_{2D}-W_{0D},\qquad
 C_Q=W_{2D}-QW_{0D}.
\]

They agree at `Q=1`, but

\[
 \boxed{C_Q-H_Q=-(Q-1)W_{0D}}.
\]

For any declared path `Q -> (Q,v(Q))`, therefore,

\[
 (D_Q C-D_Q H)_{Q=1}=-W_{0D}(1,v(1)).
\]

After division by the full partition function this becomes exactly

\[
 (D_Q c-D_Q h)_{Q=1}=-\pi_{0D}.
\]

The individual derivatives depend on the path.  Their counterterm difference
does not, because the factor `Q-1` kills every chain-rule contribution from
the derivative of `W_0D` at the endpoint.

More generally a lift change

\[
 O_Q\longmapsto O_Q+(Q-1)X_Q
\]

shifts its first tangent by `X_1`.  Higher derivatives shift as
`D_Q^n O -> D_Q^n O+n D_Q^{n-1}X` at `Q=1`.  Thus an unqualified raw
`Q` jet is an affine, section-dependent object at every order.

## 2. Critical-polynomial lift is horizontal on the exact critical relation

Along a critical torus relation

\[
 W_{2D}(Q,v_c(Q))=QW_{0D}(Q,v_c(Q)),
\]

the two sections behave very differently:

\[
 C_Q=0,\qquad H_Q=(Q-1)W_{0D}.
\]

Hence

\[
 D_Q C|_1=0,\qquad D_Q H|_1=W_{0D}|_1,
\]

or, for normalized probabilities, `D_Q c=0` and `D_Q h=pi_0D`.

This reclassifies the proposed parameter-free homology tangent.  Relative to
the critical-polynomial/periodic-TL section, it is exactly the finite
transition function between two lifts.  It is not by itself a dynamical
response or a unique logarithmic-field insertion.  That does not make it
useless: it makes it a known connection coefficient that must be subtracted
or transported before tangents from #258, #262, #263, and #275 are compared.

## 3. Exact tiny-torus oracle

`scripts/q_lift_covariance_oracle.py` exhausts square-bond FK configurations,
retaining the repository's integer torus-homology rank and the exact sparse
polynomials in `(Q,v)`.

At `Q=v=1` it finds:

| L | configurations | W0D | W1D | W2D | pi0D |
|---:|---:|---:|---:|---:|---:|
| 2 | 256 | 69 | 118 | 69 | 69/256 |
| 3 | 262144 | 75460 | 111224 | 75460 | 18865/65536 |

For both sizes:

- the sparse polynomial identity `C-H=-(Q-1)W0D` is exact;
- at fixed `v=1`, the individual normalized tangents differ between `H` and
  `C`, with difference `-pi0D`;
- on `v=sqrt(Q)`, `D_Q c=0`, `D_Q h=pi0D`, and their difference is again
  `-pi0D`.

The complete sparse state sums and exact fractions are in
`results/q-lift-covariance/latest.json`.

## 4. Consequence for the research program

The minimum claim-bearing descriptor for a `Q` derivative must include the
generic-`Q` lift, normalization, and path.  Endpoint channel semantics alone
cannot type the tangent.

A productive next object is not another raw score but a transported one:

\[
 \nabla_Q O = D_Q O-\Gamma_Q[O],
\]

where `Gamma` is fixed by an exact section choice.  The CP section supplies a
particularly concrete candidate: declare the exact critical-polynomial
projector horizontal on the critical manifold.  The homology tangent then has
known connection term `pi_0D`.  A Vasseur--Jacobsen--Saleur energy/log-pair
positive control must decide whether the remaining transported residue is the
same physical logarithmic observable across lattice representatives.

## Boundary

This result proves lift ambiguity and an exact transition coefficient.  It
does not prove that the CP-horizontal connection is the unique physical
choice, nor identify the Matching-One H4 field with a logarithmic partner.

## Sources

- Jacobsen--Scullard, critical polynomial restricted state sums:
  https://arxiv.org/abs/1211.4335
- Jacobsen, periodic-TL critical-polynomial crossing:
  https://arxiv.org/abs/1507.03027
- Arguin, critical Potts/FK torus homology sectors:
  https://arxiv.org/abs/hep-th/0111193
- Vasseur--Jacobsen--Saleur, logarithmic observables at `Q -> 1`:
  https://arxiv.org/abs/1206.2312
