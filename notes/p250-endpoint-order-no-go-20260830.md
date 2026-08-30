# P250: endpoint moments cannot identify ordered morphisms

Date: 2026-08-30. Status: **exact non-identifiability theorem and physical-semantic gate**.

The current projective-leg archives have reached the limit of endpoint expansion.  More radii, higher total degree, and higher moments can refine the commutative endpoint module, but none can distinguish `T_x T_y` from `T_y T_x`.  The next useful object must retain a state after the first operation.

## 1. Endpoint-only no-go

Let `W=<x,y>` be the free word monoid and let

```text
ab(w) = (#x in w, #y in w) in N^2.
```

An endpoint-only archive has the form

```text
O(w) = Phi(G(ab(w))),
```

where `G(a,b)` may be scalar, vector, charged, multi-radius, or a joint moment row, and `Phi` may be nonlinear.  Since

```text
ab(xy)=ab(yx)=(1,1),
```

exactly

```text
O(xy)-O(yx)=0.
```

The same proof applies to every word pair with equal abelianization and to any finite collection or function of such endpoint rows.  Increasing insertion radius, total displacement degree, moment order, charge rows, covariance precision, or sample count cannot change this identifiability statement.

The certificate checks every word through length eight, radii 1--4 and moment orders 1--6.  More importantly, the result is an algebraic factorization theorem, not a bounded numerical observation.

### Actual runner gate

The physical implementation confirms the factorization.  `charged_rows(context,field,parent_indices)` receives final parent indices, constructs one static black-NN/white-matching component index, evaluates root membership, and applies the five-fiber DFT.  It has no path or first-step argument.  The exact gate checks `xy=yx` at every parent vertex in the checked N325 runner geometry.

This is compatible with the existing conclusion that the microscopic fivefold spatial cover is flat.  It also means that an endpoint Hankel failure is not evidence for noncommuting translations: it can reflect larger commuting state, nonreduced commutative quotient, periodic/context mixtures, or simply truncation.

## 2. What the P333 positive control actually proves

At three abstract connectivity marks, take

```text
D = formal detach of mark 1,
J = formal join of marks 0 and 1.
```

On the four-dimensional endpoint radical, exact rational algebra gives

```text
rank(DJ-JD)=2,
K=D+J-DJ-JD=(D-J)^2,
rank(K)=1,
K^2=0,
H K=K^T H.
```

Thus P333 supplies two related but distinct controls:

- `DJ-JD` proves that the two ordered histories are different on the marked state;
- `K` is a Gram-compatible connected signed rectangle with a nilpotent Jordan chain.

`K` uses both orders but sums them with the same sign.  It is not itself the order-antisymmetric observable.  A future runner should retain both

```text
R_plus  = L_D + L_J - L_DJ - L_JD,
R_minus = L_DJ - L_JD.
```

The first transports the P333 positive control; the second answers the P250 ordering question directly.

## 3. Why P333 cannot yet be embedded physically

The current projective-leg state at a root is only

```text
+1  black NN rank-one component,
-1  white matching rank-one component,
 0  otherwise,
```

followed by a Z5 fiber DFT.  This is not the three-mark set-partition module used by P333.

A formal rename is invalid for five independent reasons:

1. the runner does not emit the connectivity partition of three marked legs;
2. no first morphism is applied and no intermediate state is retained;
3. formal detach does not specify which occupied site, incident edges, or cut set is removed;
4. formal join does not specify an occupied connector and may merge colour states that no percolation move can join;
5. a physical site mutation changes black NN and complementary white matching connectivity simultaneously, while P333 acts on one uncoloured partition.

Therefore the exact decision is

```text
NOT_IMPLEMENTED_NO_PHYSICAL_SEMANTICS.
```

The P333 matrix remains a valid algebraic positive control.  It is not yet a projective-leg production observable.

## 4. Minimal connected two-morphism rectangle

For each common-random-number replica and one marked triple, the smallest sufficient online record is:

```text
base mark descriptor S0,
first operation tag and microscopic support,
intermediate marked descriptor after D,
intermediate marked descriptor after J,
responses L_D, L_J, L_DJ, L_JD.
```

The marked descriptor must contain, for every leg/block:

- black or white colour;
- component ID in the appropriate NN or matching graph;
- component ambient rank;
- primitive homology line when rank one.

The two morphisms must be defined as paired primal/matching mutations on the actual site configuration.  Their support cannot be inferred from a final displacement label.

The minimal publishable rectangle is then

```text
sample_id, marked_triple, morphism_D, morphism_J,
S_D, S_J, L_D, L_J, L_DJ, L_JD,
R_plus, R_minus.
```

Storing `S_D,S_J` is the essential addition.  Endpoint responses alone can show that two final answers differ, but cannot certify that the declared physical operations realized the intended detach/join morphisms.  Full microscopic paths are unnecessary if the operation support and these typed intermediate descriptors are retained.

No production runner is added here because the repository has not yet declared physical `D` and `J` mutations.  The next protocol should freeze those mutations first; only then should sampling begin.

## 5. Consequence for the current P250 program

The radius-six/rank-eight endpoint results already answer the commutative quotient question at their observed degree.  Another endpoint shell can improve a flatness bound, but it cannot test path order.  The next mechanism-changing acquisition should instead be the connected two-morphism rectangle above.

This does not claim path memory, noncommuting microscopic translations, a physical Jordan block, or a continuum field.  It proves why the current data cannot decide those propositions and specifies the minimum additional state needed to make the question identifiable.

## Reproduction

```sh
python3 scripts/p250_endpoint_order_no_go.py
python3 tests/test_p250_endpoint_order_no_go.py
```

Output: `results/exact-p250-endpoint-order-no-go/latest.json`.  Five focused tests pass; the inherited P333 and P262 exact suites also pass.

Related issue: #250.  P333 positive-control source: commit `6c60b0e`.
