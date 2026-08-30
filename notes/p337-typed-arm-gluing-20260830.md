# P337/P334: the reverse comparison needs globally typed arms

Date: 2026-08-30. Status: **deterministic gluing theorem with exact fixtures; no exponent claim**.

The direct-birth carrier theorem gives

```text
direct theta birth  =>  at least 3 occupied NN + 3 vacant matching arms,
direct figure-eight =>  at least 4 occupied NN + 4 vacant matching arms.
```

This note determines exactly what is needed in the reverse direction.  Ordinary arm colours and landing locations are insufficient.  Once the event retains a small global homology type, a bounded local surgery gives a scale-independent finite-energy comparison.

## 1. Ordinary six arms are insufficient

Take a square torus with side `L>2R+2`.  From the vacant origin draw three occupied spokes to radius `R`: east, north and west.  Keep the remaining sites vacant.  Three matching-lattice vacant arms run through the northeast diagonal, northwest diagonal and south ray.  The six arms are separated and the radius-R box is embedded.

Before adding the origin there are three contractible occupied components.  Adding the origin joins them into one tree.  The ambient rank remains zero.  Thus even perfectly separated ordinary six arms do not imply a direct birth.

The checked fixture uses `L=13,R=4`; its before/after tuple is

```text
(components, ambient rank): (3,0) -> (1,0).
```

This counterexample persists at every scale.  No untyped six-arm probability can be reverse-comparable to the direct-birth probability by local finite energy alone.

## 2. Minimal globally typed event

Fix a lift of the candidate birth site `v`, an embedded annulus, and a fixed finite inner surgery box `Q`.  A **theta-typed six-arm event** records:

1. the exterior occupied configuration has ambient rank zero;
2. three occupied NN arms land in a fixed cyclic landing word and are separated by three vacant NN+NNN matching arms;
3. the three occupied landing sites belong to one projected exterior component `C`;
4. relative to one lift of `C`, the locally chosen lifts of the three landings have deck addresses `lambda_1,lambda_2,lambda_3` with

```text
det(lambda_2-lambda_1, lambda_3-lambda_1) != 0.
```

Only the relative address vectors are needed; absolute addresses are gauge.  The four fields above are sufficient.  Dropping any of the first three global fields admits an immediate obstruction:

- without exterior rank zero, the event need not be a `0->2` transition;
- without the common component ID, the three spokes counterexample applies;
- without the nonzero determinant, the birth creates rank at most one;
- without a fixed-width landing word, a bounded template need not be routable without a separate arm-separation input.

The last field is geometric rather than topological.  It can be replaced by any finite family of uniformly routable landing words plus a deterministic canonical choice.

## 3. Bounded theta surgery

Use `Q=[-2,2]^2`.  The birth site `(0,0)` is withheld.  Put exterior gates at

```text
g_E=(3,0), g_N=(0,3), g_W=(-3,0).
```

Force the six sites on the three two-step axial corridors from these gates to the east, north and west neighbours of `v` occupied.  Force the other 18 noncentral sites of `Q` vacant.

Each corridor attaches to the exterior at exactly one gate.  The corridors have no NN contacts with one another, so they are three pendant trees.  Adding them cannot change the zero ambient rank of the exterior.  After `v` is inserted, the new cycles have deck displacements

```text
lambda_N-lambda_E, lambda_W-lambda_E,
```

whose determinant is nonzero.  The transition is therefore a direct one-carrier theta birth.  Forced vacant sites continue the three exterior matching separators through the surgery box.

Under independent Bernoulli site measure, conditional on the exterior, the exact local pattern has probability

```text
c_theta(p) = p^6 (1-p)^18.
```

For `p in [eta,1-eta]`, this is at least `eta^24`, independent of torus size and arm radius.

### Finite landing families

Let `A_theta_sep` be a union of finitely many routable fixed-width landing words.  Choose a landing word canonically from the exterior configuration and apply its precomputed local template.  If every template forces at most `K` sites, conditional independence gives

```text
P_p(direct theta at v) >= eta^K P_p(A_theta_sep).
```

There is no many-to-one loss: the canonical word and template depend only on the exterior, and the inequality is obtained by conditioning rather than by counting images.

If `A_theta` is the unrestricted globally typed event and a uniform landing-separation estimate supplies

```text
P_p(A_theta_sep) >= delta P_p(A_theta),
```

then, together with the deterministic forward inclusion,

```text
eta^K delta P_p(A_theta)
    <= P_p(direct theta at v)
    <= P_p(A_theta).
```

The bounded surgery is proved here.  A uniform positive `delta` is a separate probabilistic arm-separation/nondegeneracy input; it is not silently assumed.

## 4. Figure-eight is a typed eight-arm subchannel

For four cyclic gates, record the two-pair component partition

```text
(g_E,g_W) in C_x,  (g_N,g_S) in C_y,  C_x != C_y,
```

and one nonzero deck difference from each pair.  Require the two differences to have nonzero determinant.  Four occupied arms and four vacant separators give the globally typed eight-arm event `A_8^fig`.

In the same `5x5` block, force eight axial corridor sites occupied and the other 16 noncentral sites vacant.  Its exact cost is

```text
c_figure(p) = p^8 (1-p)^16 >= eta^24.
```

The predecessor consists of pendant paths on two rank-zero components.  Adding `v` closes one loop in each component, with independent deck directions, and gives a direct figure-eight birth.

The exact extra-arm upper bound is

```text
P_p(direct figure-eight at v) <= P_p(A_8^fig) <= P_p(A_6 ordinary).
```

With a figure-eight landing-separation constant, the first inequality also has the bounded reverse comparison just proved.  No deterministic argument makes `P(A_8^fig)/P(A_6)` vanish.  Any asymptotic suppression of the figure-eight share needs genuinely probabilistic multi-arm input.

## 5. Exact fixtures

The certificate checks three independent pieces.

1. `L=13,R=4`: three occupied and three matching-vacant arms, but adding the origin leaves rank zero.
2. Gaussian quotient `3+i`, `N=10`, old mask `122`, birth site `7`: exact one-carrier theta, rank `0->2`.
3. Gaussian quotient `3`, `N=9`, old mask `30`, birth site `0`: exact two-carrier figure-eight, rank `0->2`.

It also checks that the theta and figure-eight local blocks have respectively `(6 open,18 closed)` and `(8 open,16 closed)` sampled sites, that their occupied corridors are pairwise NN-disconnected, and that their vacant separator seeds are valid matching-lattice paths.

## 6. Consequence for `D_N`

The missing bridge is now narrow:

- the topology and bounded local surgery are exact;
- the correct comparison object is not an ordinary six-arm event, but a globally typed theta-six event plus a globally typed figure-eight-eight event;
- the remaining nonlocal input is the frequency with which separated arms acquire the required component partition and rank-two deck landing matrix.

This global typing can be measured directly in production by storing, at the inner/outer annulus cuts, the occupied landing component IDs and two relative deck-address vectors.  Full microscopic paths are unnecessary.

No universality transfer, six/eight-arm exponent, nonzero continuum amplitude, or asymptotic `D_N` law is asserted.

## Reproduction

```sh
python3 scripts/p337_typed_arm_gluing.py
python3 tests/test_p337_typed_arm_gluing.py
```

Output: `results/exact-typed-arm-gluing/latest.json`.  Five focused tests use only the Python standard library and the preceding exact homology certificate.

Related issues: #337, #334.
