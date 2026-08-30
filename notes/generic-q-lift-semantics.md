# Generic-Q lift semantics and CP-horizontal tangent transport

Status: Issue #333 Phase B.  This is a small claim-bearing type and one exact
transport map, not a repository-wide scorer migration.

## 1. Endpoint semantics and lift semantics are a product type

Issue #146 types the observable that exists at `Q=1`: topology channel,
combination, probability coordinate, orientation order and normalization.
That remains necessary, but it cannot type a `Q` derivative by itself.

A generic-Q tangent descriptor additionally records:

```text
endpoint observable id
generic-Q lift id
sector weights as functions of Q
normalization / denominator
path Q -> (Q,v(Q))
explicit Q dependence of the lattice insertion
projector convention and its Q derivative
field normalization / counterterm convention
```

The implementation keeps this as a separate descriptor that composes with
the endpoint descriptor.  It does not enlarge the finite wrapping-channel
registry or pretend that every existing scorer differentiates in `Q`.

## 2. Registered lifts

For the common endpoint

```text
A_top(Q=1) = P_2D-P_0D,
```

the registry contains two lifts:

| lift | weights `(2D,1D,0D)` | explicit insertion derivative at Q=1 |
|---|---|---|
| `L_hom` | `(1,0,-1)` | `(0,0,0)` |
| `L_CP` | `(1,0,-Q)` | `(0,0,-1)` |

Both normalized probabilities and unnormalized restricted state sums are
registered.  Each is registered separately on fixed `v=1` and on the square-
bond critical path `v=sqrt(Q)`.  The path descriptor includes `dv/dQ`; it is
not a free-form annotation.

`L_hom` is the literal unweighted ambient-homology measurement.  `L_CP` is
the critical-polynomial / periodic-TL projector.  Their projector conventions
are not aliases even though their endpoint weights coincide.

## 3. Exact transport

At a fixed path and normalization, endpoint equality cancels all measure,
path and denominator derivatives in the *difference* of the two tangents.
Only the explicit insertion-weight derivative remains.  Therefore

```text
normalized:   D_Q L_CP = D_Q L_hom - pi_0D,
unnormalized: D_Q L_CP = D_Q L_hom - W_0D.
```

The transport routine derives this from the registered weight derivatives;
`-pi_0D` is not a special fitted constant.  Reverse transport adds the same
endpoint sector coordinate.

The Phase-A exact tori supply two checks on both declared paths:

```text
L=2: pi_0D=69/256
L=3: pi_0D=18865/65536.
```

For example, at fixed `v=1` on `L=2`,

```text
D_Q h = -39/256
D_Q c = -27/64
transport: -39/256 - 69/256 = -27/64.
```

On the exact critical path,

```text
D_Q h = pi_0D
D_Q c = 0,
```

so CP-horizontal transport removes the complete raw homology tangent.  This
does not say that every physical logarithmic tangent vanishes.  It says this
particular raw term is exactly the lift transition coefficient.

## 4. Comparison rule for #258, #262, #263 and #275

The machine gate is deliberately strict:

> Two raw Q tangents are directly comparable only when their complete
> descriptors are identical.  Equality of their Q=1 endpoint is insufficient.

For the registered homology/critical-polynomial pair, transport both objects
to the CP-horizontal section and only then compare.  Consequently:

- #275's unweighted homology tangent carries `L_hom`; subtract `pi_0D` before
  comparing it with a CP-horizontal quantity;
- #258 must state whether its measure score differentiates only the random-
  cluster measure or also an explicitly Q-dependent insertion, and must state
  the `(Q,v)` path;
- #262's representation/partition-algebra projector is a prospective `L_rep`,
  not silently `L_CP`; no map exists until the full projector derivative and
  counterterm are supplied;
- #263's continuum tangent can serve as a common comparison target only after
  its field normalization and connection section are declared.

Thus the current exact transport licenses comparison of `L_hom` and `L_CP`.
It does not manufacture transports for categorical or microscopic lifts.

## 5. Scientific boundary

CP-horizontal is a concrete comparison convention because the critical-
polynomial projector is exactly zero on the finite critical relation.  This
branch does not prove it is the unique physical LCFT connection.  A complete
VJS-type collision can add representation-projector and field-mixing
counterterms; those must become descriptor fields or a new exact transport,
not be absorbed into the existing `-pi_0D` map.

Artifacts:

- `scripts/generic_q_lift_semantics.py`;
- `analysis/generic_q_lift_transport_manifest.json`;
- `tests/test_generic_q_lift_semantics.py`.

