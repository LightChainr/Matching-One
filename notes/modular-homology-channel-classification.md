# Modular classification of Matching One torus channels

Status: exact Phase-A bridge for Issue 114.  This note classifies lattice
topology labels; it does not identify their CFT matrix elements.

## 1. Exact action

The torus homology backend records each wrapping component by an integer
winding subgroup `W` of `H_1(T^2,Z)`.  Changing the two period generators by
`M in SL(2,Z)` sends every winding vector to `M w`.  Because `M` is invertible
over the rationals,

```text
rank_Q(M W) = rank_Q(W).
```

The repository definitions then immediately give

```text
either = [rank > 0],
cross  = [rank = 2].
```

Therefore `rank`, `either`, and `cross` are modular scalars.  This statement
is independent of a particular primitive basis for `W` and includes
finite-index rank-two subgroups.

## 2. Which familiar names are not scalars

`direction_0`, `direction_1`, and `both` use the chosen period generators.
They are useful typed observables, but they are not modular scalars.

For an exact counterexample, start with the rank-one spiral `(1,1)`.  It has
both directional flags.  The determinant-one shear

```text
[[1,-1],
 [0, 1]]
```

sends it to `(0,1)`.  Rank and `either` are unchanged, while `direction_0`
and `both` change.  The lower shear sends `(1,1)` to `(1,0)` and similarly
changes `direction_1` and `both`.

This also explains why an individual primitive winding character is naturally
vector-valued: `SL(2,Z)` permutes/mixes the primitive sectors.  Summing all
rank-one directions restores a scalar rank channel.

## 3. Matching/complement combinations

Geometric relabelling acts on the black-primal and white-matching winding
subgroups by the same homology matrix.  Occupation complement is pointwise and
commutes with this relabelling.  Consequently any fixed sum or difference of
the scalar rank channels remains scalar:

```text
primal cross/either,
matching cross/either,
matching-even cross/either,
matching-odd cross/either.
```

This is a classification of the channel label, not a claim that every
microscopic orientation response is already a single homogeneous CFT field.

## 4. Elliptic stabilizer consequence

For a homogeneous first-order spin-`s` perturbation in a modular-scalar
channel, invariance under an order-`n` elliptic stabilizer requires

```text
s = 0 mod n.
```

At `tau=i`, `n=4`; at the hexagonal point, `n=6`.  A square-lattice harmonic
that also survives at the hexagonal point therefore has spin divisible by 12:

```text
H4:  forbidden,
H8:  forbidden,
H12: allowed.
```

The result applies directly to a typed scalar `cross` or `either` response.
It does not apply to `both`, a fixed primitive character, or an untyped mixture
of scalar and vector-valued channels.  This distinction is compatible with the
observed non-scalar primitive character in the Pell/hexagonal analysis.

## 5. What this advances—and what it does not

This completes the channel-scalar part of the requested Issue-114 bridge and
turns one premise of the elliptic H4 zero into an exact lattice-topology
statement.  It supplies a fail-closed rule: only `cross/either` rank channels
may inherit the scalar elliptic selection rule without further representation
data.

Still unresolved:

- which FK/Potts/CFT matrix element the finite matching observable computes;
- whether its response is an ordinary one-point insertion, an integrated
  homology correlator, a defect/twisted matrix element, or a logarithmic block;
- the overlap of the global matching observable with thermal Q4 versus the
  lower `V_(2,2)` four-leg spin-4 field;
- field-normalization and projector derivatives at `Q -> 1`.

The parent issue therefore remains open.
