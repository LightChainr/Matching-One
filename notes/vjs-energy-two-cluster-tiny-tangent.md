# Exact VJS energy/two-cluster tangent on the tiny FK torus

## The finite object

Vasseur--Jacobsen--Saleur define two neighbouring marked spins at each of two
locations and the FK probabilities `P0,P1,P2`: four distinct clusters, one
propagating cluster, or two propagating clusters. Let `P_ne` be the one-pair
probability of belonging to distinct FK clusters and

```text
G(Q)=P0(Q)+P1(Q)-P_ne(Q)^2.
```

After centering their field `psi_ab`, equations (13)--(15) of
arXiv:1206.2312 reduce in the unordered-pair diagram basis to

```text
C(Q) = 4 Q^-4 J G(Q)
     + 2 Q^-2 [P_singlet(Q)+P_[2](Q)] P2(Q).
```

This formula is exactly the bridge between #258 and #262. The probabilities
carry the FK measure derivative, while #262 supplies the only finite projector
combination:

```text
P_regular(1)=I+X-4J,
dQ P_regular(1)=X.
```

The singular derivatives of `P_singlet` and `P_[2]` are never used.

## Exact L=2 result

On the square-bond L=2 torus, take the two horizontal neighbour pairs
`(0,1)` and `(2,3)`. Exhaustive enumeration of all 256 bond configurations at
Q=1 gives

```text
P_ne = 37/256,
P0   =  1/256,
P1   =  3/128,
P2   =  9/256,
P0+P1+P2 = P(both local pairs distinct) = 1/16.
```

Along the self-dual critical manifold the score is `T=k+b/2`. Every
probability derivative is therefore an exact covariance with T. Applying the
product rule to `C(Q)` yields three separately nonzero tensor rows:

```text
measure score, holding the Q=1 field fixed
  = (39/32768) I + (39/32768) X + (40617/2097152) J

finite confluent projector derivative
  = (9/128) X

explicit insertion/colour normalization derivative
  = (-9/64) I + (-9/64) X + (1881/4096) J.
```

Their exact sum is

```text
(-4569/32768) I
+ (-2265/32768) X
+ (1003689/2097152) J,
```

which independently equals direct differentiation of the generic-Q formula.
This is the first tiny finite-graph closure in the repository where the
measure score, finite projector tangent and explicit field-definition term are
all present and none is silently set to zero.

## What remains external

This exact identity does not measure an LCFT logarithm. A one-distance L=2
graph cannot recover the VJS continuum inputs

```text
lim_Q->1 (Delta_psi-Delta_epsilon)/(Q-1)=sqrt(3)/pi,
A(1)=A_tilde(1),
logarithmic coefficient 2sqrt(3)/pi.
```

Those are scaling-limit statements. The tiny oracle instead calibrates the
finite lattice derivative bookkeeping that any larger-distance realization
must reproduce.

## Reproduction

```bash
python3 scripts/exact_vjs_collision_tangent.py \
  --output results/vjs-collision-tangent/latest.json
```
