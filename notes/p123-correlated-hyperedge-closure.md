# Correlated four-terminal closure and the square-site obstruction

Damavandi and Ziff's four-edge coordinates give the correct finite local language. For an isotropic planar four-terminal cell, the six orbit probabilities have multiplicities `(1,4,4,1,2,2)` and duality exchanges

```text
(P1,P2,P3,P4,P5,P6) -> (P4,P3,P2,P1,P6,P5).
```

Thus local selfduality is exactly

```text
P1=P4,  P2=P3,  P5=P6.
```

The Bernoulli site embeds without approximation as the all-or-none tensor

```text
A(t)=(1-t,0,0,t,0,0).
```

It intersects the local selfdual manifold only at `t=1/2`, transversely since `d_t(P4-P1)=2`. Correlated hyperedges therefore solve the old independent-bond support problem, but do not turn the homogeneous square-site threshold into a new local selfdual number.

## A minimal composition obstruction

Take two all-or-none four-hyperedges sharing two internal vertices. In one orientation they expose the boundary pairs `AB` and `CD`; average equally with the rotated `BC`/`DA` orientation to restore C4. With independent activation probability `t`, exact enumeration gives

```text
P1=(1-t)^2,
P2=t(1-t)/2,
P3=0,
P4=t^2,
P5=P6=0.
```

The duality-odd coordinates are

```text
P4-P1=2t-1,
P3-P2=-t(1-t)/2.
```

They have no common root: the first demands `t=1/2`, while the second demands `t=0` or `1`. So even the first nontrivial local composition of exact site tensors immediately creates an unpaired partial-connectivity direction.

This remains true for the whole correlated two-block family. Let, within either orientation,

```text
Prob(00)=a,
Prob(10)=Prob(01)=b,
Prob(11)=c,
a+2b+c=1.
```

Then `P=(a,b/2,0,c,0,0)`. Selfduality forces `b=0` and `a=c=1/2`: the only solution is perfect 00/11 correlation, which collapses the composite back to one all-or-none half-occupied hyperedge. This is a finite-local algebraic no-go, not a numerical observation.

## The positive correlated manifold

Allow every partition together with its dual partner. The maximal isotropic probability simplex fixed by local duality is

```text
P=(w,u,u,w,v,v),
w=1/2-4u-2v,
u,v>=0, 4u+2v<=1/2.
```

This is a genuine two-dimensional correlated-cell selfdual family. It exposes the price of closure: partial pair/triple states are independent couplings, and `P4=w<=1/2`. Moving away from the all-or-none point moves the all-connected mass downward.

There is also an exact inhomogeneous escape. Put any tensor `P` on one hyperedge colour and `dual(P)` on the other. Lattice duality exchanges the two colours, so the alternating model is globally selfdual. For all-or-none cells this alternates `A(t)` with `A(1-t)` for arbitrary `t`, precisely the type of construction described in the paper. It becomes homogeneous only at `t=1/2`.

Hence four-edge selfduality gives a structural trichotomy for square-site:

1. retain the exact homogeneous Bernoulli tensor and be fixed at `1/2`;
2. compose locally and generate partial-connectivity couplings;
3. keep arbitrary `t` only by alternating it with its dual `1-t`.

None supplies a new exact equation for homogeneous square-site `p_c`. A surviving exact route must use a larger/nonlocal projection in which `P4` is no longer literally the site-occupation marginal, or find extra integrable structure beyond local four-edge duality.

Source used for coordinate conventions and the alternating all-or-none construction: [Damavandi and Ziff, *Percolation on Hypergraphs with Four-Edges*, arXiv:1506.06125](https://arxiv.org/abs/1506.06125).
