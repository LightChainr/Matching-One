# A physical two-hole closure separates a finite tangent from a local Q1 family

**Result.** The specified C4-averaged four-port tensor has a finite,
nonzero one-insertion original-U response, but its **unrenormalized
two-insertion conditional response has a nonremovable Q1 pole**. A concrete
square-lattice exterior realizes the pole. This supplies a mechanism
decision without another N25 source score or a new stochastic block.

The completed linear result is
[`923f66b9`](https://github.com/LightChainr/Matching-One/blob/923f66b979a6b6132875f783106c041ed3c0c1a9/notes/local-four-port-transmission-result.md),
`branch_only`, V_av(25)=+0.0018155512845251097. Nothing below negates it.
The kernel is fixed by
[`9dc3c426`](https://github.com/LightChainr/Matching-One/blob/9dc3c4269c1c44a0a0a82b15f8cd13f922a8b2d4/notes/closed-source-local-four-port-pair-kernel.md).
In particular this result addresses precisely the higher-insertion
regularity left open in that kernel note, rather than redefining its
single-insertion claim.

## 1. One fixed two-copy contraction

Let P(a,b;c,d)=i P_[Q−2,2] i† and R(a,b,c,d)=P(b,c;d,a). Keep the
specified Kbar=(P+R)/2. All statements about ordinary colour spaces start
at integer Q≥4; continuation is the resulting rational diagram function.
The known projector identity gives ||P||²=||R||²=d2=Q(Q−3)/2. Direct
contraction of the 15 exact equality patterns of four colour indices gives

```text
<P,R> = Q(Q−3)(Q²−3Q+4)/[4(Q−2)(Q−1)],
G(Q)=||Kbar||²
    = Q(Q−3)(3Q²−9Q+8)/[8(Q−2)(Q−1)].                 (1)
```

Each equality pattern with b distinct colours has (Q)_b assignments;
there is no interpolation from a finite set of integer evaluations.
The [recorded symbolic calculation](../results/local-pair-two-insertion/latest.json)
retains all 15 terms. It gives

```text
G(Q) = 1/[2(Q−1)] + O(Q−1).                           (2)
```

The finite part is zero. Double poles of the individual sector terms
cancel, but this simple pole remains. The independent
[fixed-cut sector resolution](local-pair-crossing-sector-resolution.md)+explains the representation mixing behind (1); Kbar is not a single
fixed-cut projector after its spatial average.

## 2. A realizable embedded exterior, not an abstract colour wiring

Use the 17×17 square torus and two vacant holes x=(3,4), y=(11,4).
Occupy exactly the following four NN paths; all other sites remain vacant:

| Path | Occupied coordinates | Hole ports connected |
|---|---|---|
| Direct | (j,4), 4≤j≤10 | xE to yW |
| Upper | (3,5), (11,5), and (j,6), 3≤j≤11 | xN to yN |
| Lower | (3,3), (11,3), and (j,2), 3≤j≤11 | xS to yS |
| Outer | (2,4), (12,4), (1,j) and (13,j), 4≤j≤8, and (j,8), 1≤j≤13 | xW to yE |

Their lengths in vertices are 7,11,11,23. They have no mutual NN
contacts: there is at least one vacant lattice line between neighbouring
paths. Thus K=52, B=48, C_B=4, B_vac=418 and c_H=422. The full graph,
including either or both filled holes, is contained in the injective
contractible rectangle [1,13]×[2,8]. Its original ambient rank stays zero.

The paths identify x's colours (N,E,S,W) with y's (N,W,S,E). This is a
reflection, and Kbar is D4-invariant: that reflection exchanges P and R.
The two-insertion colour contraction is consequently exactly G(Q), with
418 spectator blocks multiplying it by Q^418. No artificial spectator
connection, wrap, colour alias or rank modification is involved.

At one hole alone, the other hole vacant leaves four singleton exterior
blocks; the other hole occupied makes all four equal. Both Bell4
contractions of Kbar vanish. A single insertion is therefore zero in
this entire fixed exterior, despite the nonzero two-insertion closure.

## 3. The pole survives connected normalization

Fix those exterior occupations, but sum over both holes' vacant and
occupied choices, with activities v_x,v_y>0 and independent local
couplings epsilon_x,epsilon_y. Remove the common exterior activity and
Q^418. The exact conditional partition is

```text
Z_xy(Q) = Q^4+(v_x+v_y+v_x v_y)Q
          +epsilon_x epsilon_y G(Q).                  (3)
```

All terms linear in either epsilon vanish, including those with the
other hole occupied. The connected mixed derivative at zero coupling is

```text
partial_epsilon_x partial_epsilon_y log Z_xy
  =G(Q)/[Q^4+(v_x+v_y+v_x v_y)Q],
Res_(Q=1) = 1/[2(1+v_x)(1+v_y)] > 0.                  (4)
```

A common Q-dependent partition prefactor cannot change (4). Separate
analytic single-site quadratic counterterms O(epsilon_x²) and
O(epsilon_y²) cannot change its mixed coefficient either. With uniform
site-average units, the named unordered pair contributes G/N² to the
epsilon² partition coefficient; differentiating twice adds the usual
factor two. None of these unit conversions removes the pole.

This excludes a regular finite-strength Q1 tensor family **with this
unmodified linear vertex and all physical conditional exteriors**. It
does not by itself prove that a particular fully summed homogeneous
partition or its original-U derivative diverges: summation over other
exteriors can cancel diagram residues. The counterexample already
suffices for the stronger local-family assertion, which must work in
every allowed fixed exterior.

## 4. Three objects that must now remain distinct

1. **The bounded Q1 occupation tangent.** S_av=−sum_x t_x/N remains
   well-defined and has the measured finite N25 response. One may define
   a positive occupation reweighting exp(epsilon S_av); its first
   derivative agrees with the already closed one-insertion source.
2. **The original multi-insertion local colour tensor.** Its coefficient
   at two distinct sites is the joint colour contraction, not the product
   of individually closed marks. In the exterior above both t_x and t_y
   vanish for every unperturbed choice of hole occupations, while (4)
   has a pole. Thus Cov(t_x,t_y) cannot be silently substituted for this
   two-insertion tensor response.
3. **A renormalized or confluent local Q1 field.** It would need an
   explicitly different continuation/cancellation. For example,
   Ktilde=sqrt(Q−1) Kbar has two-copy limit 1/2 in (1), but every finite
   single-insertion contraction, including the measured original-U
   response, then tends to zero. This normalization changes the
   mechanism; it does not preserve the old nonzero linear coefficient.

The next theoretical question is consequently specific: **can a stated
finite combination cancel the separated-insertion pole while retaining
the intended nonzero single-insertion original-U response?** A formal
pair-representation name or a common normalization is not such a
combination. If no such completion is proposed, the well-defined
occupation-tangent size comparison remains available on its own terms.

The [fixed-size prediction](local-pair-size-response-predictions.md)
uses W_N=N V_av(N) and R=W_(4N)/W_N. Under its declared single-field,
nonzero-loading assumptions, x=17/4 predicts R→2 and x=21/4 predicts R→1.
These linear-response alternatives do not require pretending that the
unrenormalized two-insertion field is regular. An order-one W can also
describe modulation of pre-existing anisotropy, so it is not a unique
thermal-Q4 identification.

## Delivery

[Definition](../analysis/local_pair_two_insertion_contract.json),
[script](../scripts/analyze_local_pair_two_insertion.py),
[result](../results/local-pair-two-insertion/latest.json) and
[receipt](../results/local-pair-two-insertion/run.json) preserve the
specific contraction and witness. The recorded symbolic calculation took
about0.114seconds in the managed research Python environment. It was
preceded by a direct symbolic exploration, as the receipt states; it is
an exact derivation, not a prospective statistical test. No occupation
enumeration, new random samples, root search or cloud job was performed.
