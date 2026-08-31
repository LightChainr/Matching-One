# A concrete local four-port perturbation and its original-U response

The completed seam-trace response leaves a specific microscopic question:
what does a **local** pair tensor do in the original observable? The
[fixed contract](../analysis/p337_local_pair_insertion_contract.json)
selects one tensor, one normalization and the old N25 geometries. No
equivalence to the full seam trace is presumed.

The [explicit rational kernel and all Bell4 contractions](closed-source-local-four-port-pair-kernel.md)
and [original-lattice topology proof](local-four-port-pair-insertion.md)
establish the microscopic source below before its numerical response is
evaluated.

## 1. The perturbation is an actual vertex tensor

The hypergraph colour representation has a colour on each original NN
edge and the four-port site tensor

```text
V_x(N,E;S,W)=1+v delta_(N=E=S=W).
```

On two ordered colour ports put `Pi2=i P_[Q-2,2] i^dagger`, where i embeds
an unequal unordered pair as its normalized symmetric ordered pair.
Initially group the incoming ports as N,E and the outgoing ports as S,W.
Average this four-index tensor over the four C4 rotations, giving Kbar.
The perturbation is `V_x -> V_x+(epsilon/N) Kbar` at every site, to first
order in epsilon. The added term belongs to the vacant-site summand.

The rank factor and q/E remain those of the original occupation A with
that site vacant. Projector diagrams are colour contractions, not extra
occupied NN bonds or a newly occupied site.

The unrotated Pi2 is orthogonal to both available bare vertex tensors:
the constant tensor and the all-equal tensor are annihilated on each
side. Its insertion is therefore not an activity derivative of that
fixed-Q bare vertex. C4 averaging preserves their complete scalar
contraction zeros, but recouples different input/output groupings. Kbar
need not itself be a projector or remain in a pure Pi2 block of one
chosen cut. This distinction matters for any eventual field assignment.

## 2. The contracted local source and its site-average unit

With x vacant, partition its four incident edge-nodes according to their
components in the exterior hypergraph. A singleton colour index has zero
partial sum in Pi2; identifying either incoming or outgoing pair also
annihilates it. The only nonzero external partitions for the unrotated
tensor are `(N,S)|(E,W)` and `(N,W)|(E,S)`. Each gives
`d2(Q) Q^(cH-2)`, where `d2(Q)=Q(Q-3)/2` and cH is the unmodified
hypergraph component count. The two lines can have equal colours in the
original law; equality of colours is not equality of components.

After C4 averaging, the relative occupation weight is

```text
beta_x(Q)=(Q-3)/(2Q) t_x,
t_x = I_((N,S)|(E,W))
      + (1/2) [I_((N,E)|(S,W))+I_((N,W)|(E,S))].
```

These events require x vacant, all four NN neighbours occupied, and
exactly two distinct occupied NN components, each containing two of those
neighbours. All other configurations have t_x=0. At Q1 the specified
closed occupation continuation is `beta_x(1)=-t_x`; it is bounded and
gives a positive perturbed occupation law for epsilon in a neighbourhood
of zero. Its full finite closed contraction, not the singular projector
entries evaluated colourwise at Q1, defines this continuation.

Write `S(A)=-(1/N) sum_x t_x(A)`. For any translation-invariant O(A),
homogeneity of the original finite occupation family gives exactly

```text
<O S> = -<O t_origin>.
```

This applies to O=1,q,E and to every thermal derivative of these finite
polynomials. Therefore a fixed-origin traversal supplies the entire
first-source global-U response to the **site-average** perturbation.
There is no missing factor of N. An extensive source would multiply the
answer by N. Multi-insertion responses would require cross-site moments
and cannot use this replacement for S squared.

## 3. Its relation to the closed source and a birth transition

On the support of t_x, activating x joins exactly two existing hypergraph
components. It changes K by1 and cH by-1. For the already defined source

```text
Sstar=2 cH-r+3K-2N+1,
ell=r(A union {x})-r(A),
```

the exact gain is consequently

```text
Delta_x Sstar=1-ell.                                       (1)
```

For the opposite pairing NS|EW, the two cycles closed on adding x
have algebraic intersection plus or minus1. Hence the post-addition
rank is2, while the pre-addition rank is0 or1. This includes a genuine
rank0-to-rank2 double birth at one site. Adjacent pairings have examples
at every original rank. These facts connect the local colour contraction
to a named rank-transition event without replacing its pre-birth q/E.

Thus this local two-cluster routing event is attached to a concrete birth
transition, not an arbitrary new contact descriptor. Its weight does not
equal the thermal increment ell: it distinguishes the local pair routing
and retains the pre-activation rank in q/E. In particular the occurrence
of this event alone does not justify replacing the source by a rank
current or by Sstar. The corresponding observable gains, if needed for
the transition itself, are `Delta q=ell` and
`Delta E=ell(2r-2+ell)`.

## 4. Complete original-observer transmission

For each geometry form separately normalized expectations and the
covariance response `jO=Cov(O,S)`. At the old pooled matching root set

```text
M=(q_axis+q_tilt)/2,       Y=(E_axis-E_tilt)/Delta,
Delta=1152/625,           D=M_h,
R=Y_h/D,                 A25=25^(13/8)/2.
```

The old root, D and R are imported, with `h=p/(1-p)` and
`D_h=D_p/(1+h)^2`. The exact first response is

```text
V_S/A25 = jY_h/D - Y_hh*jM/D^2
          - R*jM_h/D + R*M_hh*jM/D^2,
h_epsilon = -jM/D.                                       (2)
```

The source covariances are evaluated as
`j_q=<qS>-<q><S>` and
`j_q,h=<qS>_h-<q>_h<S>-<q><S>_h`, and similarly for E.
Unlike the completed seam trace, this local insertion is not restricted
to rank1. Its direct q/E numerators must be retained. Equation (2) is
the existing full observable functional, not a new estimator definition.

The only frozen numerical question is whether this specified V_S is zero
on the original N25 pair. A nonzero value establishes one concrete local
tensor-to-global-U route; it would not turn Kbar into the full central
seam projector or a named continuum primary. The kernel and topology
proofs are complete; the missing moments have been released for one
fixed enumeration under the contract.
