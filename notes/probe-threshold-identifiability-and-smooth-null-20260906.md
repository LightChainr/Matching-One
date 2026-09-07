# Threshold identifiability scope note and the smooth-flow discriminator checklist

**Scope:** two adjacent #599 programs — "does a bounded finite-horizon
task order constrain an infinite-volume threshold?" and "when does a rank-1
finite-size flow deserve an RG/state reading?"  Both notes are deliberately
short: their job is to keep a bounded-realization claim and a rank-1-flow claim
from silently becoming threshold evidence.

---

## Part 1. Finite task order does not constrain an infinite-volume threshold

### Exact anchor (already in the repository)

#435/#549 give two families of finite gadgets with *identical complete
unbranched survival polynomial* `S(z)` — hence identical unbranched survival
behaviour for **every** future length and **every** per-vertex opening
probability `p` — that nevertheless differ under a delayed branching/fork test.
In the parallel-bank form the fork probabilities differ by the positive amount

```
F_{k,a+1} - F_{k,a} = 1 / [2k (8k-1)^2]   (> 0).
```

So the *entire* unbranched future language (a "complete survival law", which
fixes all unbranched finite-horizon statistics and all finite linear/Hankel
data derived from them) cannot separate the two mechanisms, while a branched
test can.  Any statement of the form

> "a chosen connectivity task has uniformly bounded finite-horizon
> balanced/Hankel order as width grows and p varies, therefore the task
> constrains the infinite-volume threshold"

therefore needs an explicit bridge that crosses from the *declared unbranched
task language* to the *branching object that actually separates the
mechanisms*.  The repository itself has the counterexample: the bounded-order
unbranched description is compatible with both fork behaviours, and the fork
behaviour is where the difference lives.

### Generic no-go shape (for #594's comment-1 target)

#594's comment formulates the target as: "bounded balanced order + no uniform
coupling/observability assumption => no threshold identifiability".  The clean
mathematical shape is a *sector surgery*: take any finite task E (finite
sources, readouts, horizon, future depth), fix a microscopic family whose E
data is what you like, and alter only coordinates that E cannot see.  Because E
is finite, there is always a large space of invisible coordinates; the
question is whether they can move an infinite-volume singularity.  For
translation-invariant percolation-type families they plausibly can (the frozen
single-site task data does not pin the large-deviation tail that sets the
threshold).  We do **not** claim to have closed this for site percolation
itself; the exact anchor above is the finite-object certificate that the
general danger is real in the repository's own language.

### What would be required for a positive bridge (checklist, no claim)

critical-sector completeness; horizon growing with the correlation length;
uniform spectral approximation; positivity; renormalization composition;
exact topological balance; planar duality/matching; or an explicit map from a
finite realization bound to an infinite-volume singularity location.  Without
one of these, "low balanced order" and "E_p[X] = 0" remain two separate
programs; the homological-balance statement is not implied by any finite
finite-horizon order statistic we know.

---

## Part 2. Rank-1 full-law flow: what the null model says and what it cannot do

### The null is accepted

#582/#584 already adopted the point: for a smooth one-parameter family
`Q(s+h,u) = Q(s,u) + h d_s Q + h^2/2 d_s^2 Q + ...`, any *local* finite
difference is approximately rank 1, so a high one-direction explained variance
over one-step transitions is generic.  Nothing in this note re-opens that.

### Discriminator checklist (cheap, no new sampling)

For a claim that a measured finite-size flow is more than smooth tangent
geometry, report at least one of the following; a smooth family satisfies none
automatically:

1. **Step-size consistency of the leading direction**: the tangent estimated
   from scale words `(N -> m1 N -> m1 m2 N)` equals the tangent from the direct
   word `(N -> m1 m2 N)` (cocycle / composition law).  Smooth families do not
   satisfy this unless the flow is exactly the exponential of one vector field.
2. **Second-order curvature law**: the curvature `[Q_{e^{2h}N} - 2 Q_{e^h N}
   + Q_N]/h^2` equals the derivative of the first tangent with respect to the
   scale word; a rank-1 exponential flow predicts a *specific* curvature, not
   any small residual.
3. **Crossed context labels**: a frozen direction trained on one lineage of
   scale words predicts a held-out lineage from a *different* context/Smith/
   deck class (already the #584 spirit) *beyond* the Taylor order.
4. **Cancellation of the affine tangent**: the rank-1 direction must survive
   with the *weighted* metric used in #582 (`0.006 sigma` after whitening), not
   only the unweighted `||v_shape||` number.
5. **Non-generic spectrum of the remainder**: the structured remainder after
   removing the tangent must cluster by pre-declared labels (Smith/deck/lineage
   /topology) and must not be reproducible from a synthetic smooth family with
   the same tangent word.

### Status

No new computation is reported here; this note is the operational checklist
that keeps a "99% one-direction flow" from being quoted as RG evidence.  The
repository's #584 null-adoption already implements items 1–2 as required
reporting; items 3–5 are the ones that would actually upgrade #582's flow to a
compositional/RG reading, and none of them is implied by the PCA number.
