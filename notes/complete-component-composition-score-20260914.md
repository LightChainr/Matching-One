# Complete-component composition score: rewrite the K/B Palm Ward identities as a Bernoulli-like composition law

Date: 2026-09-14. Algebraic reparameterization plus scaling conjectures for draft PR #773 / issue #774.

For a complete winding component C in the direct two-fugacity activity, let

    K = number of occupied component sites,
    B = number of distinct external vacant boundary sites,
    T = K+B,
    p_hat_C = K/(K+B).                                    (1)

The physical activity is

    p^K q^B, q=1-p.                                       (2)

The fixed-width Palm score already proved in PR #739 is

    S = K/p - B/q.                                        (3)

## 1. Exact composition form of the score

Substituting `K=T p_hat_C` and `B=T(1-p_hat_C)` gives the exact identity

    S = T(p_hat_C-p)/(pq).                                (4)

Equivalently,

    qK-pB = T(p_hat_C-p) = pq S.                          (5)

Thus the thermal derivative of the winding intensity measures a **size-weighted composition bias** of the complete component plus its forced vacant boundary.

Each fixed shape has its activity maximized, as a function of p, at

    p_shape = K/(K+B)=p_hat_C.                             (6)

The Palm Ward cancellation near criticality can therefore be read as a saddle-point statement: the shapes carrying most activity have empirical occupied fraction close to the ambient p when occupied sites and distinct forced-vacant boundary sites are combined.

## 2. First Ward identity in composition variables

The exact score equation becomes

    pq d(log nu_w)/dp
      = E_Palm[T(p_hat_C-p)].                              (7)

Under the conditional near-critical intensity scaling

    nu_w(p_c+lambda/(a_t w^(3/4)))
      = w^-1 I(lambda)[1+o(1)],                            (8)

its right side is only `O(w^(3/4))`.

If `E T` is of order `w^(91/48)`, the **T-weighted mean composition bias** is therefore

    E[T(p_hat_C-p)]/E T
      = O(w^(-55/48)).                                    (9)

This is the invariant content behind the earlier ratio prediction `E B/E K -> q/p`.

## 3. Second Ward identity suggests binomial-scale transverse fluctuations

The exact second derivative is

    (log nu)''
      = Var(S)-E[K]/p^2-E[B]/q^2.                         (10)

Multiply by `p^2q^2` and use (5):

    Var[T(p_hat_C-p)]
      = q^2 E K + p^2 E B
        + p^2q^2 (log nu)''.                              (11)

This is exact at every fixed width.

If composition locking makes `K~pT`, `B~qT` at leading order, then

    q^2 E K+p^2 E B
      ~ pq E T.                                           (12)

The near-critical scaling makes the final term in (11) at most `O(w^(3/2))`, lower than `E T~w^(91/48)` under the stated fractal-mass hypothesis. Hence the stronger Ward prediction is

    Var[T(p_hat_C-p)]
      ~ pq E T.                                           (13)

This has exactly the variance scale of a Bernoulli composition score from an effective sample of size T.

## 4. Standardized composition variable and a stronger CLT conjecture

Define

    Z_C = [K-p(K+B)] / sqrt[pq(K+B)]
        = [qK-pB] / sqrt[pq T].                            (14)

Equation (13) suggests `E Z_C^2` should be of order one, but because T is random it does not by itself imply a unit second moment.

A deliberately stronger, falsifiable conjecture is:

> Conditional on the macroscopic complete-component shape/size amplitude, the transverse occupied-vs-boundary composition obeys a central-limit law,
>
>     Z_C => Normal(0,1),                                  (15)
>
> while the leading macroscopic randomness of `(K,B)` lies along the composition ray `(p,q)`.

This would refine the strong random-ray proposal

    w^(-91/48)(K,B) => A(p_c,q_c)                          (16)

by adding `sqrt(T)` transverse fluctuations.

There is no proof of (15). Distinct boundary sites are geometrically dependent and are not an independent Bernoulli sample. The claim is motivated only by the exact exponential-family score structure and the proposed asymptotic scale separation.

## 5. Cumulant diagnostics available from the same two-fugacity resolvent

Because the direct complete-component generating function already carries fugacities `(u,v)` for `(K,B)`, directional derivatives along the physical thermal score yield cumulants of S. Equivalently, derivatives in coordinates

    total-size direction:       log u + log v,
    composition-score direction: q log u - p log v        (17)

can separate longitudinal and transverse component fluctuations.

Issue #774 should, if the existing automatic/analytic differentiation makes it inexpensive, report third and fourth cumulants of

    W = qK-pB                                             (18)

and standardized skewness/kurtosis relative to `pq E T` or `Var W`. A Gaussian trend would support (15); persistent non-Gaussianity would kill the strong composition-CLT while leaving the first two Ward identities intact.

No new state variables are required: W is a linear combination of K and B.

## 6. Geometric interpretation and limitation

B counts **distinct external vacant boundary sites**, not boundary edges and not only the accessible outer hull. This distinction is essential. The conjectured composition law concerns the exact activity coordinates used by the complete-component partition function; replacing B by another perimeter observable changes the score.

The ratio `p_hat_C` is not a local occupation density inside an arbitrary Euclidean box. It is the occupied fraction of the forced set `C union boundary(C)`. Equation (6) is therefore an activity/saddle-point fact, not a statement that the critical cluster is spatially dense with density p.
