# A delayed second-birth source survives the moving-root E_top projection

2026-09-30. Finite author derivation with an explicit square-L4 positive
control. This supplies a complete **dynamic source -> second birth ->
M/E_top -> moving-root normal response** map. It is not identification of
the original norm-4 source or the separately normalized six-coordinate U.

## 1. What is unchanged: first birth and the occupation clock

Apply the [normalized safe-degree kick](birth-selection-intervention-20260930.md)
only at a fixed insertion count b and only when the current rank is one.
Before b the source is uniform; after this one kick it is uniform again.
Total completing probability c/m and total arrival rate remain unchanged.

Couple all sources with the same arrival times and the same prefix through b.
If a prefix is affected, its first birth has already happened. If its rank
is zero or two at b, it is not affected at all. Therefore its first-birth
count J1 is unchanged pathwise, and so is its label time T1 under the common
clock. The entire occupation-count law K(p) is also unchanged. The second
birth J2, and hence the joint birth law, can change.

At any later observation p, write P_r(theta,p) for rank probabilities and
g(p)=partial_theta P_2(0,p). Then exactly

    partial_theta P_0 = 0,
    partial_theta P_1 = -g,
    partial_theta P_2 = g,
    partial_theta M = partial_theta E_top = g,              (1)

where M=P2-P0 and E_top=P2+P0. These are **source derivatives**, not the
ordinary p derivatives. The latter generally move both birth marginals.

A single-kick response conditional on rank_b=1 must be multiplied by the
unperturbed rank_b=1 probability before making an unconditional rank claim.
The current L512735-step calculation uses a fixed count endpoint; it is not
yet the binomial count mixture that yields a numerical label-time g(p).

## 2. Moving the balance root does not erase a nonzero response

At theta=0 the finite model is ordinary iid occupied-site percolation.
Let p_* be its interior balance root, M(0,p_*)=0. Differentiating the moved
root p_*(theta) gives

    p_*'(0) = -g(p_*)/M_p,
    d/dtheta E_top(theta,p_*(theta)) |0
        = g [1 - (E_top)_p/M_p]
        = g * [-2(P0)_p / M_p].                            (2)

All p derivatives are at the baseline root. In a nondegenerate finite
torus, (P0)_p<0 and (P2)_p>0 for every interior p: each is a nonconstant
monotone Boolean event, and differentiating the product measure gives a
sum of pivotal probabilities. At least one pivotal configuration exists;
every configuration has positive weight for 0<p<1. Thus M_p>0 and the
factor in (2) is strictly between0and2.

**A nonzero one-sided second-birth source therefore has a nonzero
moving-root E_top response of the same sign.** It cannot be entirely a
thermal coordinate shift: a nonzero p shift would change P0, whereas this
source leaves P0 unchanged. This is an exact two-coordinate statement,
not a classification of a local CFT field or a result for another normalizer.
If g(p_*)=0 by cancellation, the conclusion does not assert a nonzero signal.

Equivalently the actual two-parameter response Jacobian satisfies

    det d(M,E_top)/d(p,theta) = -2(P0)_p*g.                 (2a)

When g is nonzero it has rank two: temperature and this second-birth source
are locally independent response directions. This uses an explicitly
normalized transition source, not extra derivative-jet coordinates or a
rank increase created only by relabelling the same scalar source.

## 3. An explicit symmetric finite source with g(p)>0 everywhere inside

Use the already computed square-L4 preparation

    A={0,1,4,5,8,12},  v=x+4y.

At b=6, activate the source precisely when the occupied set belongs to
A's orbit under translations and the square dihedral group. There are64
distinct such sets, all rank one, c=0, with the same continuation law by
graph symmetry. At theta=0 their total prefix probability is

    alpha = 64/binom(16,6) = 8/1001.

This source is invariant under lattice translations and dihedral symmetries;
it is deliberately constructed and history-dependent, not the physical
norm-4 Q-source. The previous exact full-lag calculation for A gives negative
single-kick survival derivatives at h=2,...,6 and zero at the other h.
It is reused here, without a new lattice or permutation enumeration.

The order statistics of the common label clock are independent of the
site-choice process. With K(p)~Binomial(16,p), setting q=1-p yields the exact
unconditional second-birth response

    g(p) = alpha sum_k binom(16,k) p^k q^(16-k) [-S'_{k-6}(0|A)]
         = (96/35) p^8 q^8 + (32/9) p^9 q^7
           + (1088/525) p^10 q^6 + (256/385) p^11 q^5
           + (16/165) p^12 q^4.                            (3)

Only k=8,...,12 contribute; for k<=6 the kick has not acted, and the
one-insertion rank response is zero. Every displayed coefficient is
positive, proving g(p)>0 for every 0<p<1, including the balance root.
Together (2)--(3) give a strict nonzero normal E_top response in the actual
unconditional finite percolation ensemble at the baseline source.

The old L4 static rank table supplies the baseline polynomial. One ordinary
floating-point root evaluation, only to show scale, gives

| Quantity | Value |
|---|---:|
| baseline p_* | 0.5906721123 |
| baseline E_top | 0.6454881965 |
| g(p_*)=delta M=delta E_top | 0.00017057649 |
| moving-root factor -2(P0)_p/M_p | 0.98801982 |
| root derivative | -0.00003425820 |
| normal E_top derivative at moving root | **0.00016853296** |

The sign theorem relies on exact positive coefficients and strict
monotonicity, not on those decimal digits. Symmetrizing the control changes
its activation probability, not its local-to-global response law.

## 4. Consequence, without changing the original-source question

This control fills a conceptual gap: a source can be invisible to current
rank and to the whole occupation clock, then acquire a genuine normal
topological response through subsequent geometric growth. “The immediate
global observable is unchanged” is therefore not a general obstruction
to temporal transmission. Here the source acts on transitions, not merely
by reweighting already conditioned static configurations.

The source is intentionally selective (activation probability8/1001 at L4).
Its existence does **not** attribute natural square-lattice anisotropy,
norm-4 residuals or the old within-K source to this mechanism. The original
U still requires its actual source, six coordinates, normalizer and
moving-root forward map. A finite constructed normal response is a positive
control for that reasoning, not a substitute for those missing columns.

For the ongoing L512 source-window job, first obtain its actual735-step
transmission. Its single-count result cannot be put into (2) by replacing
p with k/N; doing so would silently change clocks. The general identities
(1)--(2) already hold, but a numerical label-time consumer needs the proper
count mixture. No extra size, p grid or high-precision root task is launched.

[Exact coefficients and illustration](../analysis/one-sided-birth-source-20260930/RESULT.md),
[result](../analysis/one-sided-birth-source-20260930/result.json),
[short table consumer](../analysis/one-sided-birth-source-20260930/derive.py).
