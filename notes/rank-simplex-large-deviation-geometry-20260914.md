# Rank-simplex large-deviation geometry: rare-sector free energies are the canonical coordinates

Date: 2026-09-14. Exact algebra plus a fixed-width long-cylinder asymptotic interpretation. This note unifies the rank-source zero analysis, the width-two spectral bridge, and the canonical `(b,d)` rank coordinates introduced in draft PR #773.

## 1. Exact finite softmax coordinates

For any rank law with positive `P0,P1,P2`, define

    b = (1/2) log(P0/P2),
    d = log[P1/sqrt(P0P2)].

Then exactly

    (P0,P1,P2)
      = (e^b, e^d, e^{-b}) / (e^b+e^d+e^{-b}).          (1)

The three equal-weight boundaries are straight lines:

    P0=P2  <=> b=0,
    P0=P1  <=> b=d,
    P2=P1  <=> b=-d.                                    (2)

Thus the balance root is only one of three distinguished lines in the full topology simplex.

A topological source `sX`, `X=r-1`, acts as

    b -> b-s,
    d -> d.                                              (3)

Hence the source never changes the intrinsic even shape coordinate d; it translates the physical point horizontally through the softmax diagram.

## 2. Fixed-width long-cylinder limit

Suppose at fixed width and fixed occupation parameter p the long-period rank-sector probabilities obey

    P0(m)=a0 exp(-m alpha0)[1+o(1)],
    P2(m)=a2 exp(-m alpha2)[1+o(1)],
    P1(m)->1,                                            (4)

with positive rates `alpha0,alpha2`. Then (1) gives directly

    b/m -> (alpha2-alpha0)/2,
    d/m -> (alpha0+alpha2)/2.                            (5)

Conversely,

    alpha0 = lim (d-b)/m,
    alpha2 = lim (d+b)/m.                                (6)

So the canonical odd/even coordinates are not merely a convenient parameterization: in a dilute two-rare-sector regime they are exactly the half-difference and half-sum of the two topology free-energy rates.

At an open/closed cylinder crossing `alpha0=alpha2=alpha`,

    b/m -> 0,
    d/m -> alpha.                                        (7)

The merged width-two bridge #705 is an explicit model-specific realization of this statement.

## 3. Extensive topological source is the same softmax picture

Let `s=m sigma`. By (3), the three exponential weights are governed by

    r=0:  b-m sigma,
    r=1:  d,
    r=2: -b+m sigma.                                    (8)

Using (5), the rank-0/rank-1 and rank-1/rank-2 boundaries are

    sigma_- = (b-d)/m -> -alpha0,
    sigma_+ = (b+d)/m -> +alpha2.                        (9)

This reproduces the previously derived extensive-source three-region law without a separate calculation:

    sigma < -alpha0      -> rank 0 selected,
    -alpha0<sigma<alpha2 -> rank 1 selected,
    sigma > alpha2       -> rank 2 selected.             (10)

The apparent direct rank-0/rank-2 equality line is `sigma=b/m`, but when both rare-sector rates are positive the rank-1 weight lies above both there. Therefore a cylinder open/closed crossing does not create a direct extensive-source 0/2 transition at zero source.

## 4. Source zeros and phase boundaries are two analytic continuations of the same coordinates

For finite m, complex source zeros are zeros of

    e^{b-s}+e^d+e^{-b+s}.                               (11)

The fixed-source zero-free strip follows from the bounded three-state variable. In the extensive coordinate `s=m sigma`, the real parts of zero families approach the real phase boundaries in (9), while their imaginary parts are O(1/m) in sigma. Thus the source-zero pinch in the changing extensive coordinate is simply the analytic continuation of the same rare-sector softmax geometry.

This is a closure/topology chemical potential. It is not a microscopic extensive field and is not evidence of a physical Jordan block.

## 5. Consequences for the broader program

1. `b` is the natural matching-odd free-energy coordinate; root consistency is the location where the physical path crosses `b=0`.
2. `d` measures the cost of leaving the intermediate rank-1 sector relative to the geometric mean of the two extreme sectors.
3. The self-normalized curve `d=D_L(b)` contains both finite-aspect rank shape and, in long cylinders, the two sector-rate functions.
4. The distinction between balance-root convergence and full-law concentration should be studied as a property of the whole physical path `(b(p),d(p))`, not only the scalar zero `b=0`.

No claim here supplies a width-uniform spectral estimate or a new full-law theorem. It only identifies the common coordinate system in which those questions should be stated.
