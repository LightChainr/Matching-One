# A bounded homology source cannot have a Lee-Yang pinch at finite real source

Date: 2026-09-13. Completed finite theorem and fixed-width crossover analysis.
The source s couples to the bounded ambient rank X=r-1, not to total occupation K,
cluster count, or a generic-Potts Q. These couplings must not be interchanged.

## 1. Uniform zero-free strip

For ANY probability law on {-1,0,1}, including every finite Matching-One ensemble,

    Z(s)=P0 exp(-s)+P1+P2 exp(s).

If s=x+iy then

    Re Z(s) = P1 + cos(y)[P0 exp(-x)+P2 exp(x)] > 0

for |y|<pi/2. This holds without a geometry, an independence assumption, or criticality.
Thus all these partition functions are nonzero in ONE strip around the entire real
source axis, uniformly in system size and occupation probability.

For real s, e^{-|s|} <= Z(s) <= e^{|s|}. Consequently the intensive log partition
function (1/N)log Z(s) converges uniformly to zero on bounded real source intervals.
No sequence of unscaled rank-source zeros can pinch a finite real s. Along any
subsequence of the three probability coefficients there is an entire limiting
Z with the same positive-real property. There is no hidden accumulation of an
unbounded number of source zeros: u*Z has degree two in u=e^s.

This is NOT a Lee-Yang circle theorem. It is a simpler bounded-variable fact. For
a random variable in [-B,B], the same cosine proof gives |Im s|<pi/(2B).
It does NOT exclude a singular limit as a function of p, and does NOT apply to
an extensive occupation source or cluster-count fugacity whose range grows with N.
It also does not imply that finite derivatives of a root with respect to s vanish.

## 2. The complete zeros are fixed by two scalar probabilities

When all P_j>0, u*Z is

    P2 u^2+P1 u+P0.

Set s_b=(1/2)log(P0/P2), u=e^{s_b}v, and

    kappa=P1/(2 sqrt(P0 P2)).

The centered quadratic is v^2+2 kappa v+1.

- 0<kappa<1: v=e^{+/- i arccos(-kappa)}, so both source zeros have Re s=s_b.
- kappa=1: v=-1 is a double zero.
- kappa>1: v=-exp(+/- arcosh kappa), hence

      s = s_b +/- arcosh kappa + (2j+1)i*pi.

At a matching balance root, P0=P2=E/2, so s_b=0. The two centered u zeros lie on
the unit circle iff E>=1/2. Rank balance alone does not impose this condition.
For E<1/2 the u zeros are distinct reciprocal negative real numbers instead.

A double source zero is not an operator Jordan diagnostic: the static law
(P0,P1,P2)=(1/4,1/2,1/4) gives Z(s)=cosh(s/2)^2. It can be specified with a purely
diagonal static representation. No propagation operator or nondiagonalizable
state evolution is required.

## 3. Extensive topological bias is a different ensemble

Suppose for one fixed width and fixed p in (0,1) the established sector laws give

    P0(m)=a0 lambda0^m(1+o(1)),
    P2(m)=a2 lambda2^m(1+o(1)),
    P1(m)->1,

where 0<lambda0,lambda2<1 and a0,a2>0. Write alpha_j=-log lambda_j.
Set s=m*sigma, not fixed s. Then the finite-sum Laplace principle gives

    lim (1/m)log Z(m*sigma)
       = max{ -alpha0-sigma, 0, -alpha2+sigma }.             (1)

Three regions follow:

    sigma < -alpha0: tilted rank -> 0;
    -alpha0 < sigma < alpha2: tilted rank -> 1;
    sigma > alpha2: tilted rank -> 2.

The rank-1 interval is nonempty. There is no direct rank-0/rank-2 transition at
sigma=0, including at the cylinder crossing lambda0=lambda2. The slopes of (1)
are -1,0,+1, exactly the corresponding X values. At each boundary an O(1/m)
window has an elementary two-weight logistic crossover (with the actual a_j).

The two source zero families satisfy

    Re(s_-)/m -> -alpha0,
    Re(s_+)/m -> +alpha2,
    Im(s_+/-)/m = pi/m  modulo 2pi/m.

They can pinch the real SIGMA axis, because sigma=s/m is a changing coordinate.
This does not contradict the fixed-s zero-free strip. It is an extensive bias
selecting exponentially rare topological sectors, not evidence of a new bulk
percolation transition at fixed topological chemical potential.

## 4. Executed width-four controls

The PR710 all-p trace factors were reused unmodified, with SHA-256
5edc624377399ad0878c7a608e722ebb65a7be454a6f022a1b6382a0661e91b0.
The new script reconstructs normalized factor traces by Newton identities.
Probabilities and zeros are evaluated at 100 digits, with relative polynomial
residuals checked. Selected results are independently recomputed at 140 digits.

At q4=0.59141717085313848179883410173592317796...,

    alpha0=alpha2=0.168595883942974875337464554747825...

The source zeros change from a complex pair at m=4, approximately

    s=-0.0057482260 +/- 2.1520920624 i,

to negative-u-root families. At m=64 their real parts are approximately
+/-10.7900954430, and at m=256 they are +/-43.1605462894; both have Im s=pi
on the chosen principal-log branch. The source drift per row tends to +/-alpha.
These numbers come from tiny factor recurrences, not random samples.

The general formulas were also checked at p=1/2 and p=7/10. Corresponding negative
and positive extensive-source boundaries are respectively

    p=1/2: -0.0828100231812642, +0.308695661403879;
    p=7/10: -0.350400958684780, +0.0646998760378965.

The complete result records all zeros and tilted probabilities for m=4,8,16,64,256.
Their numerical precision is not represented as rational interval certification.
No Kac label, universal exponent, new p_c value, or original-U source identification
is inferred. The zero-free theorem itself is exact and independent of PR710.
