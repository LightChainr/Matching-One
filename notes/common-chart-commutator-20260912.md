# A parameter-free necessary condition for one common symmetrizing coordinate

2026-09-12. Mathematical follow-up to the completed #706 quadratic-chart
control. No real-data commutator has been computed here. This supplies a
specific falsifier rather than another fitted polynomial family. No novelty
claim is made.

## Object and lemma

Let Q_i be strictly increasing quantile functions with inverse F_i, all on the
same probability interval. Define their decreasing reflection involutions

    R_i(p)=Q_i(1-F_i(p)),        R_i o R_i = id.

For the current spin0 object, Q_i means the ALREADY COMBINED quantile function.
F_i is its inverse, not the orientation-weighted mixture of physical CDFs.
The test is therefore about the declared effective quantile laws. Nonlinear
coordinate changes need not commute with orientation weighting, so failure
cannot be silently transported to a different uncombined physical observable.

Suppose one increasing coordinate psi simultaneously symmetrizes all three
laws, possibly about different centers c_i. Precisely, on the common domain,

    psi(R_i(p)) = 2*c_i - psi(p).

Then T_ij=R_i o R_j is conjugate under psi to translation by 2(c_i-c_j).
Translations commute, hence the following NECESSARY identity:

    (R1 R2)(R2 R3) = (R2 R3)(R1 R2),
    R1 R3 = R2 R3 R1 R2.                              (1)

Products denote composition, rightmost first. This is the proof; neither a
polynomial truncation nor a fitted scaling exponent is involved. The residual

    C(p)=R1(R3(p))-R2(R3(R1(R2(p))))                   (2)

must vanish wherever both compositions are defined. Under a common increasing
coordinate change all maps are conjugated, so equality (1), or its failure,
is coordinate-independent. The NUMERICAL MAGNITUDE of (2) is not invariant;
keep it on the declared raw p axis and propagate its uncertainty there.

This is necessary, not sufficient. A zero at three points does not prove a
functional identity. Even identity can be non-identifying in degenerate cases:
R1=R2 makes (1) true for ANY R3. Nearly coinciding maps/centers can make the
commutator too small to measure. An unresolved commutator is not evidence for
one universal chart.

## Exact control and a genuine counterexample

With the shared logit coordinate, R_k(p)=k(1-p)/(p+k(1-p)), k>0. Using k=1,2,3,
(1) holds exactly at rational p=.3,.5,.7; tested with Fraction arithmetic.
Replace the third map by (1-p^3)^(1/3). It is still a decreasing involution,
but C(.5)=0.01884049384875336739614... is nonzero. Thus individual ability to
symmetrize each distribution does NOT imply a common symmetrizing chart.
The comparison kernel and three local tests are in `scripts/shape_common_chart.py`
and `tests/test_shape_common_chart.py`.

## A bounded real-data test, not yet executed

Use R1,R2,R3 for the existing N145/290/725 spin0 quantiles. Fix three reference
points p_j=Q290(u_j), u_j=.3,.5,.7. Retain the randomness of these reference
points when deleting the N290 block. Evaluate both sides of (1), each R_i's
involution control, the transformed intermediate arguments and the full joint
covariance under separate deletion of each independent size. Never pair batch
labels across different seeds. The existing source seeds/observable contract
remain unchanged.

A nested inversion is allowed; an uncontrolled interpolation is not. Compare
the existing float path with higher-precision/bracketed inversion on the pooled
points first and bound the composed numerical error. Profile one deletion before
committing CPU time. If intermediate quantiles leave the declared central range
or numerical error is comparable to statistical error, return that limitation
rather than selecting different favorable points. Store the joint result and
label every covariance reference nominal. No new Monte Carlo is necessary.

One resolved failure rules out a single common symmetrizing coordinate for
THESE effective finite-size laws on the tested domain. It does not identify a
CFT field, rule out asymptotic symmetry with finite corrections, or establish an
exponent. Zero/unresolved gives no positive identification. This distinction is
more informative than allowing higher polynomial degree after each failed fit.
