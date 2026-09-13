# Full-profile common-coordinate collapse has a noncrossing transport test

2026-09-12. Mathematical extension of #706; the real-data result must be read
from the executed job before a verdict is claimed. Fixed points remain
p=Q290(.3,.5,.7), the same central set used in the reflection controls.

Suppose ONE increasing psi and ONE profile z satisfy psi(Q_i(u))=a_i+b_i*z(u),
b_i>0, for i=145,290,725 on a common domain. Define A=Q145 o F290 and
B=Q725 o F290, where F290 is inverse of the already combined spin0 Q290.
Then psi A psi^-1 and psi B psi^-1 are positive-slope affine maps. For any
affine maps a(x)=s*x+t and b(x)=r*x+v,

    a(b(x))-b(a(x))=(s-1)*v-(r-1)*t,

which is independent of x. Thus A(B(p))-B(A(p)) is either identically zero
or has one strict sign throughout the common domain. Increasing conjugacy
preserves sign, not magnitude. Opposite signs rule out common affine
conjugacy, and hence the specified exact common-coordinate location-scale
collapse. This is stronger than mere simultaneous reflection symmetry,
which need not pin the symmetric part of the profile.

IMPORTANT: nonzero alone does NOT reject. Aff(1) is nonabelian; two affine
maps with distinct fixed points need not commute. This would be the same
mistake as using symplectic orthogonality as Euclidean orthogonality: the
correct group identity, not a familiar label, carries the inference.

The script uses separate-size nonlinear deletion, including randomness in
the Q290 reference points, full stored covariance, and 55-dps pooled
Bernstein checks against a fast beta-mixture path. A fixed Bonferroni set
of three two-sided nominal intervals (total alpha=.0027) tests whether
opposite signs are resolved. It uses no inverse covariance, no fitted
chart, no covariance-selected contrast, and no search over points. The
entire composition and every deletion must remain in central quantiles
[.05,.95]. A failed domain or precision condition is a limitation, not
permission to choose nicer points. The result is C2 on existing data.

No resolved crossing does not prove a common coordinate. A future power
claim needs a specific alternative, not repeated null-compatible tests.
Three mathematical controls cover affine noncommutation, nonlinear
conjugacy, and an actual crossing of two increasing non-affine maps.
