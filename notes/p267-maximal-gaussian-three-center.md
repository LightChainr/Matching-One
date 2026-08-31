# A maximal Gaussian core with three positive moment centers

This is a post-reveal exploratory explanation of the existing `N100`/`N400`
`D_A` rank-step profiles. It does not change the independent `N900` width target.
The fixed construction and unused readouts are in
`experiments/p267_max_gaussian_three_center_20260831.json`.

## Definition before numerical evaluation

Normalize the signed rank-step profile by its area and write its centered,
unit-variance coordinate as `u`. Let `m_r` denote its moments, so that
`m_0=1,m_1=0,m_2=1`. This normalization does not assume that every empirical
rank bin is positive. A *positive moment realization* is an additional claim
to be checked, not an assumption about the source.

For a candidate common Gaussian variance `t`, define

\[
q_r(t)=\sum_{j=0}^{\lfloor r/2\rfloor}
\frac{(-t)^j r!}{2^j j!(r-2j)!}m_{r-2j}.
\]

These are the formal moments after inverse Gaussian convolution. Form
`H_3(t)=[q_{i+j}(t)]_{i,j=0}^3`. Seek the first boundary `t_*` of its PSD
region starting at zero. Feasibility is downward closed: if the functional
at `t` is positive on squares of cubics, then the functional at `s<t` is its
Gaussian convolution with variance `t-s`; averaging translated polynomial
squares preserves positivity. Thus this is one interval endpoint, not a scan
over disconnected mixtures. Its variance condition gives `t<=1`.

If `H_2(t_*)` is positive definite and `H_3(t_*)` is flat of rank three, the
truncated moments have a unique positive three-atom realization. Its monic
orthogonal cubic is

\[
P_3(x)=x^3+c_2x^2+c_1x+c_0,\qquad
H_2(c_0,c_1,c_2)^T=-(q_3,q_4,q_5)^T.
\]

The ordered real roots are the centers, and positive weights reproduce
`q_0,q_1,q_2`. Flatness ensures reconstruction through `q_6`. Convolving these
centers with `Normal(0,t_*)` gives the declared maximal-common-variance
three-center moment model through `m_6`. A failure of positivity, flatness,
or moment reproduction stops this construction; it does not trigger another
component count or family.

The moments `m_7,m_8` are then new algebraic predictions of this construction,
not inputs. They are read directly from the same original histograms and
compared using the full common-batch leave-one-out covariance, with the
entire construction rerun inside each replicate. Their status is **unused
moments in a post-reveal analysis**, not independent held-out data.

Three atoms are not three fields. Even exact agreement of finitely many moments
would not identify physical Gaussian components, and a finite-support
rank-step function cannot literally equal an unbounded Gaussian mixture with
positive variance. The purpose is a minimal positive explanation of the
observed changing shape, with a definite next failure direction.

## Numerical result

Pending the frozen calculation.
