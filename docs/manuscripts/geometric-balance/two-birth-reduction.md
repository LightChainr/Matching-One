# Two sharp births, not an arbitrary broad threshold law

2026-09-13. A continuation of the **same** geometric-balance manuscript (#739).
The new organizing statement here is a direct application of the published
Friedgut--Kalai sharp-threshold theorem, followed by finite identities and
consequences. It is not presented as a new sharp-threshold theorem or as a
novelty certificate. No new continuum model, source dictionary, or width scan
is introduced.

## 1. The next question, and what is settled here

The main manuscript distinguishes the balance root from concentration of the
whole birth-time mixture. It gives the full-law geometric criterion
`log N / ell -> 0`, and endpoint splitting for `log N / ell -> infinity`.
The intermediate regime is not a reason to fit an unrestricted profile.

For an honest integer-period square-site torus, write

\[
 f_1(p)=\Pr_p(r\ge1)=1-P_0(p),\qquad f_2(p)=\Pr_p(r=2)=P_2(p).
\]

Both are increasing nonconstant Boolean-event probabilities on **N independent
site variables**. The translation group `Z^2/Lambda` acts transitively on
those variables and preserves either event. Translation acts trivially on
ambient homology, but preserves its rank; no rotation of the period lattice
or irreducibility assumption on its homology representation is needed.

Let `T1 <= T2` be the two birth times in the usual uniform-label coupling.
Define their unique medians by

\[
 f_1(a_N)=\tfrac12,\qquad f_2(b_N)=\tfrac12.
\]

Finite monotonicity gives `a_N <= b_N`. The mixture CDF and its median are

\[
 F_N=\tfrac12(f_1+f_2),\qquad q_N=F_N^{-1}(1/2).
\]

The main conclusion below is uniform over the period shape:

\[
 W_1\!\left(\mathcal L(T),\tfrac12\delta_{a_N}+\tfrac12\delta_{b_N}\right)
 \le \frac{C_{\rm FK}}{\log N}.                                      \tag{1}
\]

Here `T` is a fair independent selection of T1 or T2. The constant is the
universal constant in the imported theorem; it is NOT estimated by our tiny
controls. At finite precision and small N the bound may be uninformative.
Equation (1) is in the **unscaled probability coordinate p**. It does not imply
an affine closure of standardized quantile shapes or of finite-size response
jets. In particular it does not undo #706's finite nonaffine shape result.

## 2. Precisely the imported theorem

Friedgut and Kalai (1996), in the theorem stated on the first page and in the
publisher abstract [FK], prove that for a monotone event A invariant under a
transitive group on N Boolean coordinates, there is an absolute C_FK such that

\[
 \mu_p(A)>\epsilon,\qquad
 q\ge p+C_{\rm FK}\frac{\log(1/(2\epsilon))}{\log N}
 \quad\Longrightarrow\quad \mu_q(A)>1-\epsilon,                       \tag{2}
\]

when the displayed parameters lie in [0,1]. This exact formulation is also
stated as Theorem 6 of Duncan--Kahle--Schweinhart [DKS]. We read the publisher's
statement and DKS section 1.3; the original AMS PDF was not fetched successfully.
The general theorem, not DKS's model-specific conclusion about synchronized
homological thresholds, is used here.

Only product measure, event monotonicity and coordinate transitivity are
needed for this application. Site sharpness, RSW and matching duality are
**not** inputs to equations (1)--(8). They remain inputs to the separate
geometric conclusions in the main manuscript.

## 3. Concentration of each birth, with a rate

Put `L=log N`, let A be either rank event, and let theta be its median. For x>0
choose `epsilon = exp(-L*x/C_FK)/2` in (2). Starting from the median bounds
the upper tail. If `mu_(theta-x)(A)>epsilon`, (2) would give
`mu_theta(A)>1-epsilon>1/2`, a contradiction, bounding the lower tail. A tail
whose parameter is outside [0,1] is simply zero. Thus

\[
 \Pr(T_j\le\theta_j-x),\;\Pr(T_j>\theta_j+x)
 \le\tfrac12 e^{-Lx/C_{\rm FK}},\qquad
 (\theta_1,\theta_2)=(a_N,b_N).                                      \tag{3}
\]

Consequently, for every positive integer k,

\[
 \mathbb E|T_j-\theta_j|^k
 \le k!\left(\frac{C_{\rm FK}}L\right)^k.                            \tag{4}
\]

This follows by integrating the tail: `E X^k = integral k x^(k-1) Pr(X>x) dx`.
In particular the mean-median discrepancy is at most C_FK/L and the variance
is at most `2 C_FK^2/L^2`.

These estimates do **not** assert independence of T1 and T2. Their joint law
is within `2 C_FK/L` of the point `(a_N,b_N)` in W1 with the l1 metric,
by their actual common-label coupling. Correlations of rescaled fluctuations
can survive. The finite 2x3 control has covariance `1123/58800`, not zero.

Couple T to the corresponding median using the same fair selector. Equation
(4) proves (1), and similarly W_k is at most `(k!)^(1/k) C_FK/L`.
Every convergent subsequence `(a_N,b_N)->(a,b)` therefore satisfies

\[
 \mathcal L(T_1,T_2)\Longrightarrow\delta_{(a,b)},\qquad
 \mathcal L(T)\Longrightarrow\tfrac12\delta_a+\tfrac12\delta_b.       \tag{5}
\]

Conversely any subsequential weak limit of the mixture has this form, by
compactness of the two centers. Centers need not converge on an arbitrary
oscillating geometric sequence. The two atoms may coincide. This is not an
assertion about the number of physical fields or Markov states.

There is also a path statement. If
`r_det(p)=1_{p>=a_N}+1_{p>=b_N}`, then under the same uniform labels

\[
 \mathbb E\int_0^1 |r(p)-r_{\rm det}(p)|\,dp
 \le\frac{2C_{\rm FK}}L.                                             \tag{6}
\]

This is an integrated-p bound, not a uniform-in-p approximation at jumps.
For `a_N+x <= p <= b_N-x`, (3) gives `P1(p)>=1-exp(-L*x/C_FK)`.
The rank-one plateau is caused by separation of two individually sharp
transitions, not by a broad transition of either Boolean event.

## 4. A finite data reduction already available in rank-birth archives

Set `Delta_N=b_N-a_N`, `G_N=integral_0^1 P1(p) dp`, and let K1,K2 be the
first and second ambient-rank birth **occupation indices** in a uniform
random permutation of the N sites. Then

\[
 G_N=\mathbb E(T_2-T_1)
    =\frac{\mathbb E(K_2-K_1)}{N+1},\qquad
 |G_N-\Delta_N|\le\frac{2C_{\rm FK}}L.                              \tag{7}
\]

The first equality is the survival-function identity. Conditional on the
permutation, the kth uniform order statistic has mean k/(N+1), giving the
second. It handles a simultaneous rank jump. K means occupation index here,
not elapsed physical time, and neither birth is replaced by a directional
wrapping proxy.

No new source or simulation is needed to read G_N from an archive that
already stores these correctly typed birth indices. For occupation-count
censuses `c_(1,k)`, the same exact identity is

`G_N = (1/(N+1)) sum_(k=0)^N c_(1,k)/binom(N,k)`.

The mixture quartiles give another version of the same asymptotic gap:

\[
 a_N-\frac{C_{\rm FK}\log2}{L}\le Q_N(1/4)\le a_N,
 \qquad b_N\le Q_N(3/4)\le b_N+\frac{C_{\rm FK}\log2}{L},
\]
\[
 \Delta_N\le\operatorname{IQR}(F_N)
 \le\Delta_N+\frac{2C_{\rm FK}\log2}{L}.                              \tag{8}
\]

For example, `f1/2 <= F <= f1` brackets Q(1/4) between the first-birth
quarter-quantile and median. The symmetric argument brackets Q(3/4).
The factor log2 then comes directly from (2).

These are asymptotic absolute-error relations, not a license to equate finite
thermal jets, IQR and the rank gap without their finite errors. C_FK is not
a calibrated finite-sample error bar.

For completeness the mixture variance separates exactly:

\[
 \operatorname{Var}(T)=\tfrac12\operatorname{Var}(T_1)
 +\tfrac12\operatorname{Var}(T_2)
 +\tfrac14[\mathbb ET_2-\mathbb ET_1]^2.
\]

Equations (3)--(4) imply an error of at most
`C_FK*Delta_N/L + 3 C_FK^2/L^2` from `Delta_N^2/4`.
A macroscopically broad mixture need not mean either birth is noisy.

## 5. What this adds to the geometric main theorem

For every finite torus,

\[
 a_N\le q_N\le b_N.                                                   \tag{9}
\]

At a_N, `F<=1/2`; at b_N, `F>=1/2`. Thus when ell_N->infinity, Theorem A
of the manuscript puts p_c between every limiting pair of birth centers.
It does **not** make q_N their midpoint. A plateau discards the exponentially
small tail odds that can determine the finite median.

Combining the reduction with the geometric theorem in #739 gives, for
N_n->infinity,

\[
 \frac{\log N_n}{\ell_n}\to0
 \quad\Longleftrightarrow\quad \Delta_n\to0
 \quad\Longleftrightarrow\quad G_n\to0
 \quad\Longleftrightarrow\quad \operatorname{IQR}(F_n)\to0.           \tag{10}
\]

Here is the potentially missing step in the reverse direction: if ell has a
bounded subsequence, #739's forced-path packing gives T1->0 and T2->1 on that
subsequence, so Delta->1. Hence Delta->0 forces ell->infinity. Equation (9)
and root consistency then force both centers to p_c, and Theorem B yields
the geometric criterion. These geometric inputs are author-supplied results
of #739; equation (10) is not an independent validation of its corridor proof.

At extreme elongation the pair is (0,1). At finite positive log N/ell we now
know the possible shape of every subsequential unscaled law: two equal atoms,
with an undetermined pair of locations. The main remaining problem is the
location of those centers, not an unrestricted limiting profile.

## 6. Why transitivity does not synchronize the births

Translation is transitive on sites but acts trivially on homology. It sharpens
each rank event separately. It does not require its two threshold locations
to coincide. DKS explicitly separates transitivity from the extra point-group
symmetry of homology in its introduction and surjectivity argument [DKS].

A short finite comparison illustrates the extra symmetry. Suppose the actual
period lattice is invariant under a quarter-turn J. Split rational homology
lines into pairs `{l,Jl}`, choosing one line from each pair for a class A,
and its rotate for class B. Define E_A/E_B to mean the ambient image contains
a line in the respective class. These events are increasing. Their union is
`r>=1`, their intersection is `r=2`, and quarter-turn symmetry gives equal
probabilities. Harris association therefore gives, at every p,

\[
 P_2(p)\ge\Pr(E_A)\Pr(E_B)\ge\tfrac14[1-P_0(p)]^2.                  \tag{11}
\]

This is the two-dimensional elementary instance of DKS's symmetry strategy,
not a new general surjectivity theorem. At a_N it implies `P2>=1/16`.
Applying (2) with epsilon=1/32 shows

`b_N-a_N <= C_FK log(16)/log N`.

Such quarter-turn symmetry is real on Gaussian ideal quotients. It is absent
on a generic elongated or tilted period lattice. Changing a period basis does
not create it. This explains why importing a synchronized-threshold theorem
from symmetric tori would answer a different question.

## 7. A quantitative first step in the unsolved intermediate regime

For axial periods (w,0),(0,m), let w->infinity and `log m / w -> d` with
`0<d<infinity`. We can bracket every limiting center without assuming a
correlation-length exponent or fitting a new root law.

A nonzero cycle contains a simple graph cycle with at least w occupied vertices.
For a graph of degree at most D, the nonbacktracking walk count gives, when
`(D-1)p<1`,

\[
 \Pr_p(r>0)\le
 N\frac{D}{D-1}\frac{[(D-1)p]^w}{1-(D-1)p}.                           \tag{12}
\]

It is an overcount, which is appropriate for this upper bound. It also works
for the matching graph: each diagonal changes either coordinate by at most
one, so a nonzero axial period still requires at least w steps. D is 4 on NN
and 8 on matching. Fully occupied horizontal rows, on disjoint site sets,
give the other bound `P0 <= (1-p^w)^m`.

Using the exact matching-complement relation between birth medians, every
subsequential pair (a,b) obeys

\[
 \frac{e^{-d}}3\le a\le\min\{p_c,e^{-d}\},\qquad
 \max\{p_c,1-e^{-d}\}\le b\le1-\frac{e^{-d}}7.                       \tag{13}
\]

For example d=log4 gives `1/12 <= a <= 1/4` and
`3/4 <= b <= 27/28`. These are proven broad brackets, not predictions of
exact center values. The constants 3 and 7 are walk-count bounds, not measured
surface tensions. The relations do NOT establish convergence of a_n,b_n for
all sequences with the same d, and do not identify the equality case in a
putative cost-versus-volume transition.

## 8. Finite controls, not additional asymptotic evidence

The existing width-two local matrix gives

`x=p(1-p), lambda_±=p(1±sqrt(1+4p(1-p)))/2`,

`P0=(1-p^2)^m-2*x^m`,

`P2=lambda_+^m+lambda_-^m-x^m`.

The companion script uses those already established #705 formulas at
m=2,4,8,16,32,128. This is not a new state engine or a new circumference.
Selected values (full values are in the result JSON):

| m | a_m | q_m | b_m | integral P1 |
|---:|---:|---:|---:|---:|
| 4 | 0.39543672 | 0.56386499 | 0.71750418 | 0.29523810 |
| 16 | 0.20590464 | 0.56519772 | 0.86713870 | 0.63896753 |
| 128 | 0.07348861 | 0.56519772 | 0.95578173 | 0.87423968 |

Width 2 is fixed, so q_m is not approaching the infinite square-site p_c.
The example illustrates how separated sharp births and an internal balance
root coexist; it is not an extrapolation test of (10).

Independent graph-potential censuses of 2x2,2x3,2x4 (336 configurations total)
check the formulas, exact beta integrals and normalization. All 720 orders on
2x3 independently verify `E(K2-K1)/(N+1)=3/14` and retain nonzero birth-time
covariance. Seven targeted tests pass. Numerical roots/quadratures are computed
at 80 digits and selected values are recomputed at 110; they are not interval
certificates. Full repository CI has not been run for this addition.

## 9. One next research target

Keep the work in #739. The next target is whether the two centers have limits
on exponentially elongated sequences and, when they do, which **microscopic
winding cost** determines them. Start with the axial finite-d family and only
then ask about orientation dependence. Do not infer these costs from a
fixed-width eigenvalue or from a continuum cusp formula with another order of
limits. Equations (12)--(13) are the completed elementary start, not a solution
to the center-location problem.

No new production is requested here. Existing ambient-rank birth archives can
already supply a_N,b_N and G_N where both births are recorded; do not substitute
directional wrapping times. The balance-root theorem and its independent
proof review remain intact. Higher source orders, more Jordan examples, and
new generic control machinery are not needed to answer this next question.

## References and claim scope

[FK] E. Friedgut and G. Kalai, *Every monotone graph property has a sharp
threshold*, Proc. AMS **124** (1996), 2993--3002,
https://doi.org/10.1090/S0002-9939-96-03732-X . The exact transitive-event theorem
is printed in the publisher abstract:
https://www.ams.org/journals/proc/1996-124-10/S0002-9939-96-03732-X/ .
Publisher theorem statement read; full original PDF not retrieved this round.

[DKS] P. Duncan, M. Kahle and B. Schweinhart, *Homological percolation on a torus:
plaquettes and permutohedra*, arXiv:2011.11903v4,
https://arxiv.org/html/2011.11903v4 . Theorem 6 is precisely (2); section 1.1
and the section-3 surjectivity discussion distinguish coordinate transitivity
from the homology point-group argument. Those sections were read; their
model-specific simultaneous-threshold theorem is not imported to arbitrary
integer-period square-site tori.

The generic sharpness input, phase separation through multiple Boolean
events, order-statistic integral and symmetry mechanism are prior tools.
This note supplies their explicit consequences for the manuscript's exact
observable and narrows the remaining research question. No priority claim,
new critical exponent, original-U identification, or complete finite-profile
closure is made.
