# Width four: exact all-p sector traces, visible cancellations, and a cylinder crossing

Date: 2026-09-12. Completed continuation of #636 / draft #708. No new engine,
width scan, Monte Carlo or GPU work. No novelty claim for the published crossing
value or for the general linear-algebra tools.

## 1. The missing all-parameter identity is now established

Use #708's exact rank-future automaton for the axis square-site torus
`(4,0),(0,m)`, physical length m>=2. Let `t=p/(1-p)` be SITE fugacity.
Each row mask s has weight `t^popcount(s)`. An entire configuration with m
rows is normalized by `(1+t)^(4m)`, not `(1+t)^(4(m-1))`.

From its 509-state quotient, three explicit nonnegative local submatrices
B5(t), B15(t), B16(t) of sizes 5, 15 and 16 give

    (1+t)^(4m) P0 = tr(B5^m)  - tr(B16^m) + 2 t^(2m),
    (1+t)^(4m) P2 = tr(B15^m) - tr(B16^m) + 2 t^(2m),

and consequently

    (1+t)^(4m) M = tr(B15^m) - tr(B5^m).                 (1)

These are exact polynomial identities for EVERY m>=2 and every p in [0,1]
(endpoints by polynomial continuity). A block trace alone is not a probability;
P0 and P2 require the common subtraction and restoration terms displayed above.
The signed trace is NOT just the difference of two leading eigenvalue powers.

The matrices are fully specified by sparse polynomial rows in
`results/research-control-20260912/width4-parametric-traces.json`. The selector
indices in the upstream 509-state table are:

* B5: `[0,16,17,18,19]`;
* B15: `[15,139,140,141,142,143,144,145,146,147,148,149,150,151,152]`;
* B16: `[1,21,22,23,24,25,153,154,155,156,157,158,160,161,163,481]`.

Precisely, `(Bd)ij` is the sum of `t^popcount(s)` over row masks whose
upstream successor of selector i is selector j; transitions leaving the set
are omitted. This is a concrete local operator, not a diagonal list of all
configurations. The entire construction is reproducible from the pinned source
blob `50b7297deefe7c50215aea2ed534ca5810461af3` at #708 head
`f782061c1a592ed2f9fd0e9dabaa45f0e54bc4e7`.

For example, putting a=1+4t+6t^2, b=1+4t+5t^2,

    B5 = [[a,t^3,t^3,t^3,t^3],
          [b,t^2+t^3,t^3,0,t^3],
          [b,t^3,t^2+t^3,t^3,0],
          [b,0,t^3,t^2+t^3,t^3],
          [b,t^3,0,t^3,t^2+t^3]].

This is not an all-width intertwiner with Jacobsen's pTL representation. It
settles the finite width-four probability/closure question in its own explicit
representation.

## 2. Why finitely many calculations prove (1) for all m and t

First, common coefficientwise strong lumping of the 509 states gives
`3 -> 35 -> 94 -> 94` blocks. Rates into each block agree separately for every
row occupation count k=0,...,4, on ALL 509 source states. Hence this quotient
preserves P0 and P2 for every fugacity, including a prescribed varying fugacity
from row to row. It is not claimed minimal among arbitrary positive realizations;
it differs from #708's deterministic 509-class rank language.

The 94-state scalar sector sequence has a linear realization of dimension 94
over Q(t). The sequence `tr(Bd^m)` has recurrence order at most d by
Cayley-Hamilton. The monomial `t^(2m)` has order one. Thus a difference between
the proposed P0 formula and the actual sequence has order at most
`94+5+16+1=116`; for P2 the bound is `94+15+16+1=126`.

The verifier proves the first 126 differences identically zero as INTEGER
POLYNOMIALS, not at a fitted grid of parameter values. It uses radix
`T=2^512`. At length m each polynomial has degree at most 4m; an l1 coefficient
bound on the P2 difference is `34*16^m` (24 for P0). For m<=126 this is below
T. If an integer polynomial with all coefficients strictly between -T and T
evaluates to zero at T, its constant coefficient is a multiple of T and must
be zero; divide by T and repeat. Evaluation is injective under this proved
bound. This is exact Kronecker arithmetic, not numerical interpolation.

The 126 polynomial zero identities then imply every subsequent identity by
the finite-dimensional recurrence over Q(t). Matrices and weights are polynomial,
so the identities extend to all parameter specializations, including ones where
a generic minimal realization degenerates.

Small-block characteristic factors are independently checked against raw matrix
traces through orders 5,15,16; Newton identities determine their characteristic
polynomials. The factor-coefficient bounds are also below the same radix, so
these checks are polynomial identities. No assumption about numerical rank or
nearby parameter samples enters the proof.

## 3. Complete generic visible spectrum

All monic factors, with coefficients in Z[t], are in
`width4-parametric-definition.json`. Denote the quadratic factor of B5 by C2,
the two cubics and quintic of B15 by O3a,O3b,O5, and the cubic/two quartics of
B16 by S3,S4a,S4b. Exact factorizations are

    char(B5)  = (x-t^2(1+t))^2 (x-t^2(1-t)) C2,
    char(B15) = (x+t^2) O3a^2 O3b O5,
    char(B16) = (x+t^2) S3 S4a S4b^2.

For a monic f let tau_f(m) be its root power sum, multiplicities included.
Then the all-p matching spectrum is

    (1+t)^(4m) M = (-t^2)^m - 2[t^2(1+t)]^m - [t^2(1-t)]^m
                   - tau_C2(m) + 2 tau_O3a(m) + tau_O3b(m) + tau_O5(m).   (2)

Every visible weight is explicit. Repetition in a characteristic polynomial is
not a license to add a Jordan-polynomial term to this trace formula.

The generic minimal scalar recurrence orders over Q(t) are

    P0: 17,    P2: 23,    M=P2-P0: 16.

Upper bounds follow from (2) and the analogous sector factor lists; nonzero
17x17, 23x23, 16x16 physical-tail Hankel determinants at t=3 prove matching
lower bounds. Their integers are saved. The generic M recurrence is squarefree,
certified by a squarefree specialization at t=3 of its monic polynomial.

In fact the ordinary-trace formula has a stronger fixed-parameter implication:
for EVERY fixed interior p, the undifferentiated scalar tail is a finite sum
of pure eigenvalue powers. Even a defective block contributes only its
algebraic multiplicity times lambda^m to an ordinary trace; its nilpotent
part has zero trace. Thus no m*lambda^m term is required by this scalar law
at any fixed p. This does NOT imply that the 509-state operator, or even
every matrix block, is diagonalizable. It is an observable statement.

Specializations matter. At t=1 (p=1/2), `t^2(1-t)=0`, whose contribution vanishes
on physical lengths, leaving order 15, exactly #708's result. At t=2 (p=2/3),
`t^2(1-t)=-t^2`; two opposite-weight terms cancel and the order is 14. Explicit
nonzero Hankel minors prove both specialized lower bounds. A single p value
therefore does not determine the generic visible order. These are changes of
observability/cancellation, not changes of microscopic state number.

## 4. Perron crossing and its defining degree-17 polynomial

All three nonnegative blocks are irreducible and aperiodic for t>0. Positive
equitable quotients of sizes 2,5,4 respectively retain their Perron roots.
Every coefficientwise quotient identity is verified. The relevant closed
quadratic and open quintic are thus C2 and O5, not arbitrarily selected factors.

Write their Perron eigenvalues alpha_c(t), alpha_o(t). At t=0,
alpha_c -> 1 and alpha_o -> 0. At large t, alpha_o >= t^4 from its full-row
self-loop, while alpha_c <= 1+4t+6t^2+4t^3, so they must cross.
An exact 7x7 Sylvester determinant yields

    Res_x(O5,C2) = t^10 (t+1)^2 R17(t),

where, descending,

    R17 = 9t^17+198t^16+1305t^15+5205t^14+12915t^13+16534t^12+1795t^11
          -28383t^10-51259t^9-52627t^8-39092t^7-23218t^6-11396t^5
          -4503t^4-1344t^3-278t^2-35t-2.

There is exactly one positive root by Descartes' rule and endpoint signs.
The existing Perron crossing must be that root, and is unique. In physical p,

    q4 = 0.591417170853138481798834101735923177964270443192880...

It matches the ALREADY PUBLISHED width-four entry of Jacobsen 2015 Table 2.
This is neither a new infinite-square pc nor a certified interval for pc.
A rational bracket for t (160 bisections) and its p image are in the companion
certificate. Polynomial irreducibility/minimality is not asserted.

Near q4 the shared block's Perron root is strictly below both alpha_o and
alpha_c. This is certified by a positive rational Collatz vector, not inferred
from printed eigenvalues. Thus the leading weights of P0 and P2 are both
exactly one near this crossing. Equal leading weights still leave all the
subleading terms in (2).

## 5. Cancellation makes the observable relaxation rate different

At q4 the shared B16 Perron root divided by the common leading root is
approximately 0.773935047152260. It appears in BOTH P0 and P2 with coefficient
-1 and cancels from M IDENTICALLY for all t, not only at the crossing.

The largest remaining subleading root in M is the positive O3a root mu,
with coefficient two. Its ratio at q4 is

    rho = 0.2517497549919253805536005833510364513603... .

An estimate based on the slower common 0.774 mode would therefore give the
wrong finite-length rate for the signed observable. This is a concrete
square-site example of why the complete operator or either probability's
slowest correction is not automatically the relevant observable correction.

Let h(p)=log(alpha_o/alpha_c), accounting for t=p/(1-p). The common probability
normalizer cancels in this ratio. Rational intervals prove h'(q4)>0, with
numerical diagnostic `2.4643526094735922165...`. All other visible subleading
root moduli divided by mu are below theta=674/1000. Complex pairs are bounded
by rational Vieta product/discriminant intervals; they are not discarded.

After the leading cancellation at q4, the remaining trace sum is bounded below
by `mu^m [2-16 theta^m] > 0` for m>=6. Lengths 2,...,5 are checked by exact
polynomial interval signs. Hence

    M_(4,m)(q4)>0,  and p_(4,m)<q4, for EVERY m>=2.

The latter uses strict monotonicity of finite expected rank, not a fitted root
sequence. Taylor expansion of the locally separated eigenvalues gives

    p_(4,m)-q4 = -2 rho^m/[m h'(q4)] * [1+O(theta^m)+O(rho^m)].          (3)

This is a fixed-width result. It provides NO bound uniform in width and no
fixed-aspect L^-4 theorem.

## 6. A useful thermal-jet warning for #275

Differentiating the leading difference at the crossing gives

    partial_p M_(4,m)(q4)
      = m lambda_*^m h'(q4) + partial_p R_m(q4),

where lambda_*=alpha_o/(1+t)^4=alpha_c/(1+t)^4 is the normalized leading root.
The factor m comes from differentiating two crossing simple Perron branches.
The leading eigenspace of `diag(B15,B5)` is SEMISIMPLE: each block has a simple
Perron root and the blocks are uncoupled. Thus a polynomial-in-length factor
in a thermally differentiated trace can arise without a defective leading
operator.

A transparent control is `diag(lambda+e,lambda-e)` with weighted trace
`(lambda+e)^m-(lambda-e)^m`: its derivative at e=0 is
`2m lambda^(m-1)`, although the matrix is diagonal for every e, including zero.
This does not refute an LCFT mechanism. It says that a candidate-specific map
through normalization, differentiation and moving-root evaluation is needed
before such a polynomial factor can distinguish mechanisms. It strengthens
#275's existing map requirement; it does not reopen that issue with a third
post-hoc model or a new sampling order.

## Sources, execution, and boundaries

Repository source: #708 head f782061c, exact state-certificate blob stated
above. Earlier #705/#707 provide width-two/three controls; none was rerun as
a new campaign. The primary external comparison is Jacobsen,
*Critical points of Potts and O(N) models from eigenvalue identities in periodic
Temperley-Lieb algebras*, arXiv:1507.03027v1, sections 4 and 6.1, Table 2:
https://arxiv.org/html/1507.03027v1 . PRIMARY_TEXT_READ on 2026-09-12.
The published eigenvalue value and limit order are attributed to that source.

`width4_parametric_traces.py` verifies the full polynomial identities with only
the Python standard library. `width4_cylinder_certificate.py` uses mpmath to
propose root bands, then checks every load-bearing sign and inequality with
Fraction arithmetic. Actual execution times and local test results are in the
handoff validation file. No full-repository CI was run for these additions.
