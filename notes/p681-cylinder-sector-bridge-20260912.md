# #681 reviewed: an exact width-two torus/cylinder bridge

2026-09-12. Supersedes the first #705 retrieval note's proposed O1. Its useful
primary-source reading is retained below; the first version remains in Git
history. The missing small calculation is now completed, not a new task.

## 1. Exact local object, not a tautological configuration-diagonal operator

Take the axis square-site torus with periods (2,0),(0,m), m>=2. It has 2m sites;
every square has four distinct corners. Retain the two periodic edges between
the same two row vertices as different lifted edges, not a collapsed simple
edge. The occupied graph is NN, the complement matching graph NN+NNN. Use the
existing digital-Alexander observable M=P2-P0=E[r_black-1].

Put x=p(1-p), y=p^2. A row is empty, left-only, right-only or both, with weights
(1-p)^2,x,x,y. On the three NONEMPTY states (left,right,both), define

    K = [[1,0,1], [0,1,1], [1,1,1]],
    T = K diag(x,x,y) = [[x,0,y], [0,x,y], [x,x,y]].

This is a local row transfer. tr(T^m) is the total Bernoulli weight of cyclic
nonempty row words with overlap between every pair of consecutive rows.

**Lemma (event dictionary).**

    P2 = tr(T^m) - 2*x^m,
    P0 = (1-p^2)^m - 2*x^m,
    M_{2,m}(p) = tr(T^m) - (1-p^2)^m.                 (1)

Proof. Without a fully occupied row there are no horizontal occupied edges.
A nonzero longitudinal cycle then exists only if every row is the same
singleton, left or right. Consequently no horizontal cycle and no longitudinal
cycle has weight [(1-p)^2+2x]^m-2x^m=(1-p^2)^m-2x^m.
If there is a both-occupied row, a transverse cycle is present. A longitudinal
cycle exists exactly when every interface has occupied overlap: an empty row
or adjacent opposite singletons is a cut; conversely overlapping interfaces can
be joined within each full row into a closed longitudinal walk. With no full
row, the two all-singleton words are rank one, so subtract 2x^m from tr(T^m).
The cycles coexist in a component (or directly use intersection), giving rank
2. This proves the configuration classification and hence (1), at every m>=2.
No census is used in the proof.

For a fugacity convention v=p/(1-p), (1+v)^(2m) M is the corresponding signed
site-count polynomial. Do not confuse this site fugacity with an independent
FK bond variable or change the number of local cells without a geometry map.
The presentation uses normalized Bernoulli weights throughout.

## 2. The complete spectrum settles the finite-length question

The antisymmetric row vector has eigenvalue x. On the left/right-symmetric
subspace the transfer is [[x,y],[2x,y]], with characteristic equation

    lambda^2 - p*lambda - p^3(1-p) = 0.

Thus

    lambda_± = p/2 [1 ± sqrt(1+4p-4p^2)],
    M_{2,m} = lambda_+^m + lambda_-^m + x^m - lambda_c^m,
    lambda_c = 1-p^2.                                (2)

For 0<p<1, T is nonnegative, irreducible and aperiodic. Its Perron eigenvalue is
lambda_+. The two LEADING coefficients in (2) are exactly one. Yet the finite
remainder is x^m+lambda_-^m. Equal leading coefficients therefore do NOT imply
M=lambda_+^m-lambda_c^m. This refutes the first note's O1 equivalence inside the
actual square-site model, not just with an abstract matrix counterexample.

## 3. Exact cylinder crossing and root displacement

Setting lambda_+=lambda_c gives

    (p-1)(2p^3+2p^2-1)=0.

The unique interior solution q is the root of

    2q^3+2q^2-1=0,
    q = 0.56519771738363939643752801324703081609848397675955...

It agrees with Jacobsen (2015), Table 2, n=2 to the printed precision. This is
an explicit finite-width probability/spectrum bridge. We have not constructed
an all-width intertwiner with Jacobsen's augmented pTL representation, and do
not promote the numerical match alone to such an intertwiner or to novelty.

At q, lambda_-=q-lambda_c<0 and

    |lambda_-|/x = q^2/(1-q^2) < 1.

Therefore M_{2,m}(q)=x^m+lambda_-^m>0 for EVERY m>=2. Since M is strictly
increasing (the existing monotone-rank argument), its unique interior root
p_{2,m} is strictly below q for every finite m. This gives a particularly clear
failure of finite-root equality despite equal leading coefficients.

Let h(p)=log(lambda_+/lambda_c), r=x(q)/lambda_c(q)=q/(1+q). Then h'(q)>0 and

    p_{2,m}-q = - r^m/[m h'(q)]
                 * [1+(-q^2/(1-q^2))^m+O(r^m)].     (3)

Here r=0.361103080528647... and h'(q)=3.353388815848793.... Divide (2) by
lambda_c^m, expand e^(m h) around q, and use the two uniform subleading ratios
strictly below one. The root displacement is O(r^m/m); differentiating the
subleading terms and the leading exponential changes the relative remainder
by O(r^m). This also proves convergence at fixed width without exchanging limits.

Selected exact-root diagnostics (rational isolation intervals are in JSON):

| m | p_{2,m} | (p_{2,m}-q) / [-r^m/(m h'(q))] |
|---|---|---|
| 2 | 0.5411961001461969844 | 1.23450 |
| 4 | 0.5638649868188458323 | 1.05138 |
| 8 | 0.5651869150079729083 | 1.00241 |
| 12 | 0.5651975952151489695 | 1.000115 |
| 20 | 0.5651977173624504152 | 1.00000027 |

These are finite-width diagnostics, NOT estimates or bounds for the infinite
square-site p_c. Fixed width violates the expanding-geometry hypothesis of
#613; taking m->infinity here must not be confused with an all-directions limit.

## 4. General conditional lemma: what coefficients actually do

Suppose, at a fixed width near an isolated crossing p0,

    Z_o=c_o(p) lambda_o(p)^m [1+epsilon_o,m(p)],
    Z_c=c_c(p) lambda_c(p)^m [1+epsilon_c,m(p)],

with c_o,c_c positive and C2, simple positive leading eigenvalues, a uniform
spectral gap giving epsilon and its needed derivatives exponentially small,
and h=log(lambda_o/lambda_c), h(p0)=0, h'(p0)!=0. The nearby balance root obeys

    p_m-p0 = log[c_c(p0)/c_o(p0)]/[m h'(p0)]
               + O(m^-2 + rho^m/m).

Nonzero unequal coefficients can create a 1/m displacement without changing
the limiting crossing. Equal coefficients AT the crossing remove that 1/m
term; they do not remove subleading spectra. Width two is the explicit case
where those remaining terms and their coefficient are now known. An all-m
two-mode identity would require cancellation of EVERY other observable spectral
mode (including Jordan-polynomial terms), not merely equality of two prefactors.

No width-uniform estimate is established here. A general local bridge needs a
specified closure and its full observable spectrum; an absent literature
formula is not a no-go theorem. The fixed-width correction above does not
supply any n^-4 or L^-4 outer-limit theorem.

## 5. Primary reading retained, with corrections

Jacobsen, J. Phys. A 48 (2015) 454003, arXiv:1507.03027v1:
https://arxiv.org/html/1507.03027v1 — PRIMARY_TEXT_READ, §§2–4,6.1,8–9.
Eqs. (4),(9),(11)–(13) describe the signed graph polynomial and the fixed-width
cylinder eigenvalue method. Eq. (32) gives the square-site local loop operator.
Table 2 supplies the n=2 comparison above. Eq. (50) is an OBSERVED convergence
law: the text after (50) explicitly says more ingredients are needed to deduce
it from (49). It is not a theorem supplied by that CFT argument.

The spin-twist Eq. (55) is not a literal q=1 bridge: q=1 has only twist zero,
and its factor (1-1/q) vanishes, leaving 0=0. Work in the FK/loop construction
at q=1 unless an actual continuation is provided. Eq. (24)'s width-one Potts
example is not by itself a theorem for all square-site local operators.

Mertens–Ziff, PRE 94 (2016) 062152, arXiv:1603.07289v2:
https://arxiv.org/html/1603.07289v2 — PRIMARY_TEXT_READ, §II, Eqs. (20)–(21).
These give the finite matching/event identity and its all-equals-none relation;
the quoted root exponent is empirical. Jacobsen 2024 Reply remains
ABSTRACT_ONLY in the original retrieval; no claim here needs its unavailable body.

## 6. Executed checks and decision

`scripts/width2_cylinder_exact.py` derives integer coefficients by the trace
recurrence. An independent lifted-homology traversal enumerated all 5456
configurations across 2x2,...,2x6 and reproduced EVERY coefficient. It preserves
parallel periodic edges; no row-compatibility code is used by that verifier.
The 2x2 control is -1+4p^2-2p^4. A second integer matrix-power check at p=1/2
agrees exactly. Three local tests passed. Finite roots carry 100-bisection exact
rational brackets; spectral decimals are explicitly diagnostic.

Result: `results/research-control-20260912/width2-cylinder-exact.json`.
The suggested smallest #681 calculation is done and its premise corrected.
Do not commission it again. This is a useful local theorem, not a reason to
build a large transfer engine before defining an all-width scientific target.
