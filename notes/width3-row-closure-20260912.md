# Width three is the last connected-row case: exact torus/cylinder closure and a width-four memory witness

Date: 2026-09-12. Direct continuation of the reviewed width-two result in PR #705.
This is a completed mathematical analysis with small independent checks, not a
request to build a large transfer matrix. No novelty claim is made. The printed
cylinder root was already published by Jacobsen; the deliverable is an explicit
finite event dictionary, complete observable spectrum, root-displacement formula,
and a precise reason that the same occupancy-only simplification stops at width four.

## 1. Scope, weights and the local transfer

Use the square-site NN graph on the honest axis torus with periods (3,0),(0,m),
m >= 2. Preserve parallel lifted edges when m=2. Let r be its ambient rational
homology rank. The existing digital-Alexander identity gives D=r-1 and M=E[D].
Write P_j=P(r=j). Occupation probability is p, not an FK bond fugacity.

There are seven nonempty occupied subsets S of a three-site row. Define

    u=p(1-p)^2,  v=p^2(1-p),  w=p^3,
    t(S)=p^|S| (1-p)^(3-|S|),
    T[S,S'] = 1{S intersects S'} t(S').

Thus tr(T^m) is the Bernoulli probability that every row is nonempty and every
successive pair, including the seam, has occupied overlap. T is nonnegative,
irreducible, and aperiodic for 0<p<1: every state has a loop and communicates
with the full row.

## 2. Configurationwise lemma: r=F+V

Let F indicate that at least one row is fully occupied. Let V indicate that every
periodic interface has occupied overlap. Then, at width three,

    r = 1_F + 1_V.                                      (1)

The proof must use the width restriction, not only a picture of overlapping rows.

If F occurs, a full row carries the transverse generator. If V also occurs,
choose one occupied edge across each interface. Within each nonempty row the
incoming and outgoing endpoints can be joined: every nonempty subset of a
three-cycle is connected. Joining around the torus gives a closed walk with
longitudinal winding one; its class is independent of the transverse generator.
Hence r=2. If V fails, an interface has no occupied vertical edge. Cutting there
puts the graph in an annulus, so longitudinal winding is impossible and r=1.

Now suppose F does not occur. Each nonempty occupied row is a point or a single
edge, hence a contractible tree. Contract these row trees. Two adjacent proper
rows have two common occupied columns only when they are the same two-site set.
The two vertical interface edges in that case differ, after contraction, by a
contractible plaquette boundary, not by an extra ambient cycle. After removing
these homologically redundant parallel edges, the row quotient is a subgraph
of the periodic row cycle. Distinct periodic interfaces are NOT identified with
each other when m=2. If V fails this quotient is a union of paths, giving r=0.
If V holds it has exactly one independent cycle, whose longitudinal winding is
one, giving r=1. Its transverse winding may be nonzero, but it remains rank one.
This proves (1). In particular, the argument does not confuse a spiral with a
rank-two cross.

Taking expectations immediately gives the all-length identity

    M_(3,m)(p) = tr(T(p)^m) - (1-p^3)^m.                 (2)

An equivalent sector-level description uses T_proper, the six-state restriction
with the full row omitted:

    P_2 = tr(T^m) - tr(T_proper^m),
    P_0 = (1-p^3)^m - tr(T_proper^m).

The proper-row contribution cancels exactly in P_2-P_0. This cancellation is not
permission to discard every subleading eigenvalue of T. Also
F_(3,m)(p)=E[r]/2=(1+M_(3,m)(p))/2 supplies an exact finite quantile laboratory,
not an additional large-N percolation experiment.

## 3. Complete spectrum, including the repeated sector

The row permutation group is S_3 (the automorphism group of the three-cycle).
On functions constant on subsets of each size, T reduces to

           [u   2v   w]
    B =    [2u  3v   w].
           [3u  3v   w]

The singleton and doubleton permutation representations each contain a
2-dimensional standard representation. Label a doubleton by its omitted vertex.
On their standard components the reduction is

    C = [u  -v; -u  0],

and each eigenvalue of C occurs twice in the full seven-state spectrum. Thus

    tr(T^m) = tr(B^m) + 2 tr(C^m),
    f_B(L)=L^3-(u+3v+w)L^2-(uv+2uw)L+uvw,
    mu_±=(u ± sqrt(u^2+4uv))/2.                          (3)

Let beta_1 be the Perron eigenvalue of B and beta_2,beta_3 its other roots. The
exact observable spectral formula is

    M_(3,m)=beta_1^m+beta_2^m+beta_3^m
            +2 mu_+^m+2 mu_-^m-(1-p^3)^m.              (4)

The leading open and closed coefficients are both one. The leading SUBLEADING
coefficient will be two, not one. This is ordinary semisimple multiplicity:
T=K diag(t) is similar to the real symmetric matrix sqrt(diag(t)) K sqrt(diag(t)).
There are no nontrivial Jordan blocks for 0<p<1. Neither a multiplicity of two nor
an alternating finite-length correction is evidence of a continuum Jordan pair.

## 4. Cylinder crossing, exact algebra and uniqueness

Set lambda_c=1-p^3. Substituting into f_B gives

    f_B(lambda_c)=(p-1)^2 P(p),
    P(p)=p^6-3p^5-5p^4-4p^3+p+1.                       (5)

There is exactly one root q in (0,1). One convenient proof uses s=p/(1-p):

    (1+s)^6 P(s/(1+s))
      =1+7s+20s^2+26s^3+8s^4-14s^5-9s^6.

In descending coefficient order this has exactly one sign change, so Descartes'
rule bounds the number of positive roots by one. Existence follows from P(0)>0
and P(1)<0. The rational certificate in the script isolates q and the three B
roots, showing that lambda_c(q) is the Perron root, not a subleading crossing.

    q=0.58888069991785299805144269575170493372217050345701...

Jacobsen (2015), section 6.1 Table 2, already prints
0.5888806999178529980514426957517049337221 for n=3. The agreement is with a
published FINITE-WIDTH estimator. It neither finds a new value of p_c nor proves
an all-width pTL intertwiner. Equation (5) is called a defining sextic here;
no unproved assertion of polynomial minimality is needed.

## 5. Every finite root lies below q; exponential displacement with alternating corrections

At q the five distinct eigenvalues are approximately

    beta_1=lambda_c=0.7957876689642490361,
    beta_2= 0.0362139729093206924,
    beta_3=-0.1005527942557272624,
    mu_+  = 0.1788658499322440484,  multiplicity 2,
    mu_-  =-0.0793337764501209504,  multiplicity 2.

The accompanying Fraction interval calculation proves, rather than infers from
these decimals,

    beta_2>0, beta_3<0, mu_+>0, mu_-<0,
    |beta_3|/mu_+ < 3/5,  |mu_-|/mu_+ < 1/2,
    beta_2 < mu_+ < lambda_c.

Consequently M_(3,m)(q)>0 for every m>=2. For even m every remaining term in
(4) is positive. For odd m, divide by mu_+^m and bound the negative terms by
(3/5)^m+2(1/2)^m < 2. The finite M is strictly increasing by the monotone-rank
argument and has endpoint values -1,+1. Its unique root p_(3,m) therefore
satisfies p_(3,m)<q for every finite length.

Define

    h(p)=log(beta_1(p)/lambda_c(p)),
    rho=mu_+(q)/lambda_c(q)=0.22476579734522079...,
    h'(q)=2.69049781437123949... >0.

The derivative sign is also certified exactly: at q,

    h'(q)=-(q-1)^2 P'(q)/(f_B,L(lambda_c(q),q)*lambda_c(q)).

The numerator derivative is negative and the denominator positive. Analyticity
of the separated eigenvalues near q and a Taylor expansion of (4)/lambda_c^m
give

    p_(3,m)-q = -2 rho^m/(m h'(q)) * [
         1 + (mu_-/mu_+)^m
           + (beta_2/mu_+)^m/2 + (beta_3/mu_+)^m/2
           + O(rho^m)].                               (6)

All ratios on the right are evaluated at q. The O(rho^m) relative error accounts
for evaluation at the displaced root and the nonlinear leading exponential.
This is a fixed-width asymptotic, with a locally uniform spectral gap; it says
nothing about a width-uniform bound or the fixed-aspect L^-4 conjecture.

The beta_3 term supplies the largest relative spectral correction and alternates
with parity. Examples are generated, not fitted:

| m | p_(3,m), diagnostic | actual shift / leading shift |
|---|---|---|
| 3 | 0.5865114551126756357 | 0.84206024 |
| 4 | 0.5883619852843523124 | 1.09362955 |
| 8 | 0.5888800907157023131 | 1.00649648 |
| 12 | 0.5888806988874187696 | 1.00055615 |
| 20 | 0.5888806999178489730 | 1.00000506 |

Exact rational root brackets are stored alongside each diagnostic decimal.

## 6. Why the occupancy-only simplification ends at width four

At width four a proper occupied row may be disconnected. Two short counterexamples
show that merely counting full rows and successive overlaps can err in either direction.
Here integer masks use bit j for column j (columns 0 through 3).

* Rows [1,5,4,5] on the 4x4 torus: all interfaces overlap, no row is full,
  but the ambient rank is zero, not one. The middle masks 5={0,2} do not join
  the incoming and outgoing paths.
* Rows [11,14] on the honest 4x2 torus: again no full row, all interfaces
  overlap, but rank is two, not one. The two complementary paths within the
  rows create transverse winding even though neither row is full. Parallel
  length-two edges are retained.

There is a stronger continuation witness, not just failure of one guessed formula.
Take two open prefixes, ending at the identical frontier 5={0,2}:

    history A: [0,5],      history B: [7,5].

Neither prefix has ambient winding. History B connects the two frontier vertices
via row 7={0,1,2}; history A does not. Append the SAME suffix [13,0], where
13={0,2,3}, and close longitudinally through the empty row:

    [0,5,13,0] has r=0 and D=-1;
    [7,5,13,0] has r=1 and D=0.

In the second history the new route via column 3 closes a transverse cycle
against the old route via column 1. Thus ANY deterministic continuation summary
that identifies these two prefixes cannot preserve this readout. Frontier
occupancy and a current wrapping flag alone are insufficient; connectivity of
frontier vertices is genuine needed memory. This is a scoped lower bound, NOT
an impossibility theorem for local transfer operators or a proof that all full
boundary partitions are necessary/minimal.

The ordinary labelled frontier partition is not the end of the story. An even
stronger pair has the same occupation count, same frontier, same ordinary
partition, and same current winding flag:

    history A': [13,5],     history B': [7,5].

Both contain five of eight sites (identical weight p^5(1-p)^3), both connect
frontier vertices {0,2}, and neither has a transverse cycle. In A' their lifted
connecting path has horizontal displacement -2; in B' it has displacement +2.
Attach the same suffix [13,0]. The new 0-to-2 route has displacement -2, hence

    [13,5,13,0] has r=0, D=-1;
    [ 7,5,13,0] has r=1, D= 0.

The second closed route has net displacement +2-(-2)=4, one circumference;
the first has displacement zero. The script independently traverses the OPEN
prefixes and records all equal ordinary fields and the differing integer lifts,
then verifies the completed torus ranks. This is genuine topological memory
that a bare set partition and a flag recording already completed wrapping lose.
Embedded annular or suitably lifted connectivity states can distinguish it;
this is NOT a no-go for those representations.

The next useful transfer question, if pursued, is therefore a precise width-four
connectivity-and-topology closure, tested against both continuation witnesses.
Plain boundary-partition counts alone do not settle the needed state space.
It is not another occupancy-mask eigenvalue calculation or a large width scan.

## 7. Executed validation and sources

The script derived unnormalized Bernstein integer coefficients by Newton trace
recurrences in site fugacity. An independent lifted-edge traversal checked all
299,584 configurations on 3x2,...,3x6, both the pointwise identity (1) and EVERY
polynomial coefficient. This generalizes the geometric checking technique in
#705 and is independent of the row-transfer formula, not claimed as an unrelated
second research program. The existing 3x3 coefficient vector
[-1,-9,-36,-78,-90,-36,36,36,9,1] is reproduced exactly.

An independent integer 7x7 matrix-power trace at p=1/2 matches the block recurrence
for m=2,...,10. Fraction interval signs certify the q bracket, the five eigenvalue
bands, positivity of h'(q), and the inequalities used in the all-m proof. Seven
small mathematical tests passed locally. No full-repository CI was run for this
new patch; no new production, GPU, or fitted exponent was used.

Files:
- scripts/width3_cylinder_exact.py
- tests/test_width3_cylinder_exact.py
- results/research-control-20260912/width3-cylinder-exact.json

Repository inputs read: PR #705 at fc19cc748527249f0ce1c69d7cd2a88a4879427b,
merged as eb89e9422791d9e3c3a78f0e65d56912b815a7bd; original digital-Alexander
proof; the 3x3 Bernstein rung via #668/#684.

External primary source, PRIMARY_TEXT_READ on 2026-09-12:
Jacobsen, Critical points of Potts and O(N) models from eigenvalue identities in
periodic Temperley-Lieb algebras, arXiv:1507.03027v1, section 6.1 Table 2 and
sections 4,7.1. https://arxiv.org/html/1507.03027v1
The table value and limit order are literature facts; the all-width or continuum
claims not established above remain unclaimed. No broad novelty search was done.
