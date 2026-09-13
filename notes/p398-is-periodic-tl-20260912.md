# P398 is exactly a periodic identified-connectivity TL stochastic chain

2026-09-12. Correction to #715 §Q3a and its review comment, with a constructive
all-width identification. No new process is defined. No generic realization
or integrability theorem is claimed as new. #709's specific input/output
certificates remain valid and acquire a precise existing-model interpretation.

## 1. The literature comparison missed a change of representation

#715 argues that TL acts on noncrossing matchings and is join-only, whereas
P398 acts on noncrossing partitions and contains detach moves. This does NOT
separate the processes. There is a standard fattening/medial bijection from
noncrossing partitions of w cyclic points to noncrossing perfect matchings
of 2w cyclic endpoints. Under that bijection the two halves of the TL generators
are exactly the point-detach and adjacent-join maps.

The relevant published model is NOT the full faithful TL representation, not
a through-line sector sum, and not the cylinder's distinct-connectivity (DC)
representation. It is the periodic IDENTIFIED-CONNECTIVITY (IC), zero-defect,
loop-weight-one representation on disk link patterns, of dimension Catalan(w).
Pearce--Rittenberg--de Gier--Nienhuis (2002), arXiv:math-ph/0209017v2,
Eq. (2.15) together with the IC quotient discussion around (2.17), names it.
Cantini--Sportiello, arXiv:1003.3376v1, §2.2 Eq. (4), gives the exact reconnection
map without requiring interpretation of a drawing.

The source's statement that lines may join without the reverse concerns the
algebra's line/defect filtration; it is not a statement that the induced map
on w-point set partitions can never detach a point. There are 2w local TL maps,
not w. Losing this factor of two hides the detach half.

## 2. Explicit bijection and proof at every width

Use zero-based indices i in Z/wZ. Replace point i by consecutive endpoints
L_i=2i and R_i=2i+1. For each block B={i_1,...,i_k}, ordered cyclically, connect

    R_(i_j) to L_(i_(j+1)),  j modulo k.                            (1)

For a singleton this is the adjacent pair (2i,2i+1). Call the resulting matching
Phi(pi). It is noncrossing: the curves can be drawn as the oriented boundary of
a small regular neighbourhood of the disjoint noncrossing block trees.

Conversely a noncrossing perfect matching of 2w cyclic endpoints pairs opposite
parities (each arc contains an even number of endpoints on either side).
The map sigma(i)=mate(2i+1)/2 is a permutation. Its cycles give the original
blocks. In the disk picture these are precisely the components obtained by
contracting each L_i,R_i interval; noncrossing guarantees their cyclic successor
order and prevents interlacing blocks. This inverse recovers (1), proving a
bijection, not merely equality of Catalan counts.

Let e_a on link patterns join adjacent endpoints a,b=a+1. If they are already
paired, it does nothing. Otherwise remove (a,c),(b,d) and replace them by
(a,b),(c,d), all indices modulo 2w. Then

    Phi(detach_i pi) = e_(2i) Phi(pi),
    Phi(join_(i,i+1) pi) = e_(2i+1) Phi(pi).                        (2)

For the first identity, the reconnection isolates L_i,R_i and reconnects the
predecessor and successor of i in its old block. For the second, it splices
the cyclic successor lists of the two neighbouring blocks. If the two points
were already in one block, they are consecutive within it, so e does nothing.
The argument includes i=w-1 at the cyclic seam. Contractible-loop weight is 1,
which is essential for no extra scalar coefficient in the no-op case.

## 3. Generator, intervention, and boundary convention

The repository Generator acts on functions, with row convention

    G_eta f(pi)=sum_i (1+eta)[f(join_i pi)-f(pi)]
                       +(1-eta)[f(detach_i pi)-f(pi)].              (3)

Thus under Phi it is exactly the transpose of the standard column-state TL
intensity generator

    L_eta=sum_i [(1-eta)(E_(2i)-I)+(1+eta)(E_(2i+1)-I)].            (4)

At eta=0 this is the homogeneous periodic dense O(1)/TL stochastic process on
2w endpoints, unit clock per generator. A discrete-time convention averaging
by 2w only rescales time. For eta!=0 it is an alternating-rate TL chain. No
integrability claim for that whole staggered family is inferred here.
The localized #709 pulse also maps exactly:

    H=J_0-J_(w-2) <-> E_1-E_(2w-3),                               (5)

with the same function/state transpose convention. Being outside span{J,D}
does not put this perturbation outside the TL operator algebra.

P398's IC disk quotient forgets front/back annular path distinctions. It must
NOT be substituted for #708's homology-preserving lifted torus closure. The
P398 process remains different from microscopic square-site percolation, and
its parameter eta is not thereby identified with occupation probability p.

## 4. An exact half-step complement and new dictionary identities

Let rho rotate all 2w endpoints by +1 and define K=Phi^-1 rho Phi. This is a
permutation of partition states, not generally an involution. If sigma is the
cyclic block-successor permutation and c(i)=i+1, then

    sigma_K=c sigma^-1,
    K^2 = one-site cyclic rotation,
    K^(2w)=id.                                                    (6)

It is a Kreweras-type complement, with the indexing convention fixed by (6).
The old geometric reflection is a different symmetry; under Phi it reflects
endpoints a->2w-1-a.

Equation (2) immediately gives

    K detach_i = join_i K,
    K join_i = detach_(i+1) K.                                    (7)

Consequently the rate-pencils are conjugate:

    G_eta[pi,pi'] = G_(-eta)[K pi,K pi'].                          (8)

They are exactly isospectral at eta and -eta. This holds for every width, and
for an arbitrary time-dependent eta(t), with the corresponding sign-reversed
schedule after conjugacy. It is not a fit from a finite parameter ladder.

The readouts are transformed too, not silently held fixed. Write b(pi) for
number of blocks, s_i(pi) for the indicator that i is a singleton, and
c_i(pi)=1{pi_i=pi_(i+1)}. Then

    b(K pi)=w+1-b(pi),
    s_(i+1)(K pi)=c_i(pi),
    c_i(K pi)=s_i(pi).                                            (9)

The last two follow by tracking the corresponding adjacent matched pair. For
the first, the planar bipartite incidence map between the blocks and the
complementary blocks has w edges and one face; Euler's relation gives
cycles(sigma)+cycles(c sigma^-1)=w+1. Equivalently, draw each noncrossing block
as a tree in the disk: the complementary regions containing the interleaved
points are exactly the cycles of c sigma^-1, and each added tree edge increases
the number of such regions by one. There are w-b(pi) tree edges, so there are
w-b(pi)+1 complementary blocks.

In particular the original wrap(pi)=1{pi_0=pi_(w-1)} becomes singleton_0 after K.
The declared three-readout dictionary is NOT K-invariant, so this does not
license a larger observable-admissible quotient than the one already certified.

For any source mu,

    E^eta_mu[b(X_t)] + E^(-eta)_(K_*mu)[b(X_t)] = w+1.              (10)

Uniform source is K-invariant. Hence at eta=0 its expected block count is
(w+1)/2 at EVERY time, even though the uniform measure need not be stationary.
All-singleton and single-block source states are interchanged. The wrapped-pair
source generally maps to a state outside that short source list; retaining the
source transformation is essential. At stationarity for |eta|<1,

    pi_eta(pi)=pi_(-eta)(K pi),
    E_(pi_eta)b + E_(pi_-eta)b = w+1.                              (11)

## 5. A published all-width stationary description is now applicable

Cantini--Sportiello's proof of the Razumov--Stroganov correspondence applies
exactly to the link-pattern operator in (4) at eta=0. In their notation,
§2.4 Eqs. (22)--(24), H_w=sum_(a=1)^(2w) e_a and the proved equality is
H_w |s_w>=2w |s_w>, where |s_w> counts square fully-packed loops refined by
boundary link pattern. Thus the P398 stationary law is

    pi_0(pi) = FPL_w(Phi(pi))/A_w,
    A_w = product_(j=0)^(w-1) (3j+1)!/(w+j)!.                     (12)

This is an import of their theorem, not a new proof or discovery. It describes
the pushforward of a static FPL ensemble; it does not construct a Markov dynamics
on FPL configurations conjugate to P398. No corresponding formula for arbitrary
eta is claimed.

The exact independent stationary solves at widths 2,3,4,5 give primitive integer
normalizations 2,7,42,429, respectively. At w=4 the 14 weights take values 1,3,7;
at w=5 the 42 weights take values 1,4,6,14,17,42. These are consistency checks
of the identified model and theorem scope, not independent evidence for (12).
No large FPL enumeration or width extension is required to use the theorem.

## 6. Executed verification and implications

`scripts/p398_tl_fattening.py` constructs every matching at w=2..8 and verifies
its inverse partition map, all 31,040 join/detach conjugacy equalities, the seam,
all defining TL map relations, the half-step square, and each readout identity.
At w<=6 an independent restricted-growth partition enumeration gives the same
state set. Fraction stationary solves at w<=5 verify (11) also at eta=+/-1/4.
The full G and H matrices at w=4,5 were independently rebuilt from reconnections
and agree entrywise with the unchanged #709 analysis script in the supplied
archive. The previously certified double-pulse claims remain scoped as before;
we do not re-score or extend their rank certificates here.

The statement in #715 that no named process matches P398 should be withdrawn.
Its realization-theory citations and the fact that #709 is a model/observer-
specific certificate remain useful. The positive outcome is not only a wording
correction: it gives an exact literature dictionary, the stationary FPL law,
and extra rate/readout transport constraints for subsequent work.

The generic TL algebra can be nonsemisimple at this loop value, and different
representations have different Jordan structure. Neither the algebra's name
nor the stationary correspondence identifies a continuum field for #275, or
proves Jordan visibility in a particular P398 input/output channel.

Primary sources read:
- Pearce--Rittenberg--de Gier--Nienhuis, J. Phys. A 35 (2002) L661--L668,
  https://arxiv.org/html/math-ph/0209017v2, §2 Eqs. (2.15)--(2.17), §3.
- Cantini--Sportiello, JCTA 118 (2011) 1549--1574,
  https://arxiv.org/html/1003.3376v1, §§2.2--2.4 and proof statement in §3.
  Source's loop-parameter sign convention differs from the 2002 notation;
  both specialize the deleted-loop multiplier to 1.
- Actual repository maps: main@eb89e942, scripts/planar_state_operations.py;
  PR #709 head 83e011a7, existing G and H; #715 head 308f66c7 and its review.
