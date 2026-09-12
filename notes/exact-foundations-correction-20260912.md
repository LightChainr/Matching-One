# Exact foundations: corrected topology, shape, and hidden criticality

Date: 2026-09-12. Base read: `main@8b5bc840e4c499e3bfa31f61bfd6f2a713b2e1d5`.
This note supersedes the affected arguments, not the historical measurements, in
PRs #675, #690, #693, #666 and #688. No percolation production, threshold
estimation, or large census is performed. New proofs below are mathematical
consequences/constructions, not claims of literature novelty.

## 1. All-size wrapping label map: no directional conjecture is needed

Use exactly the scope of `notes/digital-alexander-duality-proof.md`: an honest
periodic square-cell torus, with the NN black graph and complementary NN+NNN
white graph. In a common homology basis let their rational ambient images be
A,C in H_1(T^2;Q). The existing embedded-reduction and subsurface lemma gives
C=A^perp for the intersection pairing omega((a,b),(x,y))=ay-bx.

**Elementary correction.** If A=span(a,b) is a nonzero line, then

    A^perp = {(x,y): ay-bx=0} = A.

In particular span(1,0)^perp=span(1,0), NOT span(0,1). PR #690 confused
symplectic orthogonality with Euclidean orthogonality; the same mistake entered
its reviews and #693. No extra face-parity or 90-degree argument is required.

**Theorem (allowed joint labels).** In the stated scope, the only possible pairs
of the rank-refined labels none, x, y, spiral, cross are

    (none,cross), (cross,none), (x,x), (y,y), (spiral,spiral).

Proof. The rank-sum theorem permits (0,2),(2,0),(1,1) only. At rank one,
C=A, so both sides have the same rational winding line. This forces the same
single-axis label or the same both-coordinate spiral label. Distinct components
of either embedded graph have zero intersection pairing. They cannot contribute
two independent homology lines on a torus; if the total image has rank two, one
component therefore carries it. The white embedded reduction preserves both
connectivity and ambient classes by the facewise replacements in the existing
proof. Thus rank two is a cross, not two disjoint transverse clusters. This
also proves the both-two label is empty. QED.

This is an **allowed-support** theorem, not a statement that every allowed cell
has positive mass at every size. It equates rational lines and labels, not the
integer subgroup, cluster counts, or a literal colour-flip operation. Degenerate
periodic quotients outside the original proof remain outside this result.

Consequences: the black-spiral and white-spiral events are the SAME configuration
set at every size in scope; their cancellation is pointwise, not an involution
between distinct configurations. Directional equality holds in the (1,1) case,
not for all configurations (all-black already refutes that unqualified wording).
For any probability measure on these configurations,

    D = 1{r_b>0}-1{r_w>0} = r_b-1,
    M = E[D] = P_2-P_0,
    F := E[r_b]/2 = (1+M)/2.

These need neither Bernoulli independence nor a census. Under the Bernoulli
family F is strictly increasing, so its median parameter is the unique matching
root: p*_N=Q_N(1/2). This says nothing about rates or a closed form for p_c.

## 2. What the finite-operator objection actually establishes

PR #675's probability/eigenvalue distinction is important but is not a general
impossibility theorem. For example, on the space indexed by configurations,
let W_p be diagonal with the Bernoulli weights and P_j the diagonal projector
onto r_b=j. Then M=tr((P_2-P_0)W_p) is an exact single-operator representation.
This is an exponentially large tautological representation, not an efficient
row transfer, not pTL, and not a leading-eigenvalue identity.

The useful unresolved question is therefore restricted: specify local weighted
propagation, closure, the sector maps and the cylinder limit before comparing
with Jacobsen's criterion. A finite probability difference need not equal a
leading-eigenvalue difference. That does NOT imply that no operator can encode
it. The black-only identity above also defeats a blanket claim that both graph
connectivities must always be stored. #636 receives no automatic large-compute
release from this correction; #681 should own the precise cylinder bridge.

## 3. Shape symmetry is not reflection about the raw parameter 1/2

PR #666 proposes C3: sup_{p in K}|M_N(p)+M_N(1-p)| -> 0 for every compact
K in (0,1), and identifies it with normalized shape symmetry. Under the
repository's location conclusion F_N(p)->0 below c and ->1 above c, this C3
forces c=1/2: if c>1/2, M_N(1/2)->-1; if c<1/2 it tends to +1. Either way,
for c!=1/2 the absolute sum at 1/2 tends to 2, not zero. Thus C3 cannot be
used as the shape hypothesis for a differently centred transition.

With symmetric anchors a,1-a and W_N=Q_N(1-a)-Q_N(a)>0, the correct diagnostic is

    A_N(u) = [Q_N(u)+Q_N(1-u)-Q_N(a)-Q_N(1-a)]/W_N
           = Z_N(u)+Z_N(1-u)-1.

Its convergence to zero is a shape question independent of the limiting centre.
It remains OPEN for the relevant percolation families. Self-duality constrains
shape but does not select a unique symmetric profile. Two exact tiny sizes,
a bond control and one site production are not a within-model asymptotic ladder.

There is a second elementary correction to #622/#624 and their descendants.
For an increasing affine map phi(p)=alpha*p+beta, ordinary pushforward gives
F_phi=F composed with phi^{-1}, hence Q_phi=phi composed with Q. These actions
DO intertwine. On a fixed [0,1] support one must track the transformed support;
forcing both endpoints fixed leaves only the identity. A nonlinear kink is
not an affine map and cannot distinguish two alleged Aff(1) actions.
The valid location-without-shape constructions are unaffected by this correction.

## 4. An irreducible exact no-go replacing the direct-sum overclaim

The old `notes/bounded-task-rank-threshold-no-go-20260907.md` conflates Markov
irreducibility and reducibility of a real linear representation. Already
G=[[-1,1],[1,-1]] is an irreducible generator commuting with the nontrivial
state swap R. Its even and odd eigenspaces are not disconnected Markov classes.
Also, a direct sum of irreducible generators is not irreducible as a chain.

**Theorem (product-chain repair).** For c in {1/2,1/3}, p in [0,1], L>=2,
let H_L(p-c) be the reflecting nearest-neighbour walk with right rate
exp(p-c), left rate exp(c-p). Set V=[[-2,2],[2,-2]] and

    G_{L,c}=V tensor I_L + I_2 tensor H_L(p-c).

Every finite chain is connected, irreducible, reversible, entrywise analytic,
and local on a two-layer ladder (maximum degree three, uniformly bounded
positive rates). Let e0=(1,0)^T, B=e0 tensor 1_L, and
C=e0^T tensor (1_L^T/L), common to both families. Since H_L*1_L=0,

    C exp(t G_{L,c}) B = e0^T exp(t V) e0 = (1+exp(-4t))/2

for all L,p,t. The exact task order INCLUDING the stationary constant is two;
the decaying contrast alone has order one. Yet the full-chain gap is

    gamma_{L,c}(p)=2*cosh(p-c)-2*cos(pi/L),
    gamma_{infinity,c}(p)=4*sinh((p-c)/2)^2.

To verify the gap formula, similarity by the square root of the birth-death
stationary weights makes each hidden off-diagonal equal sqrt(ab)=1; its
nonzero eigenvalues are -(a+b)+2*sqrt(ab)*cos(k*pi/L), k=1,...,L-1.
The product-chain rates are sums with 0 and 4. For the stated c,p,L,
the smallest hidden rate is at most 2*cosh(2/3)<4, so it is also the full gap.
Only p=c closes in the limit. QED.

Thus exact task equality can coexist with different gap-closing locations
EVEN FOR OVERALL IRREDUCIBLE CHAINS. This is an abstract spectral-threshold
counterexample, not percolation, not a universality claim, and not robustness
to arbitrary interactions that destroy the product structure. The construction
is presented without a novelty claim. The accompanying mpmath check is a
numerical diagnostic of the formula, not the proof.

A related correction: a mode's response residue is (C v)(w^T B). It vanishes
if either factor vanishes, not only when both do. In systems with repeated
modes/couplings, use the minimal realization/Kalman decomposition rather than
an unsupported if-and-only-if assertion about an arbitrary invariant subspace.

## 5. The nonnegative-rank separation survives with logarithmic growth

`notes/positive-vs-signed-slack-separation-20260906.md` and #688 state that a
regular n-gon's slack matrix has nonnegative rank n. This is false. Its ordinary
rank is three and its nonnegative rank is Theta(log n). The unbounded-versus-
bounded separation remains valid; the claimed linear growth does not.

PRIMARY_TEXT_READ, 2026-09-12: Fiorini, Rothvoss, Tiwary, *Extended formulations
for polygons*, arXiv:1107.0371v2, sections 1-3 (HTML):
https://arxiv.org/html/1107.0371v2
Theorem 1 identifies extension complexity with slack nonnegative rank;
Theorem 2 gives O(log n) for regular polygons; section 1 gives the matching
Omega(log n) bound. No polygon-to-native-cut-network realization is proved here.

## 6. Consequence for further work

Do not commission larger wrapping censuses to prove the directional label map;
it follows already. Do not test #666's raw-p C3 as though it were shape symmetry.
Do not infer a universal one-operator no-go from #675, or a percolation threshold
from the product-chain construction. Retain the valid exact data, location
lemma, finite non-scalar failures and original-U identifiability problem.

The next empirical target should be the corrected within-model shape diagnostic
and full-vector residual on existing histograms, not another fitted exponent.
P3 needs an explicit covariance-support policy before singular scores are used;
its missing N580 covariance should first be recovered from committed `_deleted`
arrays, not bought by another production.

Validation: `python -m unittest discover -s tests -p 'test_research_control_exact.py'`.
Five checks pass locally, including 54 product-chain response checks at 40 digits.
No full repository test suite was run in this connector-only checkout.
