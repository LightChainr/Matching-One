# Jordan visibility is a property of an operator AND its observable

2026-09-12. Direct mathematical review of #724/#731, with exact counterexamples
and a fully positive stochastic control. This does not settle #275 or identify
a continuum module. General perturbation theory and realization theory are not
claimed as new; the deliverable makes the relevant distinctions executable.

## 1. What remains right in the retrieval

Vasseur--Jacobsen--Saleur, arXiv:1206.2312v2, Eq (11), identifies a Jordan cell
through the transformation mixing of specified fields. Its conclusion explicitly
distinguishes that logarithmic observable from logarithms obtained by differentiating
Boltzmann weights. The narrow statement "a thermal m lambda^m factor alone does
not identify the underlying physical Jordan block" remains correct.

But #724 additionally states a converse-like diagnostic (linear split implies
semisimple), treats pointwise semisimplicity as guaranteeing analytic eigenpairs,
and labels m lambda^(m-1) a semisimple signature. Those stronger implications
are false. #731's otherwise useful citation check did not correct them.

Bamieh arXiv:2002.05001v2 section 2 EXPLICITLY assumes that eigenvectors and
eigenvalues are analytic near zero before constructing the series. This is an
extra hypothesis, not a consequence of pointwise diagonalizability. The Kato
body was not read by #724; its table of contents cannot certify a proposed iff.

## 2. Three exact counterexamples

### 2.1 Linear analytic branches through a defective matrix

    A(t) = [[1+t, 1], [0, 1-t]].

The eigenvalues are 1+t and 1-t. At t=0, A is a genuine size-two Jordan block.
A(t) is diagonalizable for t!=0. A nongeneric triangular unfolding can therefore
have linear branches despite defectiveness at the collision. "Linear split =>
semisimple" must be deleted.

### 2.2 Semisimple at zero, diagonalizable pointwise, nonanalytic branches

    B(t) = [[0, t], [t^2, 0]],    char_B(x)=x^2-t^3.

B(0)=0 is semisimple, and B(t) has two distinct eigenvalues for t!=0. Nevertheless
its eigenvalues are +-t^(3/2), which cannot be holomorphic through t=0. An analytic
root would have an integer vanishing order k, but 2k=3.

### 2.3 The problem occurs even for an affine pencil

    C(t) = [[0,t,0], [0,0,t], [t,0,1]],
    char_C(x)=x^2(x-1)-t^3,
    discriminant=-t^3(4+27t^3).

C(0)=diag(0,0,1) is semisimple. For 0<|t|<1/4 all eigenvalues are distinct, so
C(t) is pointwise diagonalizable throughout the neighborhood, including t=0.
The two roots tending to zero are not analytic: x^2(x-1)=t^3 again forces order
3/2. Thus restricting the physical matrix family to A0+t A1 does not rescue
the unsupported analytic-branch claim. The first reduced perturbation on the
degenerate eigenspace is itself defective; individual eigenvectors become badly
conditioned even though each full matrix is diagonalizable.

**Do not overcorrect.** A genuine leading square-root t^(1/2) split under a
holomorphic finite matrix perturbation is incompatible with a semisimple
coalescing eigenvalue: semisimple perturbations move locally by O(t). Higher
fractional powers such as t^(3/2), however, do not imply a defective A(0).
With self-adjointness or appropriate analytic eigenpair assumptions stronger
analytic conclusions are available. Numerical fitted powers are not certificates.

## 3. Positivity and irreducibility do not make traces identify Jordan structure

Let e=(1,1,1)^T, u=(1,-1,0)^T, v=(1,1,-2)^T and define

    P=ee^T/3,    Q=I-P,
    D=uu^T/2-vv^T/6,    N=uv^T/6.

Then N^2=0, N!=0, PN=NP=0, D^2=Q, DN=N, ND=-N. Define two row generators

    G_D(t)=-Q+tD,
    G_J(t)=-Q+tD+N/2,                 |t|<=1/4.                  (1)

Both have zero row and column sums. All off-diagonal entries are at least
1/6 for G_D and 1/12 for G_J over the entire interval. Thus BOTH are irreducible
continuous-time Markov generators with the SAME uniform stationary law. G_D is
symmetric/reversible; G_J is not asserted reversible.

On span(u,v), their matrices are diag(-1+t,-1-t) and

    [[-1+t, 1/2], [0, -1-t]].

At t=0, G_D is semisimple and G_J has a nontrivial Jordan block at -1. Their
minimal polynomials are z(z+1) and z(z+1)^2. Their characteristic polynomials
are identical for every t:

    z[(z+1)^2-t^2].                                           (2)

For T_D=I+G_D/2 and T_J=I+G_J/2 every entry is strictly positive and both row
and column sums are one. The common ordinary trace functions are

    tr(T(t)^m)=1+[(1+t)/2]^m+[(1-t)/2]^m,
    tr(exp(sG(t)))=1+2 exp(-s) cosh(st).                         (3)

Therefore ALL ordinary traces, ALL lengths, and ALL t derivatives agree in
these two families. Any common deterministic normalization, moving-root rule,
or nonlinear statistic formed solely from those identical trace coordinates
also agrees. More precision or more ordinary thermal derivatives cannot
separate this pair.

This is an abstract Markov/stochastic-matrix control. It is NOT claimed to be
square-site percolation, a candidate LCFT, or the actual source map of #275.
It shows that positivity, irreducibility and a common stationary normalization
are not sufficient extra assumptions to make the proposed trace diagnostic valid.

## 4. One ordinary probability readout distinguishes the same pair

Keep the source delta_state_1 and readout 1{state=3}. At t=0 the discrete-time
responses are exactly

    y_D(m)=(1-2^(-m))/3,
    y_J(m)=(1-2^(-m))/3 - (m/6)2^(-m).                         (4)

These are actual transition probabilities, not arbitrary signed test vectors.
The 2x2 Hankel determinant of y_D, starting at m=0, is -1/36; the 3x3 determinant
of y_J is -1/6912. Their exact minimal scalar recurrences have polynomials

    (x-1)(x-1/2),       (x-1)(x-1/2)^2.

Thus the **fixed physical operator** has a visible repeated pole in this specified
response, even though every ordinary trace in (3) is blind to it. In continuous
time the response difference is -s exp(-s)/6. A stationary uniform source instead
is blind to the difference for every readout, because the invariant source kills
N. The source matters as much as the eigenvalues.

## 5. Exact general visibility condition

For a finite matrix A, let P_lambda be its generalized spectral projectors and
N_lambda=(A-lambda I)P_lambda. For a specified row source C and column readout B,

    C exp(sA) B
      = sum_lambda exp(s lambda)
          sum_{k>=0} s^k/k! C N_lambda^k P_lambda B.            (5)

A Jordan contribution at a given eigenvalue is visible precisely when at least
one coefficient with k>=1 is nonzero (after combining all blocks at that same
lambda). This is elementary finite-dimensional functional calculus, not a new
LCFT criterion. A repeated pole in a minimal exact fixed-operator response proves
Jordan structure in that minimal realization; noisy finite approximations do not
supply such an exact conclusion without separation/error assumptions.

For an ordinary trace, every nilpotent term has trace zero. More generally, if
C is an invariant/commuting mark, [C,A]=0, then C commutes with P_lambda and
N_lambda. For k>=1, C N_lambda^k P_lambda is nilpotent, so

    tr[C N_lambda^k P_lambda]=0.                              (6)

Thus ordinary traces and invariant-sector traces do not directly read Jordan
nilpotents. A noncommuting mark or a specified off-diagonal matrix element may.
This is NOT a blanket statement about all marked traces or words containing
multiple noncommuting operators.

For #275 the practical preliminary question is: for its ACTUAL source/closure,
are the coefficients L(N_lambda^k P_lambda) represented and nonzero? The answer
cannot be inferred from the word "trace", "spin 4", or an apparent m factor.
The existing two-candidate original-U gate stays in place; no source change or
third model is authorized by this note.

## 6. Thermal jets are a different operator

For an analytic physical matrix A(t), the block matrix

    J_1(t0) = [[A(t0), A'(t0)], [0, A(t0)]]

encodes the first parameter derivative: its upper-right block in J_1^m is
partial_t A(t)^m at t0, and similarly for the exponential. This follows by
multiplying block upper-triangular matrices. Higher jets act on a truncated
polynomial ring, with derivative blocks divided by factorials.

Even A(t)=diag(1+t,1-t), which is semisimple, gives a nontrivial nilpotent part
in J_1(0). A repeated pole or polynomial length factor in the JET realization
therefore does not prove a Jordan block in the original physical A(t0).
One must identify which operator and which observable is being reconstructed.

## 7. Executed checks and citation boundaries

The symbolic program verifies both characteristic/minimal polynomials, exact
endpoint positivity of affine generator entries (hence positivity throughout),
row/column sums, the two small Hankel determinants, three perturbation examples,
and nine jet-power identities. It compares trace polynomials and their derivatives
through order four for lengths 0..12 as regressions. The all-length statements
follow from the displayed exact triangular blocks, not those 65 regressions.

A separate standard-library Fraction program reconstructs the physical matrices
independently, verifies 10 generator instances, 105 trace moments and both
Hankel determinants by rational Gaussian elimination. It imports neither SymPy
nor the producing script.

Primary texts read this round: Bamieh 2002.05001v2 section 2 and assumptions;
Vasseur--Jacobsen--Saleur 1206.2312v2 Eq (11) and conclusion. Qian--Chu--Tan,
SIAM J. Matrix Anal. Appl., DOI 10.1137/15M1053050, is ABSTRACT_ONLY and is not
used as the proof of any counterexample. No global novelty search was performed.
