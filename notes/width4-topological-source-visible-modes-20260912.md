# The intrinsic source exposes twelve modes canceled by the matching function

2026-09-12. Completed exact continuation of #337 and #636. This note uses the
all-parameter P0/P2 trace identities in the preceding delivery, whose source
is draft #708's certified rank automaton. It does not extend the state table,
change width, or identify a continuum source.

## 1. Existing source algebra versus new length complexity

The configuration variable X=r-1 belongs to {-1,0,1}; hence X^3=X. The finite
source partition function is exactly

    Z(p,s)=1+M(p)*sinh(s)+E(p)*(cosh(s)-1),
    M=P_2-P_0, E=P_2+P_0.                              (1)

The three-dimensional algebra of single-configuration observables does NOT
imply a three-state length evolution. Its coefficient functions can contain
many independently evolving spectral modes. This distinction is now exact on
the actual width-four site-percolation model, not just on P398.

Use t=p/(1-p), and unnormalized sequences

    m_m(t)=(1+t)^(4m) M_(4,m)(p),
    e_m(t)=(1+t)^(4m) E_(4,m)(p).

The preceding all-length identities give

    m_m=tr(B15^m)-tr(B5^m),
    e_m=tr(B15^m)+tr(B5^m)-2 tr(B16^m)+4 t^(2m).      (2)

The even source response therefore restores the shared block that cancels in
M. Compared with M, its additional factors are

    shared_3, shared_4a, shared_4b, positive_square,

of degrees 3,4,4,1. Thus exactly twelve additional generic modes become visible.
The names and coefficient polynomials are those of the preceding
`width4-parametric-definition.json`; they are not new fitted factors.

## 2. Exact minimum orders over Q(t)

| sequence | generic minimum constant-coefficient recurrence order in m |
|---|---:|
| m_m: unmarked matching response | 16 |
| e_m: raw second source derivative | 28 |
| (1+t)^(4m) times twice the tilted-mean numerator at exp(s)=2 | 28 |
| 2(1+t)^(4m) Z(p,log 2) | 29 |

The finite source example exp(s)=2 is fixed algebraically, not tuned against
data. The numerator in the third row is 4 P2_un-P0_un. The partition sum in
the last row is

    2(1+t)^(4m)-P0_un+2 P2_un,

which also equals the direct weighted sum with rank-state weights (1,2,4).
Its one additional mode is the total Bernoulli row weight (1+t)^4. Dividing
these sequences by (1+t)^(4m) rescales eigenvalues and not the generic order.

The normalized tilted mean is a RATIO of these sequences:

    M_s = [M*cosh(s)+E*sinh(s)]/[1+M*sinh(s)+E*(cosh(s)-1)]. (3)

No finite order-28 linear recurrence for this ratio is asserted. Likewise
higher derivatives of log Z mix moments nonlinearly; the table is explicitly
about raw source moments, the stated numerator, and the partition sum.

## 3. Why the minimum orders are proved, not guessed

Upper bounds follow from the exact all-length trace decomposition: multiply
one copy of every distinct visible characteristic factor. No new interpolation
or Berlekamp-Massey inference is needed.

For the matching lower bounds, set t=3 (p=3/4), avoiding the known specializations
t=1 and t=2. Direct occupation-weight propagation through the original 509-state
transition table generates the first 60 terms of P0, P2 and all four combinations.
Every term agrees exactly with an independently assembled Newton trace sequence.

For each claimed order d, form the d-by-d physical-tail Hankel minor

    [s_(2+i+j)]_(0<=i,j<d).

Its integer Bareiss determinant is nonzero. Independent Gaussian elimination
in two finite fields verifies the determinant residues:

| sequence | d | determinant mod 1,000,003 | determinant mod 1,000,033 |
|---|---:|---:|---:|
| M | 16 | 950378 | 23606 |
| E | 28 | 303738 | 890451 |
| twice tilted numerator | 28 | 584157 | 627512 |
| twice Z | 29 | 619482 | 991022 |

The full integer minors and determinants are saved. A nonzero specialization
makes the corresponding polynomial Hankel minor nonzero over Q(t), proving
the generic lower bound. The all-length trace formulas supply the upper bound;
checking only 60 moments would not do so by itself.

This proof reuses the certified finite operator rather than constructing
another engine. The module contains no enumeration of the 1,448 geometric
states and no new physical production.

## 4. A real perturbation reactivates canceled modes

At source s, before the common Bernoulli normalization, the mean numerator is

    e^s tr(B15^m)-e^-s tr(B5^m)
      -2 sinh(s) tr(B16^m)+4 sinh(s)t^(2m).            (4)

The B16 coefficient is identically zero at s=0 and nonzero for real s!=0.
At q4 its relative Perron rate is about .773935, slower than M's surviving
leading correction rate .251750. Therefore the faster s=0 cancellation does
not furnish a uniform description of source responses.

This is analogous in logic, but not identical in mechanism, to the P398
result: a description sufficient for one source/readout contract need not
remain sufficient when the contract changes. Here the reactivation comes from
the intrinsic topological reweighting itself, not a new Markov time program.
The finite algebra supplies its exact coefficients, with a physical normalizer.

## 5. Decision and scope

The rank source is a legitimate new probe of this finite model, but it has not
been mapped to either of #275's actual continuum candidates. Low source algebra
dimension, minimum recurrence order, deterministic future-state count, and
continuum module dimension remain different quantities.

No GPU, higher-width scan, broad retrieval job, or fresh Monte Carlo is needed
for this completed result. Record it on existing #337/#636 rather than opening
a new copy of the state-sufficiency or parametric-spectrum task.

Files: `scripts/width4_topological_source_spectrum.py` and
`results/research-control-20260912/width4-topological-source-spectrum.json`.
The input definition remains the preceding delivered artifact; it is NOT
asserted to have reached main. The source automaton is #708 at f782061c.
