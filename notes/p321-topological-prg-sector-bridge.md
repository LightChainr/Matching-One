# P321: the homology root and the periodic-TL crossing are one aspect-ratio family

Status: exact sector dictionary plus an `E3` scaling consequence.  The common
field identity remains a falsifiable mechanism hypothesis.

## 1. Exact finite-volume bridge

For a periodic basis of width `n` and length `m`, let `Z_0D`, `Z_1D`, and
`Z_2D` be the FK/site-percolation state sums with homology-image rank zero,
one, and two.  Jacobsen's graph polynomial is

\[
P_{n,m}(Q,p)=Z_{2D}(Q,p)-QZ_{0D}(Q,p).
\]

At `Q=1`, after the harmless positive normalization that turns the state sum
into a probability measure,

\[
P_{n,m}(1,p)=P_2(p)-P_0(p).
\]

Digital Alexander duality gives independently

\[
M_{n,m}(p)=\mathbb E_p[r_B]-1=P_2(p)-P_0(p).
\]

Therefore the matching/homology polynomial and the `Q=1` critical polynomial
are not merely analogous continuum observables.  They are the **same finite
topological projector**, with `Z_1D` absent because its coefficient in
`E[r_B]-1` is zero.

For fixed `n`, Jacobsen's periodic Temperley--Lieb transfer matrix decomposes
the zero-string block into `T_open` and `T_closed`.  Its Perron roots obey

\[
Z_{2D}=a_o\Lambda_o^m(1+o(1)),\qquad
Z_{0D}=a_c\Lambda_c^m(1+o(1)).
\]

Thus

\[
\lim_{m\to\infty}\{P_{n,m}(Q,p)=0\}
\quad\Longleftrightarrow\quad
\Lambda_o(n,p)=\Lambda_c(n,p).
\]

The sector meanings are explicit:

- `open`: an FK cluster propagates along the infinite cylinder;
- `closed`: a dual FK cluster propagates along the infinite cylinder;
- both are zero-string sectors because one-dimensional torus contributions
  were removed by assigning winding loops weight zero.

So `p_n^TL` is the infinite-aspect-ratio endpoint of the same family of roots
`p_{n,m}^H`, whereas the ordinary matching root uses finite aspect ratio
(usually `m/n=1`).  No speculative modular transform is needed to establish
this common origin.

### A useful finite-length correction

Before taking `m` infinite, the root equation reads

\[
m\log(\Lambda_o/\Lambda_c)+\log(a_o/(Qa_c))+o(1)=0.
\]

Hence fixed-`n` finite-length roots approach the TL root generically as
`1/m`, with a coefficient determined by the channel amplitudes and the
derivative of the log-eigenvalue ratio.  This separates aspect-ratio
convergence from the `n^-4` transverse finite-size drift.

## 2. The Q=1 spin-language obstruction

For integer `Q>1`, Jacobsen also writes the crossing as

\[
\Lambda_0^{(1)}=\Lambda_1,
\]

between a nontrivial global-spin representation in the untwisted sector and a
cyclically twisted sector.  This is an electric/magnetic equality for the same
magnetic scaling dimension.

At literal `Q=1`, however, the one-state Potts spin representation has neither
a nontrivial global-spin irrep nor a nonzero twist label; moreover the spin
formula used to extract the graph polynomial carries a prefactor `1-1/Q`.
Consequently the spin-sector dictionary must be understood by analytic
continuation `Q -> 1`.  The FK/TL `open` versus `closed` statement is the
nonsingular percolation formulation and is the correct bridge to
`P_2-P_0`.

This is a precise obstruction to saying that two literal `Q=1` spin
eigenvectors have been identified.  It is not an obstruction to the exact
topological-projector identity above.

## 3. One aspect-ratio scaling function gives the root exponent

Let `rho=m/n`.  Suppose the first matching-odd perturbation of the continuum
projector has dimension `x`, coupling `g_u`, and aspect-ratio response
`F_u(rho)`.  Let the thermal field have eigenvalue `y_t=3/4`, metric factor
`g_t`, and response `F_t(rho)`.  The dimensionless contrast has the expansion

\[
M_n(p_c+t;\rho)=
g_u n^{2-x}F_u(\rho)+g_t t n^{y_t}F_t(\rho)+\cdots .
\]

The implicit root is therefore

\[
p_n^H(\rho)-p_c=
-\frac{g_u}{g_t}\frac{F_u(\rho)}{F_t(\rho)}
n^{-(x-2+y_t)}+\cdots .
\]

For the thermal level-four spin-four candidate `x=21/4`,

\[
x-2+y_t=\frac{21}{4}-2+\frac34=4.
\]

The same calculation can be made with the cylinder free-energy difference:
its irrelevant splitting is `n^-x`, while its thermal derivative is
`n^(-2+y_t)`, again giving `n^-4`.

This sharpens the central conjecture:

> The square-torus homology root and the TL crossing are exact aspect-ratio
> versions of one projector.  They share the `n^-4` law if the same
> matching-odd `x=21/4` perturbation has nonzero response at both `rho=1` and
> `rho=infinity`.

The candidate is a Potts/internal singlet but a **spatial spin-four** field;
calling it simply an "irrelevant singlet" discards the orientation information
that made it identifiable.

If numerator and thermal response each have relative `n^-2` corrections,
their quotient automatically produces the next root term `n^-6`.  This does
not reuse the rejected scalar-`q=2` hypothesis from a different derivative
observable.

## 4. New amplitude-ratio target

The common microscopic coupling and thermal metric cancel between two aspect
ratios:

\[
\frac{C_H(\rho=1)}{C_{TL}}
=
\frac{F_u(1)/F_t(1)}{F_u(\infty)/F_t(\infty)}.
\]

This is the parameter-free object that can identify a common field more
strongly than a shared exponent.  It can be obtained either from continuum
torus/cylinder matrix elements or by a controlled sequence of rectangular
matching roots at several `rho`, followed by a frozen `rho -> infinity`
extrapolation.

The decisive failure modes are now clean:

1. `F_u(infinity)=0`: the TL `n^-4` term comes from another operator despite
   the exact common projector;
2. the exponent is four at both shapes but the frozen amplitude-ratio curve
   fails: exponent degeneracy, not a common field;
3. the ratio closes across aspect ratios: direct support for the shared
   thermal-Q4/spin-four mechanism.

## 5. Parameter-free check of the published TL sequence

`scripts/p321_fixed4_sequence_score.py` evaluates the dyadic statistic

\[
R_n=\frac{p_{4n}-p_{2n}}{p_{2n}-p_n},
\qquad R_n\to 2^{-4}=1/16,
\]

without inserting or fitting `p_c`.  On the committed Jacobsen widths:

| n | `R_n` | effective exponent |
|---:|---:|---:|
| 2 | 0.0478877 | 4.3842 |
| 3 | 0.0618990 | 4.01394 |
| 4 | 0.0548143 | 4.18930 |
| 5 | 0.0579138 | 4.10995 |

The `n=1` row is strongly pre-asymptotic.  The late rolling two-width
prediction residual (test widths 14--21) has median absolute value
`1.26e-9` for fixed power four, versus `5.44e-8`, `4.88e-8`, and `9.58e-8`
for fixed powers three, five, and six.  These alternatives are controls, not a
free exponent fit.

This confirms that the old cylinder data contain a strong parameter-free
`Delta=4` signature.  It does not by itself identify the field; the
aspect-ratio amplitude test above does.

## 6. Next calculation

The highest-information acquisition is no longer another cylinder width.
Generate roots of the same `P_2-P_0` projector at a small frozen rectangle of
aspect ratios, for example `rho in {1, 3/2, 2, 3, 4}`, using shared rank-threshold
streams.  Fit only

\[
p_n(\rho)=p_c+C(\rho)n^{-4}+D(\rho)n^{-6}
\]

and compare the normalized `C(rho)` curve with the TL endpoint.  Preserve
cross-`rho` covariance.  This directly estimates the open/closed channel
interpolation and creates the amplitude-ratio fingerprint.

### Equal-area implementation already supported by the repository

The general integer-period threshold-rank engine requires paired matrices of
equal determinant.  This is an advantage: compare several rectangles at the
same area, reusing one square reference and the exact same counter permutation.
`scripts/p321_equal_area_rectangle_design.py` freezes three scaled copies of

\[
(12,12),\ (9,16),\ (8,18),\ (6,24),\ (4,36),
\]

with aspect ratios `1`, `16/9`, `9/4`, `4`, and diagnostic `9`.  At scales
one, two, and three the common areas are `N=144,576,1296`.  The root law is
`N^-2` with a first `N^-3` correction.

Run the square as the first matrix against each rectangle using identical
seed, replica interval, and batch boundaries.  The repeated square histograms
must then be byte-identical.  The aligned delete-one-batch roots provide the
full cross-aspect covariance without adding a new Monte Carlo engine.  The
`rho=9` row is only an endpoint diagnostic, not a claim that the TL limit has
already been reached.

## Sources used

- Jacobsen, *Critical points of Potts and O(N) models from eigenvalue
  identities in periodic Temperley--Lieb algebras*,
  https://arxiv.org/abs/1507.03027
- Arguin, *Homology of Fortuin--Kasteleyn clusters of Potts models on the
  torus*, https://arxiv.org/abs/hep-th/0111193
