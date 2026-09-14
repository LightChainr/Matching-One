# Exact algebraic charge roots at widths 2, 3 and 4

2026-09-14.  Symbolic regression controls for the transparent safe transfer.

The purpose is not threshold estimation.  These small-width polynomials are exact fingerprints of the lifted homology convention and of the NN/complementary-matching charge-sector criterion.

## 1. Setup

Let

\[
K_{4,w}(p),\qquad K_{8,w}(1-p)                               \tag{1.1}
\]

be the two safe frontier kernels.  Their Perron roots are `lambda_4,lambda_8`.  The fixed-width charge root solves

\[
\lambda_4(p)=\lambda_8(1-p).                                 \tag{1.2}
\]

For small `w`, every matrix entry is an integer polynomial in `p`, so one can eliminate `lambda` exactly from the two characteristic equations.  Trivial endpoint/zero-eigenvalue factors are discarded only after symbolic factorization; the physical factor is identified by containing the Perron crossing root in `(0,1)`.

## 2. Width two

The NN safe characteristic polynomial is

\[
\chi_{4,2}(\lambda)
=\lambda^2(\lambda+p^2-1).                                   \tag{2.1}
\]

The complementary matching safe polynomial is

\[
\chi_{8,2}(\lambda)
=(\lambda+p^2-p)
(\lambda^2-p\lambda+p^4-p^3).                               \tag{2.2}
\]

Eliminating `lambda` gives endpoint powers times

\[
\boxed{F_2(p)=2p^3+2p^2-1.}                                  \tag{2.3}
\]

Its unique root in `(0,1)` is

\[
p_2^{ch}=0.5651977173836394\ldots.                           \tag{2.4}
\]

The lifted width-two convention is essential here: the two physically distinct horizontal bonds are retained rather than collapsed into one quotient edge.

## 3. Width three

The NN characteristic polynomial collapses to

\[
\chi_{4,3}(\lambda)
=\lambda^6(\lambda+p^3-1).                                   \tag{3.1}
\]

The matching polynomial factorizes into one squared quadratic block and one cubic Perron block.  Eliminating the NN Perron factor with the matching cubic yields the physical factor

\[
\boxed{
F_3(p)
=p^6-3p^5-5p^4-4p^3+p+1.}                                  \tag{3.2}
\]

The other nontrivial resultant factor has no physical Perron crossing in `(0,1)`.

The unique physical root is

\[
p_3^{ch}=0.5888806999178535\ldots.                           \tag{3.3}
\]

## 4. Width four

The 19-state NN kernel has a Perron eigenvalue contained in a quadratic characteristic factor.  The 19-state complementary matching kernel has its Perron eigenvalue in a quintic factor.  Their exact resultant contains endpoint powers and the degree-17 physical polynomial

\[
\boxed{\begin{aligned}
F_4(p)={}&28p^{17}-98p^{16}-34p^{15}+286p^{14}+122p^{13}
-320p^{12}-362p^{11}\\
&+117p^{10}+377p^9+144p^8-134p^7-174p^6-5p^5\\
&+37p^4+14p^3-10p^2+p+2.
\end{aligned}}                                               \tag{4.1}
\]

It has exactly one real root in `(0,1)`:

\[
p_4^{ch}=0.5914171708531392\ldots.                           \tag{4.2}
\]

## 5. Relation to the semi-infinite graph-polynomial sequence

These three algebraic numbers agree with the first corresponding values of the square-site `n x infinity` eigenvalue-identity sequence reported by Jacobsen.  This is another exact/near-exact bridge between the Bernoulli homology-safe transfer and the graph-polynomial topological sectors.

The present calculation does not imply that the symbolic polynomials `F_w` are the minimal graph-polynomial factors in Jacobsen's variables; changes of local weight variable and elimination can add or remove algebraic factors.  The invariant statement is the physical Perron crossing.

## 6. Why these are useful regression controls

Any proposed rewrite of the safe transfer into

- no-zero-block type-B noncrossing states,
- dilute periodic TL states,
- a more compact canonical annular basis,

should reproduce (2.3), (3.2), and (4.1) exactly after the declared change of variables.

This is stricter than comparing decimal roots: an incorrect multiplicity, seam gain, matching diagonal, or forbidden-winding convention usually changes the symbolic factor immediately.

## 7. Claim boundary

The polynomials above come from exact symbolic characteristic/resultant elimination of the transparent safe kernels.  They are finite-width identities, not asymptotic CFT claims and not new threshold estimates.
