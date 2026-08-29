# Exact neutral-area covector and the bounded mode-front test

Status: exact finite-`N` identity, exact `N=10` oracle, and a retrospective
bounded-order diagnostic on committed threshold histograms.

## Exact covector

For one threshold pair define

\[
m_k=\begin{cases}-1&k<K_-\\0&K_-\le k<K_+\\+1&k\ge K_+,\end{cases}
\qquad U_k=1-m_k^2.
\]

Marginal threshold histograms suffice because

\[
q_k=\mathbb E U_k
=\Pr(K_-\le k)-\Pr(K_+\le k).
\]

Let `c_r=E_{p0}[q_K H_r(K)]` use the repository's positive-score
orthonormal `Bin(N,p0)` Krawtchouk basis.  The exact finite-size generating
function and termwise integration give

\[
\boxed{
\frac{\mathbb E(K_+-K_-)}{N+1}
=\sum_{r=0}^N \ell_r(N,p_0)c_r,
}
\]

where

\[
\ell_r(N,p_0)=
\sqrt{\binom Nr}\,
\frac{(1-p_0)^{r+1}-(-p_0)^{r+1}}
{(r+1)[p_0(1-p_0)]^{r/2}}.
\]

Thus canonical neutral width is a known left covector of the **complete**
score vector.  It is not, however, automatically a well-conditioned
functional of a low-order local jet.

At half filling the covector reduces exactly to

\[
\ell_{2j+1}=0,\qquad
\ell_{2j}=\frac{\sqrt{\binom N{2j}}}{2j+1}.
\]

The exhaustive self-matching `N=10` cross-channel oracle has neutral area
`1/7`, mean gap `11/7`, zero odd modes, and successive nonzero area
contributions

\[
\frac5{16},\ -\frac5{16},\ \frac3{16},\ -\frac5{112}
\quad (r=0,2,4,6),
\]

which sum exactly to `1/7`.

## Frozen bounded diagnostic

For the size-local intrinsic center define

\[
T_R=\frac{\mathbb E G}{N+1}-\sum_{r=0}^R\ell_r c_r.
\]

The retrospective diagnostic froze `R<=16` and relative tolerances `0.05`
and `0.10`.  The proposed Hermite turning-point picture suggested that the
first adequate order might scale as `N^(1/4)`.

That proposal does **not** survive the committed data at bounded order.  For
every analyzed size `N=65,85,130,170,185,265,325,425`, neither tolerance is
reached by `R=16`.  Even the best partial sum over `R<=16` has relative tail
larger than one; high-order partial sums rapidly amplify cancellation and
sampling error.

The reason is structural.  The canonical area integrates over all
`p in [0,1]`, while the intrinsic Krawtchouk basis is localized near
`p0 approximately 0.592746`.  The exact covector therefore performs a severe
global extrapolation from local score coordinates.  Its coefficients grow
rapidly and the full degree-`N` cancellation is essential.

## Consequence for transfer models

The exact identity remains a valuable complete-vector regression test, but
the rank gap should not be eliminated as a redundant `r<=6` or `r<=16`
coordinate.  Current data reject the specific low-order `N^(1/4)` mode-front
closure in the intrinsic basis.  A future low-rank transfer model may use the
rank gap as a separate global observable, or change to a basis designed for
the uniform-`p` area functional; it should not infer global neutral width by
truncating the local Krawtchouk jet.

This is one correlated coordinate analysis of existing histograms, not a new
evidence row and not a failure of the exact neutral-area identity.
