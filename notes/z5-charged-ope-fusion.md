# Z5 charged OPE fusion after the chiral GLS freeze

This is the smallest charged three-point continuation of the norm-five
handed-response experiment.  It uses the exact deck group, not a fitted
latent charge.

## Exact algebra

Write the characters of the norm-five deck group as

\[
\chi_r(k)=\zeta_5^{rk},\qquad r,k\in\mathbb Z/5\mathbb Z.
\]

Their fusion and invariant tensors are

\[
\chi_r\otimes\chi_s=\chi_{r+s},\qquad
g_{rs}\ne0\Longrightarrow r+s=0,
\qquad C_{rst}\ne0\Longrightarrow r+s+t=0.
\]

Among triples with three nontrivial charges and no neutral proper subset,
there are only four unordered channels:

\[
(113),\ (122),\ (244)=\overline{(113)},\
(334)=\overline{(122)}.
\]

Thus a fixed deck generator leaves exactly two primitive chiral cubic
channels, `A=113` and `B=122`.  Every one-point factor and every disconnected
two-point factor in these channels vanishes by exact charge conservation.
Their raw third moment is already connected; no noisy disconnected
subtraction is needed.

## What can and cannot be normalization-free

For charged fields, the phase of a single three-point coefficient is not a
universal number until a charged basis has been fixed.  Under independent
rescalings \(O_r\mapsto\lambda_r O_r\), the exponent matrix of

\[
(C_{113},C_{122},C_{244},C_{334})
\]

with respect to \((\lambda_1,\lambda_2,\lambda_3,\lambda_4)\) is

\[
E=\begin{pmatrix}
2&1&0&0\\
0&2&1&0\\
1&0&0&2\\
0&0&2&1
\end{pmatrix},\qquad \det E=-15.
\]

The nonzero determinant is a useful no-go: no nonconstant monomial made only
from these four cubic coefficients removes all charged-field rescalings.
Two-point normalization removes magnitudes but leaves charged basis phases.

The phase-gauge-free charged analogue of a squared normalized structure
constant is instead

\[
\mathcal I_{rst}=
\frac{C_{rst}C_{-r,-s,-t}C_{000}}
{C_{0,r,-r}C_{0,s,-s}C_{0,t,-t}}.
\]

Every \(\lambda_a\) cancels exactly.  Under reflection positivity this is the
normalized modulus squared.  A genuinely complex OPE phase requires either
the repository's exact transported deck basis or a larger closed fusion loop,
such as a four-point crossing invariant.  This is why the current chiral GLS
basis transport is scientific content rather than bookkeeping.

## Frozen next prediction

Use the same N325 parent and the same two hands as the existing chiral GLS.
At one randomly translated fixed-shape three-anchor pattern, measure the two
primitive channels in both hands:

\[
(C_{A,+},C_{A,-},C_{B,+},C_{B,-}).
\]

The first score is the zero-parameter complex closure

\[
\boxed{C_{A,+}C_{B,-}-C_{A,-}C_{B,+}=0}. 
\]

It has two real degrees of freedom and avoids division by a noisy channel.
It says the two allowed fusion channels share one handed eigenphase.  Failure
would show that charge fusion resolves more than the current one-field chiral
GLS, even if one channel separately happens to prefer H4.

Conditioned on a single spin-\(s\) Hecke eigenfield dominating the same local
angular insertion on all three legs,

\[
C_{A,+}=q_{3s}a,\quad C_{A,-}=a,\qquad
C_{B,+}=q_{3s}b,\quad C_{B,-}=b,
\]

where

\[
q_{3s}=\left(\frac{2+i}{2-i}\right)^{3s}.
\]

A joint eight-real-component GLS fits only the two complex amplitudes
\(a,b\), leaving four degrees of freedom.  The H4, H8 and H12 hypotheses make
the distinct frozen predictions \(q_{12},q_{24},q_{36}\).  This is the direct
cubic closure of the already-frozen linear handed phase; it does not refit an
exponent or allow a separate phase per charge channel.

Also report \(\mathcal I_{113}\) and \(\mathcal I_{122}\) with their joint
covariance.  These are the clean normalization-free OPE magnitudes.  Do not
call a raw complex cubic phase universal without the declared transported
deck convention.

## Minimal measurement

Reuse the N325 common Bernoulli field and charged landing-marked local row.
For each replica, translate one fixed three-anchor pattern uniformly, evaluate
the three local rows, and accumulate `113` and `122` for both hands plus their
reflected conjugates.  Save the full 8x8 covariance of the four complex
primary channels.  Nonneutral triples are exact labeling-error nulls.

This costs roughly three marked insertions per replica and asks a new
operator-algebra question: whether the charged sectors share one fusion
eigenphase.  It is not another estimate of the existing global H4 mean.
