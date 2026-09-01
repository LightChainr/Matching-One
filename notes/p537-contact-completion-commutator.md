# Contact fusion as a noncommuting completion operator

Status: synthesis of the exact N25 contact-stage tensor and the frozen held-out
N65 full-root result.  The finite statements below are results; the continuum
operator and limiting matrix are explicitly labelled as mechanism hypotheses.

## 1. The finite tensor measures a commutator

Let the row space be the two essential-birth transitions

\[
V_b=\operatorname{span}\{|01\rangle,|12\rangle\}
\]

and let the column space contain exchange-even single contact `s` and fused
double contact `d`.  Write the complete pooled-root Schur allocation as

\[
K_N=\begin{pmatrix}
K_{01,s}&K_{01,d}\\
K_{12,s}&K_{12,d}
\end{pmatrix},\qquad
\rho_r={K_{r,d}\over K_{r,s}}.
\]

There is an exact algebraic identity

\[
\det K_N=K_{01,s}K_{12,s}(\rho_{12}-\rho_{01}).             \tag{1}
\]

At N25 and N65 the ratios are

| N | `rho_01` | `rho_12` | `rho_12/rho_01` |
|---:|---:|---:|---:|
| 25 | `+1.0185467` | `-1.1638779` | `-1.1426848` |
| 65 | `+0.5845570` | `-1.1958736` | `-2.0457776` |

Thus the negative determinant is a stage-dependent fusion eigenvalue, not a
second small marginal effect.  If single contact is represented by the birth
raising map `B` and double contact by `F B`, then

\[
K_{r,s}=b_r,\qquad K_{r,d}=f_{r+1}b_r,
\]

and

\[
\det K_N=b_0b_1(f_2-f_1),\qquad
[F,B]|r\rangle=(f_{r+1}-f_r)b_r|r+1\rangle.                \tag{2}
\]

Equations (1)--(2) give the minimal operator content of the finite result:
contact fusion does not commute with passage through topological completion.
A scalar connected-contact coefficient, `K_(r,d)=lambda K_(r,s)`, would give
`det K=0` and cross-ratio one.  It is excluded by both the exact N25 tensor
and the independent N65 production.

## 2. Two-scale fingerprint

The N25-to-N65 effective powers, written in the site count `N`, split sharply:

| cell | effective `q` in `K_N proportional to N^-q` |
|---|---:|
| entry x single | `3.04044` |
| completion x single | `2.97846` |
| completion x double | `2.95007` |
| entry x double | `3.62157` |

Three cells therefore share an `N^-3` backbone, while entry-double is already
numerically at `N^-29/8` (`29/8=3.625`).  Equivalently, after the common
`N^3` rescaling, the N65/N25 cell ratios are approximately

\[
\begin{pmatrix}
0.962&0.552\\
1.021&1.049
\end{pmatrix}.                                               \tag{3}
\]

The exceptional exponent is not only a rare-carrier effect.  Use the positive
source-normalized exposure `E` and conditional signed density `D=K/E`.  For
entry-double,

\[
\begin{array}{c|cc|c}
&N=25&N=65&\text{effective power}\\ \hline
E&3.2694090\,10^{-2}&2.7268969\,10^{-4}&N^{-5.00948}\simeq N^{-5}\\
D&-8.9841554\,10^{-5}&-3.3839408\,10^{-4}&N^{+1.38791}\simeq N^{11/8}\\
K=ED&-2.9372879\,10^{-6}&-9.2276578\,10^{-8}&3.62157\simeq29/8
\end{array}                                                   \tag{4}
\]

The other three cells have exposure powers between `4.65` and `5.19` and
conditional-density growth powers between `1.67` and `2.15`, giving their
common net `N^-3` behavior.  Entry-double retains the same `N^-5` occurrence
law but receives only `N^(11/8)` rather than `N^2` conditional amplification.
Thus the exceptional cell factorizes numerically as

\[
N^{-29/8}=N^{-5}\,N^{11/8},                                 \tag{5}
\]

and its extra `N^-5/8=L^-5/4` suppression relative to the backbone lies in
the conditional Schur/OPE strength, not in event frequency.  The exponent
`5/4` is the percolation thermal scaling dimension, making a thermal insertion
the minimal candidate for the missing conditional amplification.  Calling
the effect a rarer contact event is therefore directly contradicted by the
frozen exposure decomposition.

## 3. Continuum mechanism hypothesis

The exact support rule says that a kernel-changing diagonal edge has no
contact-free radius-one representative; one-arm contact occurs only in the NN
source orbit, while every non-NN source orbit requires double contact.  This
identifies a connected two-arm collision as the local vertex.  The sign
rotation in (1), however, requires an additional completion action.  The
minimal proposed operator is

\[
\boxed{
\mathcal O_{537}=
\Pi_{H_4}^{\rm odd}\,
\mathcal C_{2\mathrm{arm}}^{\rm conn}\,
\Pi_{\mathrm{rank}=2}\,B .}
                                                                  \tag{6}
\]

Here `Pi_H4^odd` supplies the spin-four/matching-odd sector,
`C_2arm^conn` is the contact OPE vertex, and the rank-two projector is the
topological completion gate.  This is primarily a noncommuting rank-birth
transfer.  A Q4 Jordan pair may dress its scale evolution, but Jordan mixing
is not needed to produce either the finite support rule or the fixed-scale
negative determinant.  A Jordan identification would require its independent
affine-log and modulus fingerprints; rank two alone is not that evidence.

The bridge to the reduced `35/8` response uses the same operator matrix
element but a different normalization.  Let
`P_N=M'_N(p_N)` be the total pivotal/birth metric.  Under the two explicit
assumptions that the root-Schur projection removes the scalar thermal tangent
and that the surviving completion state is L2-normalized in this metric,

\[
G_4(N)=N^{-2}P_N^{-1/2}\,
\langle H_4|[F,B]|\mathrm{source}\rangle+o(N^{-2}P_N^{-1/2}). \tag{7}
\]

Since `P_N proportional to N^(3/8)`, equation (7) gives

\[
G_4\sim N^{-35/16}=L^{-35/8},\qquad
\Xi_N=N^2\sqrt{P_N}\,G_4\longrightarrow\mathfrak c_{FB}.    \tag{8}
\]

The observed `Xi` near eight is therefore naturally read as the dimensionless
completion-commutator matrix element.  Equation (8) does not identify the raw
contact tensor with a Q4 one-point function: H4 fixes the spin/parity sector,
whereas `35/8` is the exponent of the pair- and pivotal-normalized composite
coordinate.

## 4. Falsifiable triangular limit

Equations (3)--(5) motivate the sharper asymptotic hypothesis

\[
N^3K_N\longrightarrow
\begin{pmatrix}
-a&0\\
-b&c
\end{pmatrix},\qquad a,b,c>0.                               \tag{9}
\]

It makes parameter-free structural predictions:

\[
\rho_{01}=O(N^{-5/8})\to0^+,
\qquad \rho_{12}\to-c/b<0,
\qquad \det K_N\sim-ac\,N^{-6}<0.                          \tag{10}
\]

Consequently the projective cross-ratio diverges through the negative real
axis rather than approaching the scalar-OPE value one.  Failure of the zero,
sign, or `5/8` separation in a genuinely held-out larger size would reject
the completion-projector form (9) without inventing another contact label.

This note does not recommend a new descriptor, a contact-radius scan, or a
repeat of the carrier production.  The next useful work is an analytic
derivation of (6) from the Alexander/birth filtration, or scoring (9)--(10)
only when an independently motivated future block already contains the same
four frozen cells.
