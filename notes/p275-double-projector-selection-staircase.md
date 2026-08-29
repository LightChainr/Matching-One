# Issue #275: exact rank-plane selection and the double-projector staircase

Status: one exact zero, one conditional continuum-parity theorem, and one
previously proved Potts-projector zero.  No lattice amplitude is fitted here.

## 1. The normalized rank plane has only one even and one odd line

Use the ordered rank-sector vector

\[
 P=(P_0,P_1,P_2)^T,
 \qquad
 C=\begin{pmatrix}0&0&1\\0&1&0\\1&0&0\end{pmatrix},
\]

where the Alexander/complement involution exchanges ranks zero and two.
The two observable covectors are

\[
 A_{\rm top}=(-1,0,1),\qquad E_{\rm top}=(1,0,1).
\]

Probability tangents obey `(1,1,1) delta P=0`.  Solving this equation together
with `C delta P=+/-delta P` gives the complete normalized tangent space:

\[
 \begin{array}{c|c|c}
 \text{parity}&\delta P& A_{\rm top}\delta P\\ \hline
 C\text{-even}&a(1,-2,1)&0\\
 C\text{-odd}&b(-1,0,1)&2b.
 \end{array}
\]

This is not a model fit.  It is the entire two-dimensional rank-plane
representation.  Once a perturbation's parity is known, its rank-sector
direction is fixed up to one amplitude.

## 2. Vacuum/KdV has an exact zero in `A_top`

Arguin's restricted-sector identity is

\[
 Z_{\rm cross}(Q,\tau)=Q Z_{\rm trivial}(Q,\tau).
\]

At `Q=1`, for every modulus,

\[
 P_2(\tau)=P_0(\tau),\qquad P_0+P_1+P_2=1.
\]

Let `W_tau` be any componentwise linear modulus/Ward operation which commutes
with the rank exchange and annihilates the normalized constant.  The existing
vacuum operator

\[
 K_4=D_2D_0=(\delta-E_2/6)\delta
\]

has exactly these properties.  Applying it to the two identities gives

\[
 (\delta P_0,\delta P_1,\delta P_2)_{K_4}
   =a(1,-2,1),
 \qquad
 \boxed{\delta A_{\rm top}^{K_4}=0}.
\]

Thus the vacuum/KdV spin-4 response is Alexander-even and is killed by the
matching projector at arbitrary `tau`.  This is stronger than reading parity
from a finite-size exponent.  It also clarifies the primitive-sector result in
`p231`: primitive rank-one responses may be nontrivial, but their sum fixes
`K4 P0=K4 P2=-(K4 P1)/2` and cannot generate `A_top`.

## 3. The correct thermal grading makes `Q4 epsilon` odd

The square-site graph alone is not self-matching, so the statement must be
made in the doubled primal/matching continuum family.  Choose the transverse
thermal coordinate `eta` so primal/matching exchange sends `eta -> -eta`:

\[
 P(\eta,\tau)=C P(-\eta,\tau).
\]

Differentiation at the self-dual continuum point gives

\[
 \partial_\eta P|_0=-C\partial_\eta P|_0.
\]

The thermal tangent is therefore exactly the odd rank line,

\[
 (\delta P_0,\delta P_1,\delta P_2)_{\epsilon}
  =b(-1,0,1).
\]

A componentwise Virasoro/Ward descendant does not act on the finite rank
index and hence commutes with `C`.  Therefore

\[
 \boxed{C(Q_4\epsilon)=-(Q_4\epsilon)}
\]

in this doubled continuum grading.  Equivalently, the rank-one response is
zero and the rank-zero/two responses are opposite.

This is a conditional theorem only in the following precise sense: the
finite square-site/matching pair must couple to the declared transverse
continuum thermal coordinate.  It is not the invalid assertion that the
square-site graph possesses an internal self-duality.

### The Jordan partner has no parity loophole

Let `q=Q4 epsilon` and `q_tilde=Q4 epsilon_tilde` satisfy

\[
 (D-21/4)q=0,\qquad (D-21/4)\widetilde q=q.
\]

Suppose the continuum exchange `J` is involutive, commutes with dilations, and
`Jq=-q`.  Write the most general action

\[
 J\widetilde q=bq+c\widetilde q.
\]

The commutator `[J,D]=0` forces `c=-1`; then `J^2=1` forces `b=0`.
Consequently the whole inherited thermal Jordan pair is odd:

\[
 Jq=-q,\qquad J\widetilde q=-\widetilde q.
\]

Adding a bottom state to the top state cannot evade the parity selection.

## 4. Ordered elimination

The Potts half of the selection rule was proved on the Issue #257/#262
branches (`9320649`, `d006f9c`): a regular unlabelled invariant one-insertion
covector annihilates the non-singlet `[2]` four-leg sector.  Charged/twisted
insertions, singular `Q -> 1` normalizations and `Q`-derivative residues remain
separate typed loopholes; none is the ordinary global one-insertion channel.

The resulting staircase is

| candidate | dimension | finite-size power | independent gate | outcome |
|---|---:|---:|---|---|
| vacuum/KdV spin 4 | `4` | `L^-2` | Alexander even | exact `A_top=0` |
| `V_(2,+/-2)` four-leg | `17/4` | `L^-9/4` | Potts `[2]` | ordinary unlabelled overlap zero |
| thermal `Q4 epsilon`/Jordan | `21/4` | `L^-13/4` | singlet and Alexander odd | first listed candidate allowed by both |

This turns the observed `N^-13/8=L^-13/4` law from an isolated exponent match
into an ordered selection statement.  The final field identification still
requires the already-derived `Q4` modulus/Jordan fingerprint.

## 5. Exact boundary: the one missing matrix element

Parity proves that the thermal `Q4` overlap is allowed, not that it is nonzero.
The remaining object is the lattice-to-continuum coupling

\[
 g_{A_{\rm top},Q_4\epsilon}
 =\langle A_{\rm top}\mid Q_4\epsilon\rangle_{\rm lattice\to CFT},
\]

and, for a logarithmic readout, its top-partner analogue.  A nonzero value or
an amplitude-free modulus fingerprint must come from the lattice/CFT bridge;
it cannot be inferred from symmetry alone.

The executable consequence is nevertheless zero-parameter: any proposed
vacuum response must lie on `(1,-2,1)`, while any proposed transverse thermal
response must lie on `(-1,0,1)`, at each modulus before amplitude fitting.

Reproduce with:

```bash
python scripts/p275_double_projector_selection.py \
  --output results/exact-double-projector-selection/latest.json
python -m unittest discover -s tests -p 'test_p275_double_projector_selection.py'
```
