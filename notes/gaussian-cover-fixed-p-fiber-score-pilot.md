# Fixed-\(p\) Gaussian-cover fiber-score pilot

Issue: #226.  Status: exact convention plus a 20,000-replica smoke test, not
production evidence.

## Exact fiber convention

For a primitive Gaussian parent \(z=a+ib\), the cyclic site label is

\[
j_p=ax+by\pmod N,\qquad N=a^2+b^2.
\]

Multiplication by \(1+i\) gives the raw child
\(z(1+i)=(a-b)+i(a+b)\), of order \(2N\), with label

\[
j_c=(a-b)x+(a+b)y\pmod {2N}.
\]

The multiplier \(t\pmod N\) solving
\(t(a-b)=a\) and \(t(a+b)=b\) gives the exact child-to-parent map

\[
j_p=t(j_c\bmod N)\pmod N.
\]

Consequently the two sites in every fiber are exactly \(j\) and \(j+N\).
For the two frozen norm-two lineages this gives:

| parent | raw child | \(t\pmod {65}\) |
|---|---|---:|
| \((8,1)\) | \((7,9)\) | 29 |
| \((7,4)\) | \((3,11)\) | 24 |

The raw children are the existing canonical \((9,7)\) and \((11,3)\)
geometries up to the already-known dihedral convention.  Keeping the raw
representatives makes the cover map and the common \(+45^\circ\) rotation
literal rather than inferred from vertex-array order.

## Scores and observables

For one Bernoulli child field \(X_0,\ldots,X_{129}\), the two exact
likelihood-score numerators are

\[
S_+=\sum_{j=0}^{129}(X_j-p),\qquad
S_-=\sum_{j=0}^{64}(X_j-X_{j+65}).
\]

Thus \(S_+\) is deck-trivial and \(S_-\) is the nontrivial \(\mathbb Z/2\)
character.  For any child observable \(O\),

\[
\partial_\epsilon E_{p+\epsilon h}[O]\big|_{0}
=\frac{E_p[O S_h]}{p(1-p)}.
\]

The primitive observations are matching-odd cross-wrap values
\(M_1,M_2\) on the two raw children.  The archived analysis transforms them
to

\[
M_{\rm global}=(M_1+M_2)/2,\qquad H_4=M_1-M_2,
\]

for both score directions, retaining the full four-by-four batch covariance.
The code deliberately does not construct a parent Bernoulli configuration:
there is no canonical deterministic collapse of the two child bits, and none
is needed for the likelihood response.

## Exact and smoke gates

The tiny \((2,1)\to(1,3)\), \(N=5\to10\) oracle enumerates all \(2^{10}\)
states at \(p=2/5\).  For matching-odd cross wrap and a directional
matching-odd \(H_4\) readout, in both trivial and detail directions, the score
derivative equals an independent degree-ten polynomial derivative reconstructed
from eleven rational finite-difference nodes.  All four equalities are stored
in the JSON.

The frozen smoke command is:

```bash
python3 scripts/gaussian_cover_fiber_score_pilot.py \
  --samples 20000 --batches 20 --workers 4 \
  --p 0.592746050790 --seed 22620260829 \
  --output results/local-20260829/P226-gaussian-cover-fiber-score-smoke/fiber_score_smoke.json
```

At smoke resolution the global/trivial response is
\(11.0051\pm0.1063\).  The global/detail and both orientation responses are
within one standard error of zero.  This is a kernel and covariance gate only;
it is not evidence for absence of a detail channel.  A production comparison
must freeze sample size independently and should add parent/child scale
combinations before interpreting an RG eigenvalue.

