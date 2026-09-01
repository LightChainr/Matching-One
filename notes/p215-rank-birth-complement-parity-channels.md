# Issue #215: complement parity of the canonical rank-birth insertion

Status: exact complement/gate algebra, exhaustive coefficient oracles, and a
conditional mechanism inference. This is an increment over `3881e88`, not a
repeat of its Russo proof.

## 1. Complement reverses the insertion endpoints

Fix all sites other than `v`. For the black primal insertion write

\[
 r_B^0=r_B(v=0),\qquad r_B^1=r_B(v=1).
\]

The matching insertion must itself run from white `v=0` to white `v=1`.
Because white occupancy is the black complement, these are the complements of
the black `v=1` and `v=0` states, respectively. Digital Alexander rank
complementarity therefore gives

\[
 r_W^0=2-r_B^1,\qquad r_W^1=2-r_B^0.
\]

Using the gate definitions from `3881e88`, this immediately yields

\[
 \boxed{I_{01}^{B}=I_{12}^{W}},\qquad
 \boxed{I_{12}^{B}=I_{01}^{W}}.
\]

A direct `0->2` black insertion maps to a direct `0->2` matching insertion and
has both gates equal to one.

## 2. The canonical even and odd insertion channels

Diagonalize the two-gate exchange:

\[
 S=I_{01}+I_{12},\qquad D=I_{12}-I_{01}.
\]

Under the combined primal/matching, `p <-> 1-p`, insertion-reversal map,

\[
 S_G(p)=S_{\widehat G}(1-p),\qquad
 D_G(p)=-D_{\widehat G}(1-p).
\]

Their integrated meanings are exact. Since

\[
 r=\mathbf1_{r\ge1}+\mathbf1_{r=2},
\]

the even density is

\[
 \boxed{\sum_v\mathbb E[S_v]=f_{01}+f_{12}=M'(p)}.
\]

Meanwhile

\[
 P_1=P(r=1)=P(r\ge1)-P(r=2),
\]

so the odd density is

\[
 \boxed{\sum_v\mathbb E[D_v]=f_{12}-f_{01}
       =-\partial_pP_1(p)=\partial_p(P_0+P_2)}.
\]

Thus the local parity split exactly mirrors the canonical global rank plane:

```text
S: derivative of A_top=P2-P0=M, exchange-even after the thermal derivative;
D: derivative of E_top=P0+P2=1-P1, exchange-odd after the thermal derivative.
```

This also locates the persistence-clock interpretation. `S` is the local
translation/total-birth susceptibility; `D` is the local derivative of the
rank-one lifetime mass.

### Direct `0->2` is automatically sorted

A simultaneous birth contributes

```text
S=2,
D=0.
```

It belongs entirely to the even susceptibility. The odd channel has no
simultaneous-birth ambiguity.

## 3. Primitive-line refinement and a computable rotation character

For a nonsimultaneous birth, complement exchange maps `0->1` to `1->2` while
preserving the canonical rank-one plateau line `ell`. The tiny oracle also
finds identical black/white endpoint indices `iota=1`; index equality is
reported only as an exhaustive small-size result, not promoted to a theorem.

Let `(x,y)` be the physical lifted period vector of `ell`. Because `ell` is
projective, the lowest square-lattice angular character that is insensitive to
its sign is

\[
 \chi_4(\ell)=\frac{(x+iy)^4}{(x^2+y^2)^2}
 =\cos4\theta+i\sin4\theta.
\]

It obeys exactly

\[
 \chi_4(-\ell)=\chi_4(\ell),\qquad
 \chi_4(R_{\pi/2}\ell)=\chi_4(\ell),\qquad
 \chi_4((1+i)\ell)=-\chi_4(\ell).
\]

Since complement preserves `ell`, line refinement does not change complement
parity:

```text
chi4(ell) S  is matching-even and spin-4,
chi4(ell) D  is matching-odd  and spin-4.
```

Every nonzero `D` contribution has a non-null `ell`: direct `0->2` insertions
cancel from `D`, and strict `0->1`/`1->2` insertions carry the plateau line.
Consequently the odd channel admits a complete canonical projective-H4
decomposition. The even channel instead contains an additional unpolarized
simultaneous-birth sector with no canonical line.

This is a real acquisition advantage, not only a representation label.

## 4. Local landing H4 has the same parity split

The existing landing mark is symmetric under exchanging its open and closed
component lists while simultaneously swapping NN/matching adjacency. It is
therefore preserved by the same complement-reversal map. Multiplying it by
`S` or `D` produces local matching-even and matching-odd H4 fields without
changing the unmarked identities.

At `p=1/2`, the exact rooted results are

| geometry | `S` | `D` | `Re chi4*S` | `Re chi4*D` | `landingH4*S` | `landingH4*D` |
|---|---:|---:|---:|---:|---:|---:|
| axis `L=2` | `3` | `0` | `2` | `0` | `0` | `0` |
| Gaussian `(2,1)` | `25/8` | `-5/8` | `-7/8` | `7/40` | `0` | `0` |
| axis `L=4` | `4209/1024` | `-17/16` | `795/256` | `-109/128` | `467/1024` | `-49/1024` |

The Gaussian control also has nonzero imaginary/projective components:

```text
Im chi4*S = 3,
Im chi4*D = -3/5.
```

Every primal/matching gate, line and local-mark relation is checked
environment by environment. The resulting even/odd, line-H4 and landing-H4
polynomials are then independently compared after the exact substitution
`p -> 1-p`; every coefficient residual is zero.

## 5. Mechanism inference

The exact statement is now sharper than “plus/minus pivotal counts”:

- `S_H4` is the canonical matching-even local derivative field whose unmarked
  integral is exactly `M'`.
- `D_H4` is a canonical matching-odd local field whose unmarked integral is
  the derivative of rank-one persistence, not the matching susceptibility.
- `chi4(ell)D` is fully defined on every nonzero odd insertion and has the
  correct sign reversal under the Gaussian `1+i` orientation change.

Conditional mechanism hypothesis: if the anomalous global matching-odd H4
response is carried by the essential-H1 lifetime mode, the line-resolved or
landing-resolved `D` channel should couple to it more directly than the old
untyped pivotal mark. In clock language this predicts a closer association
with lifetime/width `W` than with center/translation `C`.

This hypothesis is deliberately strong and cheap to falsify. A common-field
production stream need only retain, at each rank birth,

```text
gate sign (+1 for second, -1 for first),
ell and optionally iota,
local landing H4.
```

It then obtains `S_H4`, `D_H4` and `chi4(ell)D` from one correlated raw block.

## Boundary

- Gate exchange and `S/D` parity follow exactly from digital rank
  complementarity.
- Coefficient identities are exhaustive on the declared controls.
- Equality of integral indices is only a tiny-oracle result.
- `chi4(ell)` is a global winding-direction character, distinct from the local
  landing geometry.
- No nonzero continuum overlap, exponent, or Q4/Jordan identification is
  asserted by the exact algebra alone.

Reproduce with:

```bash
python scripts/rank_birth_parity_channels.py \
  --output results/rank-birth-parity-channels/latest.json
python scripts/rank_birth_parity_channels.py --format markdown \
  --output results/rank-birth-parity-channels/latest.md
python -m unittest discover -s tests -p 'test_rank_birth_parity_channels.py'
```

