# Matching involution as an RG parity projector

Status: structural scaling argument separated from the more specific `x=21/4` operator identification.

## 1. The exact lattice involution

For a site-matching pair `G, Ghat`, combine model exchange with complementation:

\[
\mathcal T:(G,p)\longleftrightarrow(\hat G,1-p).
\]

At the infinite-lattice critical pair,

\[
p_c(\hat G)=1-p_c(G).
\]

For a finite topological observable `R`, define

\[
S=\frac12\left[R_G(p)+R_{\hat G}(1-p)\right],
\]

\[
D=\frac12\left[R_G(p)-R_{\hat G}(1-p)\right].
\]

Thus `S` is even and `D` is odd under `T`. The matching function is `M=2D` for the wrapping representations used in this repository.

This parity decomposition is more fundamental than any particular correction exponent.

## 2. RG scaling fields can be chosen with matching parity

Near the common continuum fixed point, linearize the RG and the matching involution simultaneously. Because `T^2=1`, irrelevant scaling directions can be decomposed into eigenfields

\[
\mathcal T u_a^{(\eta)}=\eta u_a^{(\eta)},
\qquad \eta=\pm1.
\]

Let `y_a=2-x_a<0` be the RG exponent of one such field. At the critical pair (`z=0`), the contribution of that field to the two microscopic realizations has equal or opposite amplitude according to `eta`.

Therefore, at the critical center,

\[
S(0,L)\supset\frac{1+\eta}{2}a_aL^{y_a},
\]

\[
D(0,L)\supset\frac{1-\eta}{2}a_aL^{y_a}.
\]

Hence the exact selection rule:

\[
\boxed{S(0,L)\text{ sees matching-even fields; }D(0,L)\text{ sees matching-odd fields.}}
\]

This does not require an Ising analogy or an assumed CFT family.

## 3. Add the square-lattice rotation projector

A scalar observable on a microscopic `C4` lattice can receive angular corrections from continuum spins `4j`. For the leading spin-4 component,

\[
\delta R_a\propto L^{y_a}\cos(4\theta)
\]

(up to the sine component when reflection symmetry allows it; the current paired geometries use reflection-related conventions in which the cosine design is the relevant real harmonic).

The same-N projector

\[
P_4[X]=\frac{X(\theta_1)-X(\theta_2)}
{\cos4\theta_1-\cos4\theta_2}
\]

therefore isolates the spin-4 part without changing physical scale or torus modulus.

The combined matching/rotation projectors classify four sectors before any exponent is fitted.

## 4. Away from the center: thermal parity

Under complementation the relevant thermal coordinate reverses direction. After a common metric normalization one may take

\[
z\mapsto-z.
\]

For one scaling field with matching parity `eta` and scaling function `F(z)`, the two lattice realizations contribute schematically

\[
R_G\supset u_aL^{y_a}F(z),
\]

\[
R_{\hat G}\supset\eta u_aL^{y_a}F(-z).
\]

Thus

\[
S\supset\frac12[F(z)+\eta F(-z)]u_aL^{y_a},
\]

\[
D\supset\frac12[F(z)-\eta F(-z)]u_aL^{y_a}.
\]

This gives the parity table:

| matching parity | in `S` | in `D` |
|---|---|---|
| `eta=+1` | thermal-even part | thermal-odd part |
| `eta=-1` | thermal-odd part | thermal-even part |

At `z=0`, this reduces to the exact center selection rule above.

In practice the repository should define thermal reflection through the intrinsic matching curve (`Mbar(p_-)=-u`, `Mbar(p_+)=+u`) rather than rely on disputed last digits of `p_c` or on equal microscopic metric factors.

## 5. Root-bias exponent identifies the dimension of the first odd field

Suppose the first matching-odd spin-4 field surviving in `D(0,L)` has scaling dimension `x_o`. Then

\[
M_L(p_c)\sim L^{2-x_o}.
\]

The thermal derivative obeys

\[
M'_L(p_c)\sim L^{y_t},
\qquad y_t=1/\nu=3/4.
\]

Linearizing the root gives

\[
p_L^*-p_c\sim L^{2-x_o-y_t}.
\]

If the observed root bias is `L^-w`, then

\[
\boxed{x_o=2+w-y_t=x_t+w},
\]

because `x_t=2-y_t`.

For percolation, `x_t=5/4`; for the observed `w=4`,

\[
\boxed{x_o=5/4+4=21/4.}
\]

This is why a level-4 thermal-family descendant is structurally natural: the correction dimension inferred from the root is **exactly the thermal primary dimension plus the square-lattice spin level**.

This relation is independent of the detailed Virasoro construction.

## 6. Current empirical assignment

The 100M P31 data support two separate H4 laws:

\[
P_4[S]\sim A_I N^{-1}=A_IL^{-2},
\]

and

\[
P_4[M]\sim A_M N^{-13/8}=A_ML^{-13/4}.
\]

The first has the dimension expected of a matching-even `x=4` spin-4 field. The second has the dimension inferred above for a matching-odd `x=21/4` spin-4 field.

The specific CFT-family labels remain hypotheses:

- matching-even `x=4`: identity-family square anisotropy (`T^2+Tbar^2` / its `c=0` logarithmic completion);
- matching-odd `x=21/4`: level-4 thermal-family spin-4 quasiprimary / logarithmic completion.

## 7. Independent predictions

The parity-projector mechanism predicts more than the existing fits:

1. A self-dual `C4` critical model should suppress the central duality-odd sector while permitting the even `x=4` anisotropy (issue #42).
2. Gaussian `pi/4` doubling must reverse the spin-4 sign; if the odd field has `x=21/4`, the fixed magnitude ratio is `2^-13/8` (issue #38).
3. The root amplitude must be normalized by `Delta cos4theta`:

\[
-\frac{N^2\Delta p^*}{\Delta\cos4\theta}\to A_M/B.
\]

4. The thermal-even part of `P4[D]` contains the central root-moving field; P33 already supports this qualitative parity assignment.
5. A new held-out same-N pair should satisfy **both** the `N^-1` even-sector law and the `N^-13/8` odd-sector law with amplitudes frozen before simulation (issue #43).

## 8. What would falsify the structural mechanism

The broad matching-projector mechanism is weakened if

- the center `D` receives a robust leading matching-even contribution under an independently verified parity convention;
- the `S/D` center selection does not survive exact small-system and full-curve checks;
- a self-dual control exhibits the same nonzero odd central amplitude without an additional odd microscopic coupling;
- held-out angular data require spins not compatible with `C4` harmonics.

The **specific** `x=21/4` identification is more easily falsified: prospective fixed-ratio and held-out exponent tests can fail even while the general matching-parity projector remains correct.
