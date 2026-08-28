# Matching-involution thermal-parity derivative spectrum

Status: falsifiable selection rule for the current two-spin-4 model.

## 1. Linearized involution statement

Near the common scaling fixed point, let an irrelevant eigenfield have scaling dimension `x`, rotation spin `s`, and matching-involution parity

\[
\eta=\pm1.
\]

For a fixed torus modulus, write its contribution to one microscopic member of the matching pair as

\[
\delta R_G=uL^{2-x}H_s(\theta)F(z),
\qquad z\propto(p-p_c)L^{y_t},\quad y_t=3/4.
\]

In a parity eigenbasis the matching partner contributes

\[
\delta R_{\hat G}=\eta uL^{2-x}H_s(\theta)F(-z).
\]

Define

\[
S=(R_G+R_{\hat G})/2,\qquad D=(R_G-R_{\hat G})/2.
\]

Then

\[
\delta S=\frac{uL^{2-x}H_s}{2}[F(z)+\eta F(-z)],
\]

\[
\delta D=\frac{uL^{2-x}H_s}{2}[F(z)-\eta F(-z)].
\]

Expanding `F(z)=sum f_n z^n` yields the derivative selection rule

\[
S^{(n)}(p_c)\neq0\quad\Longleftrightarrow\quad(-1)^n=\eta,
\]

\[
D^{(n)}(p_c)\neq0\quad\Longleftrightarrow\quad(-1)^n=-\eta,
\]

for that individual parity eigenfield, modulo accidental zero amplitudes.

Every `p` derivative also contributes `L^(n y_t)`, so an allowed derivative scales as

\[
L^{2-x+n y_t}.
\]

This statement is stronger than the center-only parity rule because it predicts an alternating tower of observables.

## 2. Apply the spin-4 projector first

The raw `D=M/2` contains the orientation-independent universal crossover function, e.g. `D'(pc)~L^(3/4)`.  To isolate lattice irrelevant fields, first project two same-N orientations:

\[
P_4[X]=\frac{X(\theta_1)-X(\theta_2)}{\cos4\theta_1-\cos4\theta_2}.
\]

The leading orientation-independent scaling function cancels exactly in this same-shape difference.

## 3. Frozen predictions for the two-field model

### A. matching-even identity-family spin 4

Working assignment:

\[
x_I=4,\qquad \eta_I=+1.
\]

Allowed projected derivatives at the center are

\[
P_4[S]\sim L^{-2}=N^{-1},
\]

\[
P_4[D']\sim L^{-5/4}=N^{-5/8},
\]

\[
P_4[S'']\sim L^{-1/2}=N^{-1/4},
\]

with alternating even derivatives in `S` and odd derivatives in `D`.

The first line is already supported by the 100M P31 sector data.  The second line is a new independent parity test.

### B. matching-odd thermal-family spin 4

Working assignment:

\[
x_T=21/4,\qquad \eta_T=-1.
\]

Allowed center derivatives are

\[
P_4[D]\sim L^{-13/4}=N^{-13/8},
\]

\[
P_4[S']\sim L^{-5/2}=N^{-5/4},
\]

\[
P_4[D'']\sim L^{-7/4}=N^{-7/8},
\]

with even derivatives in `D` and odd derivatives in `S`.

The first line has passed the prospective Gaussian-doubling test.  `P_4[S']~N^-5/4` is a new independent matching-parity prediction.

## 4. Why this is useful

The center amplitudes alone could in principle be mimicked by unrelated powers.  The derivative spectrum forces the **same two field assignments** to predict four additional exponent/parity combinations.

In particular:

- `P4[D']` should be controlled by matching-even fields, not by the `x=21/4` center field;
- `P4[S']` should be controlled by matching-odd fields, not by the `x=4` center field.

This cross-channel alternation is difficult to fake with a single effective power law.

## 5. Logarithmic caveat

At `c=0` the thermal Kac operator belongs to a logarithmic multiplet.  A logarithmic partner can multiply the predicted powers by polynomials in `log L`; it does not change the matching/thermal parity alternation itself.  Therefore model comparisons should hold the power fixed first and compare no-log versus log-amplitude variants on held-out sizes.

## 6. Numerical protocol

Extend the threshold-rank analyzer so the separate `K_plus` and `K_minus` distributions reconstruct the primal and matching probabilities and derivatives individually, not only their difference.

For each frozen same-N pair compute at an intrinsic center (preferably the `Mbar=0` center from the full curves):

- `P4[S]`, `P4[D]`;
- `P4[S']`, `P4[D']`;
- if stable, `P4[S'']`, `P4[D'']`.

Fit the four primary powers only on training sizes and score held-out sizes:

```text
P4[S]   : N^-1
P4[D]   : N^-13/8
P4[D']  : N^-5/8
P4[S']  : N^-5/4
```

Do not fit four independent free exponents first.  Free-exponent/log alternatives are secondary falsification models after the frozen power test.
