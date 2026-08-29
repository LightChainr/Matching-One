# Exact rectangular Ward/Hecke ratio for the thermal Q4 candidate

Status: exact modular-form calculation plus a conditional lattice-observable
prediction.  The mathematical ratio is exact; its application to the matching
observable requires the explicit bridge stated below.

## 1. Exact weight-4 calculation

Let `E4` be the normalized Eisenstein series.  Its normalized prime-2 Hecke
eigenvalue is

\[
T_2E_4=(1+2^3)E_4=9E_4.
\]

At the square CM point this reads

\[
8E_4(2i)+\frac12\left[E_4(i/2)+E_4((1+i)/2)\right]
=9E_4(i).
\]

The weight-4 modular transformation law gives

\[
E_4(i/2)=16E_4(2i),
\qquad
E_4((1+i)/2)=-4E_4(i).
\]

Writing `x=E4(2i)/E4(i)` therefore gives `16x-2=9`, hence

\[
\boxed{\frac{E_4(2i)}{E_4(i)}=\frac{11}{16}}.
\]

This derivation is exact rational arithmetic; it does not use a truncated
`q`-series or any percolation data.

## 2. The area-normalized spin-4 shape

For coordinate periods `ell*(1,tau)`, the area is

\[
A=\ell^2\operatorname{Im}\tau.
\]

The Ward identity for the thermal level-4 quasiprimary contains
`ell^-4 E4(tau)`.  At fixed area this is

\[
A^{-2}\widehat E_4(\tau),
\qquad
\widehat E_4(\tau)=(\operatorname{Im}\tau)^2E_4(\tau).
\]

Consequently the exact physical-frame ratios are

\[
\boxed{
\frac{\widehat E_4(2i)}{\widehat E_4(i)}=\frac{11}{4}
},
\qquad
\boxed{
\frac{\widehat E_4((1+i)/2)}{\widehat E_4(i)}=-1
}.
\]

The second identity is the spin-4 phase under the modular basis change that
maps `(1+i)/2` back to the square torus.  It is a useful sign/convention
control, but it is not an independent physical shape measurement.

## 3. Conditional matching-observable bridge

Let `D4_N(tau)` be a fully typed, same-modulus H4 projection of the central
matching-odd response, and let

\[
K_N(\tau)=\bar M'_N(p_{0,N};\tau)
\]

be the intrinsic center slope.  If `D4` is dominated by the ordinary
`c=0,h=5/8` thermal `Q4` descendant and `K` supplies the same thermal-primary
block and metric, then

\[
\mathcal W_N(\tau)
=N^2\frac{D4_N(\tau)}{K_N(\tau)}
\longrightarrow C\widehat E_4(\tau).
\]

The unknown lattice coupling `C` cancels, giving the conditional no-fit score

\[
\boxed{
\mathcal W(2i)-\frac{11}{4}\mathcal W(i)=0
}.
\]

This is a hypothesis about the lattice-to-CFT observable bridge.  The exact
modular identity does **not** establish that the matching observable obeys it.
Additional fields, a defect/vector-valued response, a logarithmic block, or an
incorrect H4/channel projector can all violate the conditional score without
contradicting the modular-form calculation.

## 4. Minimal experiment

Use the general integer-period backend to construct a rectangular `tau=2i`
pair with two microscopic square-lattice embeddings at the same modulus and
area.  Freeze the topology channel, primal/matching parity, embedding order,
period-basis convention and exact H4 response matrix before reading targets.

Recompute `p0`, `D4`, `K` and `W` inside every synchronized delete-one
replicate.  Score the `11/4` residual with the square-torus source covariance
and rectangular target covariance.  The modularly equivalent diagonal view
may be used as a sign/implementation check, but must not be counted as an
independent scientific result.

The rectangular test is complementary to the Pell--hexagonal zero: its H4
signal is enhanced rather than annihilated, so it is a cheaper first check of
the proposed full finite-modulus Ward fingerprint.

## 5. Reproduction

```bash
python3 scripts/derive_rectangular_thermal_q4_hecke.py
python3 -m unittest tests.test_rectangular_thermal_q4_hecke -v
```

The machine-readable freeze is
`predictions/rectangular_thermal_q4_hecke_20260829.yaml`.
