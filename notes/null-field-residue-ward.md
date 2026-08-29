# Null-field residue and a direct thermal-Q4 versus four-leg Ward gate

The level-two thermal null is strong enough to distinguish four mechanisms
without fitting a radial exponent.  The same algebra also gives a direct
positive-mode test separating the thermal Q4 descendant from the new
`V_(2,2)` four-leg primary adversary.

## Parent null residue

At \(c=0,h=5/8\), set

\[
N_2=L_{-2}-\frac23L_{-1}^2.
\]

For a fixed-\(c\) Jordan top normalized by

\[
(L_0-h)|\psi\rangle=|\phi\rangle,
\qquad L_n|\psi\rangle=0\quad(n>0),
\]

the null no longer vanishes homogeneously.  With
\(R=N_2|\psi\rangle\), exact Virasoro commutators give

\[
L_1R=-\frac83L_{-1}|\phi\rangle,
\qquad L_1^2R=-\frac{10}{3}|\phi\rangle,
\qquad L_2R=0.
\]

This is an inhomogeneous Ward identity, not a log fit.

There is a sharply different parameter-derivative mechanism.  The level-two
Kac curve is

\[
c(h)=\frac{2h(5-8h)}{2h+1},\qquad c'(5/8)=-\frac{40}{9}.
\]

Differentiating \(N_2(h)|\phi_h\rangle=0\) along this curve gives

\[
N_2|\partial_h\phi\rangle=-\frac{16}{27}L_{-1}^2|\phi\rangle,
\]

so the first two rows remain \((-8/3,-10/3)\), but

\[
L_2R=-\frac{20}{9}|\phi\rangle.
\]

The ratios are therefore

```text
fixed-c Jordan:       b/a = 5/4, d/b = 0
Kac-curve derivative: b/a = 5/4, d/b = 2/3
ordinary thermal:     a=b=d=0
```

For any ordinary primary of chiral weight \(h_*\), applying the thermal null
operator gives

\[
a=\frac{5-8h_*}{3},\qquad b=2h_*a,\qquad d=0.
\]

Thus a negative-weight primary has the exact sign pattern \(a>0,b<0\).
The spin-\(+4\) four-leg primary has \(h_*=33/8\), hence

\[
(a,b,d)=(-28/3,-77,0),\qquad b/a=33/4.
\]

Its opposite chirality, \(h_*=1/8\), gives \((4/3,1/3,0)\).

## Direct spin-4 operator gate

The repository Q4 convention is

\[
Q_4=40L_{-2}^2-60L_{-3}L_{-1}-9L_{-4}.
\]

Commuting positive modes through this state gives three scalar contractions:

\[
\left(L_2^2,L_1L_3,L_4\right)Q_4|\epsilon\rangle
=(40,-60,30)|\epsilon\rangle.
\]

The primitive ratio is

\[
\boxed{4:-6:3}.
\]

An ordinary spin-4 primary, including `V_(2,-2)` with \(x=17/4\), is killed
by all positive Virasoro modes and gives `(0,0,0)`.  This is the cleanest
exact Q4-versus-four-leg gate available in the current program: it asks
whether the field is a descendant, not which of two nearby powers fits
better.

## Basis-independent Jordan residue after Q4

Let

\[
W=(W_{22},W_{13},W_4)
 =(L_2^2,L_1L_3,L_4)Q_4.
\]

For the fixed-\(c\) Jordan top,

\[
W|\psi\rangle=(40,-60,30)|\psi\rangle
 +(864,-546,48)|\phi\rangle.
\]

The two covectors

\[
J_A=3W_{22}+2W_{13},\qquad J_B=W_{13}+2W_4
\]

annihilate the ordinary vector.  They also remove the arbitrary Jordan basis
shift \(|\psi\rangle\mapsto|\psi\rangle+\lambda|\phi\rangle\).  Their fixed-c
residues are

\[
(J_A,J_B)=(1500,-450),\qquad
\boxed{J_A/J_B=-10/3}.
\]

For a Kac-curve parameter derivative, the central-charge tangent changes the
residue to

\[
(J_A,J_B)=(-4820/3,1850/3),\qquad
J_A/J_B=-482/185.
\]

The ordinary Q4 descendant has both residues zero.  These are four distinct
outcomes:

```text
four-leg primary:      W = 0
ordinary thermal Q4:   W proportional to (4,-6,3), J_A=J_B=0
fixed-c Q4 Jordan top: W has residue J_A/J_B=-10/3
Kac-curve log tangent: W has residue J_A/J_B=-482/185
```

## Level-10 boundary

The thermal singular levels begin `2,10,16`.  The next ordinary spin-4 pair
after `(4,0)` is `(7,3)`, of total level 10 and dimension `45/4`.  It produces
the familiar relative correction `L^-6=N^-3`, but it remains in the ordinary
level-two null quotient.  It cannot create a nonzero inhomogeneous residue.
This prevents a level-10 ordinary correction from imitating the Jordan gate.

## Frozen lattice prediction

Calibrate one source-frozen lattice-to-Virasoro Ward projection on an
exact-critical control before looking at the target.  Save the three Q4 mode
contractions and the three parent-null contractions in one same-stream
covariance block.

The risky prediction for the leading local/torus matching-odd spin-4 row is:

1. its direct mode vector follows `4:-6:3`, rejecting the zero vector of the
   four-leg primary;
2. after removing the ordinary Q4 line, `J_A/J_B=-10/3` and the parent null
   jet has `d/b=0`, selecting a fixed-c Jordan lift.

Score the linear constraints before displaying ratios.  Transport the same
Ward covectors over norm 2 and the existing norm-5 handed pair; none of the
coefficients may be refit by cover.  This remains conditional until a lattice
observable implementing the mode projection is calibrated.
