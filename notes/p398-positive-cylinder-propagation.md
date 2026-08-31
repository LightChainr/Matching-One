# P398: two ordinary propagating responses in a positive width-four cylinder

**Result:** the existing rooted AP and landing charge-one readouts have two
distinct, nonzero propagation eigenvalues in a concrete positive probability
measure. Their finite-separation response is not one common ray. This is a
physical finite-lattice response calculation, not another formal jet closure.

Parent: `afc619c35188af718d450c19935f4d8dd8d39ff1`. No Monte Carlo, remote
machine, fitted amplitude, new mark, or continuum identification is used.

## 1. Probability measure and the two already-defined observables

Use independent square **bond** percolation, Q=1, with horizontal and vertical
occupation probabilities both 1/2 on a width-four cylinder. A state is the
partition of the four frontier vertices by connectivity through the past.
The existing circular noncrossing basis has 14 states. Columns of matrices
are source states. One full row is

\[
T=HV,\quad V=\prod_{i=0}^3(I+D_i(1))/2,\quad
H=\prod_{i=0}^3(I+J_i)/2.
\]

The vertical layer is applied first, then the new horizontal layer. Readouts
are taken after that horizontal layer. The generators within a layer commute.
This is exactly the sum over its eight independent bond bits, not a signed
Gram or a derivative accumulator. The small focused check compares all
14 columns with the direct 256-bond-mask sum.

Solving `T pi=pi`, `sum pi=1` gives, in the committed codec order,

```text
pi = (656,88,88,64,52,88,4,64,88,52,52,4,52,41)/1393.
```

The stationary probability is unique and strictly positive. There is a
positive self-loop at every state, so this is the stationary limit of the
finite-strip construction, not an arbitrary weighting of connectivity states.

The old Gram-derivative columns simplify to bounded configuration functions:

\[
A=(1_{0\sim1}-1_{2\sim3})+i(1_{1\sim2}-1_{3\sim0}),
\]
\[
L=(1_{0\text{ singleton}}-1_{2\text{ singleton}})
+i(1_{1\text{ singleton}}-1_{3\text{ singleton}}).
\]

For AP, joining its adjacent pair changes the block count by minus one
exactly when those sites were disconnected. For landing, wiring the other
three vertices leaves one extra block exactly when the selected vertex was
singleton. Constants cancel in the charge-one differences. These identities
are checked against the inherited columns on all 14 states.

Both means vanish. Their charge is under cyclic permutation of the four
frontier labels; it is **not a continuum spin assignment**.

## 2. The entire charge-one sector closes, and is diagonalizable

Write `f=(A,L)` as a row of functions. The backward transfer satisfies

\[
T^t f=fB,\qquad
B=\begin{pmatrix}
1/16&(-1-i)/32\\
(-1+i)/64&1/32
\end{pmatrix}.
\]

This is not a fitted two-state approximation. The four real columns
`Re A, Im A, Re L, Im L` are independent and exhaust `ker(R^2+I)`, whose exact
real dimension is four. Hence each of the two conjugate cyclic-charge
sectors has complex dimension two.

\[
\det(\lambda I-B)=\lambda^2-\frac3{32}\lambda+\frac1{1024},
\quad \lambda_\pm=\frac{3\pm\sqrt5}{64}.
\]

The eigenvalues are approximately 0.08181356215 and 0.01193643785, with decay
lengths 0.3994707408 and 0.2258274545 rows. Explicit eigenfunctions are

\[
\Psi_\pm=A+\frac{(1\mp\sqrt5)(1-i)}4L.
\]

Thus this positive realization has two ordinary exponential responses,
not a Jordan block. The fast/slow ratio is
`(7-3 sqrt(5))/2 = 0.145898...`: an apparent single response at long separation
is expected even though the exact finite-separation rank remains two.

## 3. Neutral two-point matrix and ordered cross term

Define `C_ab(d)=E[O_a(X0) conjugate(O_b(Xd))]_connected`. All means vanish.
With the convention above,

\[
C(d)=C(0)\overline B^d,\qquad
C(d+2)=\frac3{32}C(d+1)-\frac1{1024}C(d).
\]

The definition baseline is

\[
C(0)=\frac1{1393}
\begin{pmatrix}912&560(-1+i)\\560(-1-i)&768\end{pmatrix}.
\]

Its positive determinant is not the headline: equal-time positive
definiteness follows from full support and independent readouts. The new
propagation statement is

\[
C(1)=\frac1{1393}
\begin{pmatrix}149/2&46(-1+i)\\47(-1-i)&59\end{pmatrix},
\quad
\det C(d)=\frac{73216}{1940449}\,1024^{-d}>0
\]

for **every finite integer d>=1**, not only the explicitly tabulated d=1..8.
An invertible C(0) and the two distinct nonzero eigenvalues also rule out
all C(d) being scalar multiples of one fixed matrix.

| d | C_AA | Re C_AL | Re C_LA | C_LL | det C |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.0534816942 | -0.0330222541 | -0.0337401292 | 0.0423546303 | 3.68471421e-05 |
| 2 | 0.00437455133 | -0.00270324838 | -0.00277054917 | 0.00343234027 | 3.59835372e-08 |
| 3 | 0.000357885970 | -0.000221181241 | -0.000226789640 | 0.000280419957 | 3.51401730e-11 |
| 4 | 2.92797869e-05 | -1.80958503e-05 | -1.85559143e-05 | 2.29374762e-05 | 3.43165752e-14 |
| 5 | 2.39548201e-06 | -1.48048866e-06 | -1.51814271e-06 | 1.87654078e-06 | 3.35122805e-17 |
| 6 | 1.95982896e-07 | -1.21124083e-07 | -1.24204869e-07 | 1.53525819e-07 | 3.27268364e-20 |
| 7 | 1.60340586e-08 | -9.90959311e-09 | -1.01616452e-08 | 1.25604862e-08 | 3.19598012e-23 |
| 8 | 1.31180345e-09 | -8.10739116e-10 | -8.31360423e-10 | 1.02761802e-09 | 3.12107434e-26 |

`Im C_AL=-Re C_AL`, while `Im C_LA=Re C_LA`. All exact fractions, stationary
probabilities, the full transfer and readouts are in the JSON artifact.

The directional difference is nonzero already at one row:

\[
C_{AL}(1)-\overline{C_{LA}(1)}=(1-i)/1393.
\]

The past-connectivity frontier process is not reversible under these two
readouts. This does not violate reflection invariance of independent bonds:
reversing the cylinder turns a past-connected readout into a future-connected
one. It is not evidence of intrinsic temporal memory or a noncommuting
microscopic spacetime action.

## 4. What changed, and the finite boundary

The formal charge-one completion now has an explicit positive-measure
realization with two nonzero propagating mixtures. It does not require a
Jordan mechanism to make the two old readouts physically visible. The common
single-ray model fails exactly here, although long-distance data would rapidly
look approximately rank one.

This neither eliminates a Jordan interpretation in a different physical
family nor identifies site-Matching, thermal Q4, a continuum scaling dimension,
or a field multiplicity. Width is fixed at four; the geometry and bond
measure are explicit. The 23-state signed retained-response module was not
used as a stochastic transfer.

If a subsequent scientific comparison is desired, the useful question is
whether the two mixtures and their propagation survive a change of width or
microscopic observer. More projected Q-jets would not answer it. No such
extension is performed here.

## Reproduction

```sh
python3 scripts/p398_positive_cylinder.py
python3 -m unittest discover -s tests -p 'test_p398_positive_cylinder.py'
```

Python standard library only; three focused tests pass. No broad repository
suite was run. The JSON pins the protocol, runner, and inherited input hashes.
