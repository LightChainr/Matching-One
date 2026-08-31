# P398: positive anisotropy cannot collide the two modes, but can reveal the fast one

Parent: `e38fe7634354b0cb2201fa55fd9b4d37ccedeef2`.

**New result:** the complete two-parameter positive square-bond family has two
distinct ordinary charge-one propagation eigenvalues everywhere in its open
probability domain. On the self-dual line, anisotropy greatly widens the
row-distance window in which the faster mode is visible. The endpoint
eigenvalue collision is scalar, not Jordan; a genuine exceptional point exists
only after analytic continuation outside the probability domain in the
example below.

No new Monte Carlo, server, terminal mark, or formal jet is used. The row word,
readout slice and AP/landing definitions are unchanged from the parent.

## 1. Exact two-parameter propagation

Let `h,v` be horizontal/vertical occupation probabilities. Columns are source
states. Vertical edges are added first; horizontal edges of the new row follow:

\[
T(h,v)=\prod_i[(1-h)I+hJ_i]\prod_i[vI+(1-v)D_i(1)].
\]

The observables are the same `f=(A,L)` of the width-four parent. The backward
action on this row of functions is `T^transpose f=f B`. Put
`q=(1-h)^2 v^2`. Then

\[
\boxed{B(h,v)=q\begin{pmatrix}
1+h-2hv&-(1+i)(1-v)\\
-(1-i)hv&v
\end{pmatrix}.}
\]

The script derives the entries directly by summing each row's 256 bond masks
on all 14 inherited connectivity states. It expands the exact Bernstein
weights and solves the four fixed real readout coefficients. Thus the formula
is a polynomial identity, not an interpolation or a fitted closure ansatz.

Writing `t=1+h+v-2hv`,

\[
\operatorname{tr}B=qt,\quad \det B=(1-h)^5v^5,
\]

\[
\boxed{\Delta=q^2\{(1+h-v-2hv)^2+8hv(1-v)\}.}
\]

The second summand is strictly positive when `0<h,v<1`. Therefore the roots

\[
\lambda_\pm=\frac q2\left[t\pm\sqrt{t^2-4v(1-h)}\right]
\]

are distinct. Trace and determinant are positive, so both roots are positive.
This is stronger than a numerical failure to find a collision.

There is also a constructive reason. With

\[
S=\operatorname{diag}(1,\sqrt{hv/(1-v)}),
\]

`S^-1 (B/q) S` is Hermitian, with conjugate off-diagonal entries
`-(1±i)sqrt(hv(1-v))`. This **specific sector** is diagonally similar to a
Hermitian matrix throughout the positive domain. We are not asserting that
positivity forbids Jordan blocks in arbitrary stochastic matrices.

## 2. The collision boundary and a real signed Jordan point

In the closed probability square the discriminant vanishes only on `h=1` or
`v=0`, where `B=0`, and at `(h,v)=(0,1)`, where `B=I`. Every such collision
is scalar/semisimple. No rank-two Jordan point lies even on this boundary.

For fixed `0<v<1`, the nontrivial analytically continued collision curves are

\[
h=-\frac{1-v}{(1+\sqrt{2v})^2},\qquad
h=-\frac{1-v}{(1-\sqrt{2v})^2}.
\]

Both require negative h. At `v=1/2` the second branch is at infinity; the
finite exceptional point is `h=-1/8`. There

\[
\lambda=243/1024,\qquad
B-\lambda I=\frac{81}{4096}
\begin{pmatrix}4&-8(1+i)\\1-i&-4\end{pmatrix},
\]

whose square is zero but which is itself nonzero. This is an exact size-two
Jordan block **with a negative horizontal bond weight**, not a physical
percolation result. It distinguishes a genuine exceptional point from the
scalar collisions of the positive family.

## 3. Exact self-dual line and the strongly anisotropic limit

On `v=1-h`,

\[
\boxed{\lambda_\pm=(1-h)^4[1-h+h^2\pm h\sqrt{2-2h+h^2}].}
\]

Equivalently, let

\[
\gamma=\operatorname{acosh}\left(1+\frac{h^2}{1-h}\right),\qquad
\lambda_\pm=(1-h)^5e^{\pm\gamma}.
\]

The fast/slow ratio is `exp(-2 gamma)`. It decreases strictly from one to zero
as h runs from zero to one, since
`d[h^2/(1-h)]/dh=h(2-h)/(1-h)^2>0`.

Although the two row eigenvalues approach one as `h->0`, this is a fine-time
limit, not a Jordan limit:

\[
\frac{B-I}{h}\longrightarrow
\begin{pmatrix}-5&-1-i\\-1+i&-5\end{pmatrix},
\]

with distinct generator eigenvalues `-5±sqrt(2)`. At fixed rescaled distance
`s=h*d`, the modes converge to `exp[-(5∓sqrt(2))s]`. Their two masses remain
separated. The symmetrizer stays well conditioned on this limit (`S->I`).

## 4. An explicit fast-mode readout, not a two-exponential fit

The same existing A/L data can form the parameter-fixed eigenobservables

\[
\Psi_\pm=A+\frac{h\mp\sqrt{1+(1-h)^2}}2(1-i)L.
\]

They obey the exact stationary identity

\[
\frac{\langle\Psi_\pm(0)\overline{\Psi_\pm(d)}\rangle}
{\langle|\Psi_\pm|^2\rangle}=\lambda_\pm^d.
\]

Both means vanish by cyclic frontier charge. Full-support stationary measure
and independent readout columns ensure both variances are nonzero. The table
also explicitly checks their numerical sizes using exact rational stationary
probabilities; only square roots and displayed decimals are floating point.

| h | v | lambda+ | lambda- | fast/slow | ratio at d=8 | Var(Psi-) | fast-mode C(8) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| .05 | .95 | .83199009 | .71964431 | .86496741 | .31332710 | .09887654 | .00711273 |
| .10 | .90 | .68532023 | .50878177 | .74240005 | .09227923 | .09930013 | .000445864 |
| .20 | .80 | .44897279 | .23915521 | .53267196 | .00648155 | .09856171 | 1.05474e-6 |
| .50 | .50 | .08181356 | .01193644 | .14589803 | 2.05303e-7 | .07546804 | 3.11000e-17 |
| .80 | .20 | .00264935 | .0000386510 | .01458887 | 2.05198e-15 | .01643510 | 8.18575e-38 |

The improvement is not just a ratio of two vanishing amplitudes: the fixed
fast-mode variance stays near 0.1 at h=.05/.10/.20. At h=.05, its d=8
correlation is directly macroscopic on this bounded-observable scale. The
relative-decay factor remains above 10% through d=15, versus only d=1 at
h=.5. These are exact-model predictions, not power or sample-size claims.

Thus a future physical readout of this finite model should use a vertical-
strong anisotropy such as `h=.05,v=.95` and the explicit Psi- combination,
rather than try to detect its isotropic d=8 tail. No new sampling is performed
or authorized here. More rows are required for a matched rescaled distance;
we have not measured a variance-per-CPU gain.

## 5. Scientific boundary and reproduction

The result rules out a positive-weight ordinary-to-Jordan collision for the
**same width-four square-bond charge-one propagation family**. It does not
rule out Jordan behavior for a different width, representation, observable,
dilation operator or Q-lift, and it does not identify a site-Matching field.
The width remains four, and cyclic frontier charge is not a continuum spin.

```sh
python3 scripts/p398_anisotropic_cylinder.py
python3 -m unittest discover -s tests -p 'test_p398_anisotropic_cylinder.py'
```

The result contains the full polynomial coefficients, factorizations,
nonphysical Jordan witness, exact stationary C0 at five declared points and
input hashes. Two focused arithmetic tests; no full repository suite.
