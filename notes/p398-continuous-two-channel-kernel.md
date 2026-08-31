# P398: a two-channel kernel with amplitude-free and metric-free fingerprints

Parent: `b35e100a3903c706dceba57c4667386eb4510ac3`. This is the continuous
distance limit of the same positive width-four square-bond cylinder, not a
new field identification. Zero Monte Carlo and no server are used.

**Main result:** the self-dual limit gives a fully determined two-channel
correlation kernel, including its amplitudes. A one-parameter departure from
self-duality changes a measurable mass ratio and mixing invariant; it is not
just a metric change. One parameter-free relation links those two invariants.

## The surviving parameter is a rate ratio, not an extra time unit

Set `h=epsilon`, `v=1-kappa epsilon`, with `kappa>0`, and hold
`s=epsilon*d` fixed. Expanding the already-derived polynomial B gives

\[
B=I-\epsilon M_\kappa+O(\epsilon^2),\qquad
M_\kappa=\begin{pmatrix}
3+2\kappa&\kappa(1+i)\\1-i&2+3\kappa
\end{pmatrix}.
\]

The full 14-state probability generator is
`G=sum_i(J_i-I)+kappa sum_i(D_i-I)`: horizontal joins have rate one,
vertical disconnections rate kappa. Only kappa=1 is the self-dual line.
For kappa!=1, the O(epsilon) displacement from `h+v=1` survives as a finite
join/detach-rate detuning in this limit.

More generally, `h=a epsilon`, `v=1-b epsilon` produces `a M_(b/a)`.
Changing a rescales distance/masses; changing kappa=b/a changes the physical
generator after that common scale has been removed.

Let `R=sqrt(kappa^2+6kappa+1)` and `m0=5(1+kappa)/2`. Then

\[
m_{s,f}=m_0\mp R/2.
\]

In the fixed phase/metric gauge

\[
D=\operatorname{diag}(1,e^{-i\pi/4}/\sqrt\kappa),\quad
H_\kappa=D^{-1}M_\kappa D
=\begin{pmatrix}3+2\kappa&\sqrt{2\kappa}\\
\sqrt{2\kappa}&2+3\kappa\end{pmatrix},
\]

the mixing angle satisfies
`theta=.5 atan2(2 sqrt(2 kappa),1-kappa)`. This angle refers to these explicitly
normalized AP/landing coordinates; it is not invariant under arbitrary mixing
of the two observables. Both masses remain positive and distinct.

There is an exact exchange relation
`H(kappa)=kappa sigma_x H(1/kappa) sigma_x`. Hence mass ratios cannot separate
kappa and 1/kappa without keeping the named channel orientation.

## Full self-dual kernel: no free amplitudes remain

The exact stationary law of G at kappa=1 gives

\[
C(0)=\frac17\begin{pmatrix}6&-4+4i\\-4-4i&6\end{pmatrix},
\qquad C_{ab}(s)=E[O_a(0)\overline{O_b(s)}],\quad O=(A,L).
\]

The two charged means vanish. Define

\[
\psi_s=(A-e^{-i\pi/4}L)/\sqrt2,\qquad
\psi_f=(A+e^{-i\pi/4}L)/\sqrt2.
\]

Their complete connected kernel is diagonal:

\[
\boxed{C_\psi(s)=\operatorname{diag}\left(
\frac{6+4\sqrt2}{7}e^{-(5-\sqrt2)s},
\frac{6-4\sqrt2}{7}e^{-(5+\sqrt2)s}\right).}
\]

The small fast-channel amplitude is nonzero. Thus neither an unknown lattice
amplitude nor a fitted two-exponential model is needed to state this finite
control. Stationary covariance and propagation diagonalize together at the
self-dual point. The script also finds nonzero lag-covariance skew at
kappa=1/4,1/2,2,4; Hermitian similarity of a transfer matrix should not be
confused with reversibility in an arbitrary stationary readout metric.

## Remove source amplitudes before comparing propagation

For any kappa let `C0` be the stationary two-point matrix and define

\[
U(s)=C_0^{-1}C(s)=e^{-\overline{M_\kappa}s}.
\]

This convention follows the parent row of functions `f=(A,L)`. Under any
fixed invertible readout change `f->fD`, U changes by
`U->conjugate(D)^-1 U conjugate(D)`; its trace, determinant and spectrum are
therefore independent of field normalization or even of a change of basis.

The complete normalized kernel is

\[
U(s)=e^{-m_0s}\left[I\cosh(Rs/2)
-\frac{2(\overline M-m_0I)}R\sinh(Rs/2)\right].
\]

Consequently the following quantity removes both unknown amplitudes and the
common distance unit, without fitting a slope:

\[
\boxed{I_m=\frac{2\operatorname{acosh}
 [\operatorname{tr}U/(2\sqrt{\det U})]}{-\log\det U}
=\frac{R}{5(1+\kappa)}.}
\]

This is invariant under general similarity and a common metric rescaling.
It determines `m_fast/m_slow=(1+I_m)/(1-I_m)` and is constant for every s>0.

## A channel cross-ratio and a linked, parameter-free fingerprint

Keeping the named AP/landing channels, the ratio

\[
X(s)=\frac{U_{12}U_{21}}{U_{11}U_{22}}
=\frac{a\sinh^2(Rs/2)}{1+a\sinh^2(Rs/2)},\quad
a=\frac{8\kappa}{R^2},
\]

is unchanged by arbitrary nonzero complex rescaling of either readout. It is
not invariant under arbitrary mixing of A and L. Its time-independent form is

\[
\boxed{I_c=\frac{4U_{12}U_{21}}{(\operatorname{tr}U)^2-4\det U}
=\frac{8\kappa}{R^2}=\sin^2(2\theta).}
\]

Eliminating kappa gives a genuinely model-specific no-fit relation:

\[
\boxed{25 I_m^2(2-I_c)=2.}
\]

This relation, constancy with separation, and `U(s+t)=U(s)U(t)` together give
a sharper finite-model challenge than a generic two-mass fit. The signed
coordinate `(U11-U22)/sqrt[(tr U)^2-4det U]=(kappa-1)/R` preserves which named
channel is which and resolves the kappa-versus-inverse ambiguity. The separate
identity `I_c+signed_coordinate^2=1` is generic 2x2 algebra, not extra physics.

At the self-dual point the parameter-free targets reduce to

```text
I_m = sqrt(2)/5,    I_c = 1,
X(s) = tanh(sqrt(2)*s)^2.
```

Even without the C0 inverse, in the fixed Hermitian gauge the **raw**
cross-ratio is `C12 C21/(C11 C22)=tanh(sqrt(2)s+acosh(3))^2`. Its offset is
fixed by the exact stationary amplitudes, not fitted; subtracting the
`atanh(sqrt(X_raw))` value at s=0 removes that amplitude offset.

## Metric versus true relative splitting

| kappa | slow mass | fast mass | I_m | I_c |
|---:|---:|---:|---:|---:|
| 1/4 | 2.32460947 | 3.92539053 | .2561249695 | .7804878049 |
| 1/2 | 2.71922359 | 4.78077641 | .2748737084 | .9411764706 |
| 1 | 3.58578644 | 6.41421356 | .2828427125 | 1 |
| 2 | 5.43844719 | 9.56155281 | .2748737084 | .9411764706 |
| 4 | 9.29843788 | 15.70156212 | .2561249695 | .7804878049 |

`1/5<I_m<=sqrt(2)/5`; the fast/slow mass ratio lies strictly above 3/2 and at
most `(5+sqrt(2))/(5-sqrt(2))`. Self-duality maximizes the relative splitting
and gives equal-angle mixing. Varying kappa changes these dimensionless
numbers, whereas varying only the common update rate does not.

These invariants are **metric/gauge independent within this model**, not
universal continuum CFT data. Width is still four; no site-Matching, thermal
Q4, field multiplicity or Jordan identification follows. The continuous
probability generator and a Q-derivative/Jantzen generator remain different
objects.

## Reproduction

```sh
python3 scripts/p398_continuous_kernel.py
python3 -m unittest discover -s tests -p 'test_p398_continuous_kernel.py'
```

The h jet and stationary distributions are rational-exact; hyperbolic
functions are numerical displays of the derived closed kernel. Two focused
checks suffice; no full repository suite or stochastic production is used.
