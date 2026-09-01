# A regular homogeneous counterterm cannot remove the shared-line Q interaction

**Result.** Within the real, entry-regular family
`K_c=K2+c(Q)K0`, with unit coefficient of K2 and
`c(Q)=1+alpha(Q-1)+O((Q-1)^2)`, the same completion at two marked
vertices has the four-line activated Gram coefficient

\[
\boxed{g(\alpha)=\frac{13}{8}-\frac\alpha2+\frac{\alpha^2}{2}
 =\frac32+\frac12(\alpha-\tfrac12)^2\ \ge\ \frac32.}
\tag{1}
\]

Every higher Taylor coefficient of c drops out of this first Q
derivative. After summing the two interior holes' vacant/occupied states
in the physical contractible four-path exterior, the connected marked
coefficient is `g(alpha)/[(1+v_x)(1+v_y)]`, still strictly positive.
Thus no choice of this common regular counterterm makes the activated
interaction additive over the two marks. No alpha is fitted here.

The [kernel algebra](https://github.com/LightChainr/Matching-One/blob/7e46c74ce149d5a0a06d1085eb36eebb1bbe6bdb/notes/local-pair-two-insertion-algebra.md)
and [canonical two-site susceptibility](https://github.com/LightChainr/Matching-One/blob/2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb/notes/regular-pair-two-site-q-susceptibility.md)
are inputs. The [completed original-U activation](https://github.com/LightChainr/Matching-One/blob/2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb/notes/regular-pair-interaction-result.md)
is not recalculated. This note derives a counterterm-robust conditional
mechanism statement, not an unconditional global-U prediction or a
universal field identification.

## 1. Use a regular basis before taking the Q1 jet

All kernels are the prescribed C4 averages on the same four colour
ports. Ordinary colour inner products initially refer to integer Q>=4;
their rational closed-diagram expressions define the continuation.
Set

```text
delta=Q-1,
Kreg=K2+K0,
L=delta K0.
```

Both Kreg and L have entrywise regular equality-diagram expansions at
Q1. L is explicitly the averaged unequal-pair tensor divided by Q.
Both vanish on the one-colour set. Thus a finite closed diagram
containing one or more of these regular kernels has value zero at Q1,
although its Q derivative need not vanish.

The exact input pairings are

```text
<K2,K2> = Q(Q-3)(3Q^2-9Q+8)/[8(Q-1)(Q-2)],
<K2,K0> = (Q-3)/[4(Q-1)],
<K0,K0> = (2Q^2-4Q+3)/[2Q(Q-1)].
```

Their regular-basis combinations simplify to

\[
\begin{aligned}
\langle K_{\rm reg},K_{\rm reg}\rangle
 &=\frac{(Q-1)(3Q^3-12Q^2+20Q-24)}{8Q(Q-2)},\\
\langle K_{\rm reg},L\rangle
 &=\frac{(Q-1)(5Q-6)}{4Q},\\
\langle L,L\rangle
 &=\frac{(Q-1)(2Q^2-4Q+3)}{2Q}.
\end{aligned}
\tag{2}
\]

For example, the middle identity follows from
`<Kreg,K0>=(5Q-6)/(4Q)`; its apparent K0 pole has already cancelled
against the K2 cross term. Each expression in (2) has a factor Q-1.
The activated Gram matrix is therefore

\[
\boxed{\mathsf G
 =\left.\partial_Q
 \begin{pmatrix}
 \langle K_{\rm reg},K_{\rm reg}\rangle&
 \langle K_{\rm reg},L\rangle\\
 \langle L,K_{\rm reg}\rangle&\langle L,L\rangle
 \end{pmatrix}\right|_{1}
 =\begin{pmatrix}13/8&-1/4\\-1/4&1/2\end{pmatrix}.}
\tag{3}
\]

Its determinant is3/4 and both diagonal entries are positive, so this
real two-dimensional activated pairing is positive definite. This
specific algebraic result is not an assertion of reflection positivity
or a physical Hilbert norm for arbitrary noninteger-Q diagrams.

## 2. Higher counterterms cannot change this first activated pairing

Entry regularity in the stated family requires c(1)=1. Write the full
regular function as

```text
c(Q)=1+delta a(Q),          a(1)=alpha,
K_c=Kreg+a(Q)L.
```

For possibly different functions at two sites,

```text
H_xy(Q)=<K_cx,K_cy>
 =<Kreg,Kreg>+(a_x+a_y)<Kreg,L>+a_x a_y<L,L>.
```

Every basis pairing on the right vanishes at Q1. Terms involving
`a_x'(1)`, `a_y'(1)` or higher derivatives consequently contribute
only at order delta squared or later. Precisely,

\[
\boxed{H'_{xy}(1)
 =(1,\alpha_x)\mathsf G(1,\alpha_y)^T
 =\frac{13}{8}-\frac{\alpha_x+\alpha_y}{4}
       +\frac{\alpha_x\alpha_y}{2}.}
\tag{4}
\]

At Q1, `partial_Q` and `partial_logQ` agree. No extra factor from
`Q=m^2` is included in this convention.

With a homogeneous completion, `alpha_x=alpha_y=alpha`, equation(4)
is the squared activated norm in (3), giving (1). The lower bound does
not require choosing the minimizing alpha=1/2 as a new model. It holds
for every real alpha. At this derivative order a common alpha suffices;
the higher Taylor coefficients at the two sites need not agree.

The common-alpha assumption is essential to a universal positive lower
bound. A cross pairing of distinct vectors in a positive Gram metric
need not be positive. Indeed, with
`mu=(alpha_x+alpha_y)/2`, `d=(alpha_x-alpha_y)/2`,

```text
H_xy'(1)=3/2+(mu-1/2)^2/2-d^2/2.
```

For example `(alpha_x,alpha_y)=(0,13/2)` makes it zero, and `(0,7)`
makes it -1/8. These are exact boundary examples, not proposed fitted
counterterms. A claim allowing unrelated site-dependent completions
cannot use (1).

The scope also fixes the K2 coefficient to one. A vanishing overall
rescaling changes the coupling normalization and can remove the
insertion itself. Adding other regular tensor directions or an explicit
two-site counterterm is a larger family, not excluded by this
two-dimensional Gram statement.

## 3. The connected coefficient survives physical hole normalization

Use the [contractible four-path exterior](local-pair-two-insertion-obstruction.md#2-a-realizable-embedded-exterior-not-an-abstract-colour-wiring)
on the 17x17 torus: four disjoint occupied NN paths join the four ports
of holes x and y. All four choices of filling the holes have ambient
rank zero. The exterior reflection gives the same Frobenius pairing
because these C4-averaged kernels are also reflection invariant.

Fix the exterior occupations but sum both holes' vacant and occupied
states, with activities v_x,v_y>0. At a vacant hole insert
`1+lambda_x K_cx` or `1+lambda_y K_cy`. Remove the common exterior
activity and spectator-colour factor. The unmarked partition is

```text
Z0(Q)=Q^4+(v_x+v_y+v_x*v_y)Q,
Z0(1)=(1+v_x)(1+v_y).
```

When the other hole is vacant, one mark has four free exterior colours.
The all-free contraction of K2 is zero and that of K0 is Q(Q-1).
When the other hole is occupied, all four colours are identified and
both kernels vanish. The full conditional partition is therefore

\[
Z_{xy}(Q;\lambda_x,\lambda_y)
 =Z_0(Q)+\lambda_x C_x(Q)+\lambda_y C_y(Q)
              +\lambda_x\lambda_y H_{xy}(Q),
\tag{5}
\]

where

```text
C_j(Q)=c_j(Q)Q(Q-1),
C_j(1)=0,          C_j'(1)=1,
H_xy(1)=0.
```

For independent marks the connected coefficient at zero coupling is

```text
partial_lambda_x partial_lambda_y log Z_xy
 = H_xy/Z0 - C_x*C_y/Z0^2.
```

The product `C_x C_y` starts at order delta squared. Differentiating
the other denominator adds a factor H_xy(1)=0. Consequently

\[
\boxed{\left.\partial_{\log Q}\partial_{\lambda_x}
 \partial_{\lambda_y}\log Z_{xy}\right|_{Q=1,\lambda=0}
 =\frac{H'_{xy}(1)}{(1+v_x)(1+v_y)}.}
\tag{6}
\]

For a common completion it is at least
`3/[2(1+v_x)(1+v_y)]`, strictly positive at finite positive activities.
It is not a bound uniform as either activity tends to infinity. If
activities themselves vary regularly with Q, only their Q1 values
enter (6); denominator derivatives still multiply zero.

More explicitly, the Q-activated logarithm relative to the unmarked
conditional partition is exactly

\[
\left.\partial_{\log Q}\log\frac{Z_{xy}}{Z_0}\right|_{1}
 =\frac{\lambda_x+\lambda_y+
         H'_{xy}(1)\lambda_x\lambda_y}
        {(1+v_x)(1+v_y)}.
\tag{7}
\]

There are only two marked vertices, so this identity has no omitted
higher lambda terms. With both holes fixed vacant it reduces to the
same formula without the denominator, recovering13/8 for the canonical
completion alpha=0. With a common mark lambda, the second lambda
derivative gives twice the coefficient in (6). Site-average units
`lambda_j=epsilon_j/N` instead divide the independent-mark coefficient
by N squared.

## 4. The robust mechanism consequence

An additive Q-activated source at these two marks, or a product of
their independently closed one-mark weights, has zero mixed
`partial_lambda_x partial_lambda_y` coefficient at first order in Q-1:
the one-mark closed coefficients themselves vanish at Q1, and their
product is order delta squared. Equations (1) and (6) exclude that
replacement for **every common real regular counterterm in the declared
K2+c(Q)K0 family**. The shared four-colour-line sewing leaves an
unavoidable interacting activation.

This conclusion is stronger than the nonzero canonical coefficient13/8:
adjusting alpha cannot remove it, even though the completed one-mark
original-U response changes as

```text
W_alpha = W_canonical-alpha V_old.
```

That previously derived affine relation can change or even annul one
finite-size original-U coefficient; it cannot annul the homogeneous
conditional four-line interaction in (6). No old Q1 source was rescored
to reach this conclusion.

The statement is local and coefficientwise in the named physical
exterior. Occupation sums, other geometries, or different mark locations
can produce additional cancellations and different original-U
responses. Nothing here asserts a positive unconditional correlator,
a nonzero homogeneous global-U mixed derivative, a long-distance
amplitude, or a universal scaling field. Those require their own
specified sums and predictions. The result already removes one precise
simplification: a homogeneous regular scalar singlet counterterm does
not turn this two-mark Q-activated local colour interaction into an
additive first-insertion occupation source.
