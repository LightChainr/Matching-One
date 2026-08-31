# A removable twist-jet interface for the original global U

This note supplies a finite, directly substitutable interface at colour
`Q=1`. Its inputs are the **complete** local trace combinations `T,I,R`
of `closed-source-hypergraph-rc-twist-projection.md`, continued by that
note's occupation-sector prescription. It does not differentiate the
number of literal integer-colour twists. The new content is the jet order,
the two normalized trace directions, and their exact coupling functional
to the original moving-root observable.

The regular unlabelled selection identity in
`weak-q-paths-and-regular-selection.md` remains in force: its contracted
colour projector `P_[2]` is identically absent along the regular family,
including all regular Q derivatives. A finite trace/confluent contribution
has to enter the interface below through an actual occupation-compatible
trace combination. A nonzero formal projector trace alone is insufficient.

## 1. Remove the twist zero before differentiating

For each geometry separately, put `m=sqrt(Q)` and

```
T = L0 + L1 + L2,             I = L0,
R = (m-1)[(m+1)L0 + L1],      J = R/(m-1).
```

Here `J` means the removable twist quotient, not a source statistic or
the all-ones diagram. The positive real-Q continuation is fixed by

```
Lr(p,Q) = sum_(A:rank=r) [y Q^(-5/2)]^K Q^(B+CB),
y=p/(1-p).
```

The equivalent hypergraph chart is allowed; its thermal coordinate must
be transformed consistently. With `Z=T+(m-1)J`, the exact original
observers at every positive Q are

```
q = [T-J-m(m-1)I]/Z,
E = [T-J+m(m+1)I]/Z.                              (1)
```

Let `epsilon=Q-1`. Use *derivatives*, rather than Taylor coefficients,
in the definitions

```
T0=T|1, I0=I|1, T1=partial_Q T|1, I1=partial_Q I|1,
r1=partial_Q R|1, r2=partial_Q^2 R|1.
```

Since `m-1=epsilon/2-epsilon^2/8+O(epsilon^3)`, the removable jets are

```
J0 = 2 r1,
partial_logQ J|1 = r2 + r1/2.                     (2)
```

In the m chart these same quantities are `J0=R_m|1` and
`partial_logQ J|1=R_mm|1/4`. Equation (2) includes a non-optional `r1/2`
change-of-parameter term. **A first Q derivative of R does not determine
the first Q response of the original q/E.** R must be known through
quadratic order, or J supplied directly through linear order.

Define the baseline and its two normalized trace directions by

```
j = 2 r1/T0,                       i = I0/T0,
chi = (r2+r1/2-j T1)/T0 - i/2,
eta = (I1-i T1)/T0.                               (3)
```

All quantities in (2)-(3) are functions of the same thermal coordinate
p. At Q=1,

```
q=1-j,                 E=1-j+2i.                 (4)
```

The subtraction `i/2` in chi removes the explicit derivative of the
coefficient `(m+1)` in J. Consequently, for the *unprojected* occupation
law with sector weights Lr,

```
partial_logQ q_unprojected|1 = -chi,
partial_logQ E_unprojected|1 = -chi+2 eta.         (5)
```

This is a useful separation: chi and eta are genuine normalized trace
responses, with the purely kinematic multiplicity already removed.

## 2. The original rank projection contributes a known finite term

The original weights are proportional to `Lr Q^(-r/2)`. Thus its complete
Q tangents, denoted here by `u_q,u_E`, are

```
v = E-q^2,                  c = q(1-E),
u_q = -chi - v/2,
u_E = -chi + 2 eta - c/2.                         (6)
```

Indeed `v=Var(q)` and `c=Cov(E,r)` in each geometry, since configuration
q takes values -1,0,1 and E=q^2. Differentiating the rank factor gives
`-Cov(O,r)/2`. Equivalently, direct differentiation of (1) gives

```
u_q = -partial_logQ(J/T)|1 - i/2 - q j/2,
u_E = -partial_logQ(J/T)|1 + 2 eta + 3i/2 - E j/2,
```

which reduces to (6). In particular the rank correction survives even
when the unprojected sector probabilities have zero Q tangent. It cannot
be discarded on the ground that colour dependence was represented by
local traces.

For the specified closed-source path, (6) is also consistent with

```
partial_logQ <O>|1
  = Cov(O, CB+B-(5/2)K-r/2)
  = (1/2) Cov(O,Sstar).                          (7)
```

The harmless additive constant in Sstar drops from covariance. The K
term is a common thermal-coordinate direction and disappears from U,
while the rank term and local B term have no such automatic exemption.

## 3. A direct finite-jet functional, with separate normalizations

Use the original fixed angular weights

```
w_f=1/Delta,     w_s=-1/Delta,     P4 O=sum_g w_g O_g,
bar O=(O_f+O_s)/2,                A_N=N^(13/8)/2.
```

Each T0 in (3) normalizes **its own geometry**. Define, at Q=1,

```
calQ(p)=bar q(p),    Y(p)=P4 E(p),
D(p)=calQ'(p),      b(p)=Y'(p)/D(p).
```

The common root p0 solves `calQ(p0)=0`, equivalently `bar j(p0)=1`.
The original observable is `U=A_N b(p0)`. It does not require, or in
general have, a separate zero of q in each geometry. Assume D(p0)>0.

Here is the promised substitutable scalar:

```
Phi(p) = sum_g [(b/2-w_g) chi_g + 2 w_g eta_g]
         + (1/2)[b bar v - P4 c],

partial_logQ U|1 = (A_N/D(p0)) Phi'(p0).          (8)
```

Everything on the right is explicitly defined by the finite jets (2)-(4).
The last line of Phi is the mandatory rank-projection term. This supplies
both a value and a nonzero criterion: **the complete jet changes the
original U exactly when Phi'(p0) is nonzero**, not when an individual
trace, projector derivative, or Phi(p0) is nonzero.

For direct implementation without a moving-root solve at nearby Q,

```
Phi' = sum_g [(b/2-w_g) chi'_g + 2w_g eta'_g + (b'/2)chi_g]
        + (1/2)[b' bar v + b bar v' - P4 c'],
b' = (Y''-b calQ'')/D,
v'_g=E'_g-2q_g q'_g,
c'_g=q'_g(1-E_g)-q_g E'_g.                       (9)
```

Thermal differentiation acts on the normalized ratios in (3), including
T0 and the coefficients b. All p derivatives are taken **before** setting
p=p0. Replacing `bar q` by zero as a functional identity before taking
these derivatives would remove root-motion terms incorrectly.

For completeness, (1) provides its thermal derivatives directly: if n_O
is its numerator,

```
O'=(n'_O-O Z')/Z,
O''=(n''_O-O Z''-2O'Z')/Z.                       (10)
```

At Q=1 these use (4). Baseline T0,I0,r1 through two thermal derivatives,
and Q jets T1,I1,r2 through one thermal derivative, suffice for (8)-(9).
No third thermal derivative and no second response of U are needed.

To see explicitly why (8) includes the moving root, its tangent is
`p0_dot=-bar u_q/D`. Substitution into the derivative of `A_N Y'/D`
gives `(A_N/D) partial_p[P4 u_E-b bar u_q]`, precisely (8). This last
identity establishes the interface; the independently usable content
is the finite trace formula (3), (6), (8), and (9).

## 4. Net finite-collision coefficients: what can actually couple

Suppose a candidate finite trace/confluent contribution leaves the
Q=1 baseline unchanged, and supplies the additional complete jets

```
t_g = delta(partial_logQ T_g)|1,
z_g = delta(partial_logQ J_g)|1,
i_g = delta(partial_logQ I_g)|1.
```

The symbol `i_g` in this display is an *additional raw jet*, not the
baseline ratio i in (3). To avoid any ambiguity, write the latter as
`i0_g` and the baseline j as `j0_g` in the next formula. The exact net
loading functional is

```
Lambda(p) = sum_g (1/T0_g) {
  (b/2-w_g) z_g + 2w_g i_g
  - [(b/2-w_g)j0_g+2w_g i0_g] t_g },

delta(partial_logQ U)|1 = (A_N/D) Lambda'(p0).    (11)
```

Thus the required residue/finite-part data are a **three-component
landing vector in each geometry**, not a lone scalar colour trace.
Equation (11) quantifies every omitted term: insert its missing vector
and differentiate. For example, omitting r2 in (2), while retaining all
other inputs, misses exactly

```
(A_N/D) partial_p sum_g [(b/2-w_g) r2_g/T0_g]|p0. (12)
```

The residue is invisible to U if and only if the full Lambda derivative
vanishes. This includes several different possibilities:

- A geometry-specific common partition factor has
  `(t,z,i)=a_g(p)(T0,J0,I0)` and Lambda is identically zero, separately
  in each geometry. Normalization must precede pooling.
- A common thermal reparametrization gives normalized tangents
  `delta u_q=a(p)q'`, `delta u_E=a(p)E'` for both geometries, and again
  Lambda is identically zero, including the derivative of a(p).
- Lambda may be nonzero but locally constant at p0, so this particular
  slope observable has zero response even though a transverse profile
  deformation exists.

In the hypergraph chart the common prefactor `Q^(1/2-N)` is of the first
kind. The `3v partial_v/2` coordinate term along `v=y Q^(3/2)` is of the
second kind for the root/slope U. Either can be eliminated only as a
consistent normalization/coordinate transformation of all observers;
one cannot delete selected terms of an unnormalized trace derivative.

For an explicit trace-only algebraic example, take

```
f(Q)=(Q-1) tr P_[2] = (Q-1) Q(Q-3)/2,
f(1)=0,       partial_logQ f|1=-1.
```

If its actual landing is `delta(T,J,I)=f(Q)(a_T,a_J,a_I)`, then (11)
uses `(t,z,i)=-(a_T,a_J,a_I)`, in each geometry. The number -1 supplies
only this sign and scale. A nonzero response still requires the thermal
derivative of the complete weighted expression (11). The landing must
be justified in the specified occupation completion; this example does
not itself establish such a contribution to the physical model.

There is a further removable-zero trap. If the same f is inserted into
**R**, rather than J, as `delta R=f(Q)a_R(p)`, then

```
delta J0=-2a_R,
delta q|1=delta E|1=2a_R/T0                     (13)
```

when T0 and I0 are unchanged. Although the added R trace vanishes at
Q=1, its removable quotient changes the baseline observable. It is not
a tangent-only activation of the original model. A baseline-vanishing
J contribution instead requires its R landing to start at order
`(Q-1)^2`, or an explicit cancellation of the linear R coefficient in
the complete trace combination. In that baseline-preserving case the
additional z in (11) is simply `delta R_QQ|1`.

## 5. Pole cancellation and the regular selection boundary

The colour-projector symbol `P_[2]` here is not the ambient rank-2 event.
In a collision calculation write each **fully contracted** contribution
as a Laurent series in epsilon=Q-1. The opposite poles

```
epsilon P_[0] -> +2 J_diagram,
epsilon P_[2] -> -2 J_diagram
```

must cancel in the complete physical T,I,J combinations. Their combined
finite coefficient fixes the baseline; their combined epsilon coefficient
fixes the jets entering (11). When a pole multiplies a Q-dependent
endpoint, operator, measure, or amplitude, that epsilon coefficient can
contain the multiplier's *second* Taylor coefficient. Keeping only the
projector derivative does not recover it. If R is assembled directly,
its zero constant term and its quadratic coefficient must additionally
be retained, as (2) shows; working with the already-regular J avoids this
extra external order.

For the regular invariant endpoint, the complete contraction satisfies
`ell(Q) P_[2](Q)=0` identically. In particular

```
ell' P_[2] + ell P_[2]' = 0,
```

with its removable contracted meaning at Q=1. This produces a zero
landing vector for that channel, so (11) is zero. Computing only
`ell P_[2]'` discards its exact cancelling endpoint term. Neither the
new twist interface nor the nonzero polynomial `tr P_[2]` reopens an
ordinary regular unlabelled four-leg Q activation.

A genuinely nonregular trace/confluence calculation can be tested
without guessing a field identity: form its complete finite jets in
the fixed occupation family, insert their landing vector into (11),
and ask whether Lambda'(p0) is nonzero. If the contribution changes the
Q=1 baseline itself, the baseline, root, and all coefficients must first
be recomputed from (1)-(4); the baseline-preserving variation (11)
cannot be used across two different baseline models.

## Scientific card

The interface converts the apparent `R/(sqrt(Q)-1)` singularity into a
finite two-direction normalized trace response plus a fixed topological
projection term. It fixes both the necessary quadratic R jet and the
exact net coefficient a finite collision must supply to alter original U.
It does not calculate a new trace spectrum, establish a surviving
collision, infer a continuum field, or supply independent stochastic
evidence. No enumeration, simulation, coupling grid, or new source was
used. The prescribed real-Q occupation completion and a simple pooled
root are essential assumptions.
