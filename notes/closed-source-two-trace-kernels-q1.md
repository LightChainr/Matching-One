# Two actual closure kernels give an explicit Q1 trace landing

**Finite algebra.** The two rank-one connection types relevant to the
N25 colour trace have the same `[Q-2,2]` trace numerator but different
colour degeneracies. Their relative weights are therefore

```text
type A: two essential components, |first winding|=1,
        beta_A(Q)=(Q-3)/(2Q);
type B: one essential component, |first winding|=2,
        beta_B(Q)=(Q-3)/2.                                    (1)
```

These are actual integer-colour closure coefficients for every Q>=4,
continued by the stated rational formulas. They are not an interpolation
of the old Q4 result or an arbitrarily appended factor Q-1. The
[N25 packing theorem](n25-winding-packing-and-pair-continuation.md)
proves that these types exhaust the trace on both fixed geometries.
The [stable character derivation](closed-source-stable-colour-character-continuation.md)
also establishes the all-integer validity range and the finite-Q aliases
that this packing excludes. The algebra below supplies their complete
finite landing.

## 1. Identity versus exchange on the same two-colour carrier

Write V for the Q-dimensional permutation representation. Let
`d2(Q)=Q(Q-3)/2`. In `V tensor V`, `[Q-2,2]` occurs once, in the
symmetric unequal-colour summand. The interchange S of the two colour
ports acts as +1 on that summand.

Type A has two independent once-winding essential components. Its seam
class function is `Fix(pi)^2=Tr[(rho(pi) tensor rho(pi)) I]`.
The identity-seam colour count is Q², while its `[Q-2,2]` contraction
is d2(Q). Hence its relative coefficient is d2(Q)/Q².

Type B has one primitive component crossing the selected period twice.
Its seam class function is instead
`Fix(pi²)=Tr[(rho(pi) tensor rho(pi)) S]`. The identity-seam count
is Q. Because S acts as +1 in the same `[Q-2,2]` summand, its contraction
is again d2(Q). Its relative coefficient is d2(Q)/Q. Reversing winding
replaces pi by pi^-1 and gives the same class function.

Contractible spectator components, the occupied-edge weight and the
stipulated `Q^-r/2` multiply both the original and filtered weights by
the same factors. None is omitted in (1). The argument is a character
contraction of a fixed connection pattern, not an assertion that two
chosen ports span every state of the full transfer matrix.

The explicit endpoint values are

| Coefficient | Q4 value | Q1 value | Q derivative at1 |
|---|---:|---:|---:|
| beta_A | 1/8 | -1 | 3/2 |
| beta_B | 1/2 | -1 | 1/2 |

Thus the Q4 trace treats the two types differently by a factor4. At Q1
their relative trace weights agree, while their Q jets still differ.
The negative Q1 coefficient is a finite analytic trace weight, not a
negative number of physical colours or a probability.

## 2. The specified positive occupation family and its trace polynomial

Keep the original reduced partition convention throughout this note:

```text
Z(y,Q)=sum_A y^K Q^(-(K+g)/2),
g=4K-2B-2C_B+r,           y=p/(1-p).                           (2)
```

This is exactly the closed-source family with the common factor
`Q^(N+1/2)` removed. The original positive occupation weights are
well-defined for all real Q>0. The signed trace is the separate function

```text
F(y,Q)=sum_(A in type A or B) y^K Q^(-(K+g)/2) beta_type(Q).    (3)
```

No change of base measure or thermal root is implied by decomposing its
partition in this way. The trace's raw q and E numerators are zero for
all Q, because both types have topology rank1. At Q1 the two necessary
raw packets are explicitly

```text
F0(y)=-sum_(A or B) y^K,
F1(y)=partial_logQ F|1
     =sum_A y^K (K+g+3)/2 +sum_B y^K (K+g+1)/2.               (4)
```

The fixed local source score from (2) is `-(K+g)/2`. Multiplying it by
beta_1=-1 and adding the explicit beta derivative gives (4). Discarding
the measure derivative would keep only the constants3/2 and1/2 and
would be a different packet.

The primary physical perturbation considered here is the trace
coefficient epsilon at Q1:

```text
w_epsilon(A)=w(A)[1+epsilon beta_1(A)],
beta_1=-1_typeA-1_typeB.                                     (5)
```

It is positive in an open interval around epsilon=0. Its original-U
response is the established rank1 transmission functional applied to
F0/Z. In particular this trace is **already nonzero at Q1** as a partition
component; it is not a baseline-vanishing Q-activated endpoint.

## 3. Its landing in the removable twist interface is explicit

Use the earlier unprojected occupation sectors Lr and the removable
combination `J=(sqrt(Q)+1)L0+L1`, with `T=L0+L1+L2`, `I=L0` and
`R=(sqrt(Q)-1)J`. Because (3) is a projected rank1 contribution,
varying its coefficient changes the unprojected rank1 sector by
`delta L1=sqrt(Q) F`. Therefore the complete finite landing is

```text
delta(T,J,I)=(m F,m F,0),
delta R=m(m-1) F,             m=sqrt(Q).                     (6)
```

The observer normalizer in those twist coordinates is
`D_twist=T+R=Q Z`; its variation is QF. Its q/E numerator variations
are exactly zero. Formula (6) makes the trace-to-observable map concrete;
there is no free landing amplitude to fit.

At Q1,

```text
delta T0=delta J0=F0,
delta T_Q=delta J_Q=F1+F0/2,
delta R_Q=F0/2,
delta R_QQ=F1+F0/4.                                         (7)
```

Thus `2 delta R_Q=F0` and
`delta R_QQ+delta R_Q/2=F1+F0/2`, exactly as required by the previously
derived removable-quotient rule. In particular the linear R zero carries
a genuine change in the Q1 trace-coefficient baseline. It must not be
mistaken for a tangent-only activation with unchanged q/E.

## 4. Keep a physical trace response separate from attribution of a Q derivative

Let `T_U[F]` denote the existing linear rank1 normalization transmission
functional, using F/Z and the full moving pooled root. The primary
response is `V_trace=T_U[F0]`. It is invariant under a common nonzero
partition prefactor.

One can also compute `J_trace=T_U[F1]` as the trace packet's additive
contribution to the Q derivative in the fixed reduced convention (2).
This is useful bookkeeping, but its attribution is not invariant under
redistributing a common Q-dependent prefactor among the components.
For example replacing both Z and F by Q^c times themselves gives

```text
V_trace -> V_trace,
J_trace -> J_trace+c V_trace.                                (8)
```

The full physical log-Q response remains unchanged; the other partition
components compensate. Restoring the original source prefactor uses
c=N+1/2. Thus J_trace is neither a probability share nor a uniquely
defined amount of the physical Q response carried by a continuum field.
It is also not the mixed derivative `partial_Q partial_epsilon U`.

The [completed fixed score](p337-q1-closed-trace-transmission-result.md)
makes this distinction numerical: `V_trace=-0.001904836180602413` and
`J_trace=+0.03826094250721058` in (2). Restoring the N25 prefactor gives
`J_trace+25.5 V_trace approximately -0.010312380098151`, even reversing
the sign of this component attribution without changing the physical
Q response. This is just (8) applied to the two published numbers;
it is not an additional source score.

The two fixed calculations in the
[contract](../analysis/p337_q1_trace_continuation_contract.json) retain
these different meanings. They use the saved Q1 root and the existing
integer histograms, with no new configuration or Q scan. The regular
unlabelled endpoint identity stays zero; this is the separate closed-trace
continuation dictated by (1)-(3).

## Scientific consequence

The finite trace does not need to be guessed from an isolated Q4 colour
character. The packing and character proofs fix its two geometric
connection types, Q1 weight and required quadratic R jet. The single
prescribed Q1 score is now complete and its trace-coefficient response
is strictly negative, establishing actual transmission to original U.
The remaining continuum question is the local pair-to-cut intertwiner
and its scaling content, not an absent finite partition/normalization
map. This result does not identify an activated sqrt(N) field: the
closed trace already has a nonzero Q1 value.

A subsequent [explicit local four-port construction](local-four-port-transmission-result.md)
now supplies one concrete microscopic End(pair) interaction and its
strictly positive Q1 original-U response. Its outside-connectivity
support is different from this full seam trace in both directions;
the two-row and contractible-path witnesses rule out identifying the
two configurationwise. Thus a local route is available, while the
claimed continuum field and its relation to the global isotypic trace
remain separate scale-dependent questions.
