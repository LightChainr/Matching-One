# A rank-one trace enters original U only through a specific transmission functional

**Result.** A nonzero torus trace in a two-cluster colour channel is not
by itself a nonzero coupling to original U. For a contribution supported
in ambient rank1, its raw q and E numerators are exactly zero. It can
nevertheless act through the geometric normalizers. After retaining the
moving pooled root and thermal slope, that entire normalization route
reduces to three specified coefficients in (8) below. A common, thermally
constant fractional partition contribution cancels exactly.

This gives a concrete interface for a finite trace or a *completed*
confluent insertion; it introduces no new observable, fitted source or
candidate power. It starts from `2690f665` and the fixed Q occupation
completion. The current overview's regular endpoint exclusion is retained:
`ell P_[2](Q)=0` throughout that rational family. The representation label
`[2]` is not the ambient topology rank2 used in the occupation law.

## 1. The exact perturbation being considered

For geometry g=f,s, use the unnormalized weights of the already named
closed-source family in a common thermal coordinate z=logit(p). Define

```text
Z_g=sum_A w_g(A),
B_g=sum_A q(A) w_g(A),   C_g=sum_A E(A) w_g(A),
q_g=B_g/Z_g,             E_g=C_g/Z_g,
q=r-1,                  E=q².
```

Suppose a specified finite insertion or derivative packet, with formal
coefficient epsilon, changes only the rank-one part of this sum:

```text
delta Z_g=F_g(z),        delta B_g=delta C_g=0.                 (1)
```

The zeros follow because q=E=0 on rank1, not because the trace itself
vanishes. F_g may be signed if it is one term in a finite representation
decomposition; the baseline Z_g is positive. Equation (1) can equally
describe an actual source derivative supported on rank1. Any singular
colour channels must be combined into a finite F_g **before** (1) is
applied. No infinite individual projector derivative is an input here.

Put f_g=F_g/Z_g. Separate normalization gives the exact first variations

```text
delta q_g=-f_g q_g,           delta E_g=-f_g E_g.              (2)
```

This is not a claim that the full colour-[2] transfer sector has support
only in topology rank1. It is the transmission rule for the rank1-supported
part, including the explicit two-essential-cluster closure. Other supports
have their own raw numerator variations and cannot be silently included
in (1).

## 2. Keep the original root, both normalizers and the thermal denominator

Let the fixed angular denominator be Delta!=0 and write

```text
M=(q_f+q_s)/2,       Y=(E_f-E_s)/Delta,
M(z0)=0,            D=M_z(z0)>0,
R=Y_z/D,            U=A_N R,       A_N=N^(13/8)/2.             (3)
```

All quantities in the following formulas, including the derivatives of
F_g/Z_g, are evaluated at the same z0. The perturbed root displacement is

```text
z_epsilon = mean_g(f_g q_g)/D.                                (4)
```

The complete derivative of the ratio at that moving root is

```text
delta U/A_N = {
 -P4[(f E)_z]+R mean[(f q)_z]
 +[mean(f q)/D] (Y_zz-R M_zz)
}/D,                    P4(v)=(v_f-v_s)/Delta.                (5)
```

To verify (5) directly, vary Y_z/D along the root, use
`delta Y=-P4(fE)` and `delta M=-mean(fq)`, and retain
`z_epsilon(Y_zz-R M_zz)`. This derivation uses only (1)-(3); a fixed-root
ratio or an unnormalized angular numerator would omit genuine terms.
The chart z is the same in both geometries. A density drift may be removed
only in the complete physical derivative, not separately in selected
trace components.

## 3. The common and geometric parts have different roles

Define

```text
f_c=(f_f+f_s)/2,      f_d=(f_f-f_s)/2,
a=(q_f-q_s)/2,        e=(E_f+E_s)/2,
H=Y_zz-R M_zz.                                                (6)
```

At the pooled root the two geometric means are q_f=a, q_s=-a;
neither is assumed to vanish. Before root substitution,

```text
mean(fq)=f_c M+f_d a,
P4(fE)=f_c Y+2f_d e/Delta.                                   (7)
```

Substitute (7) into (5). The two terms proportional to f_c itself cancel
because Y_z=R D. The resulting *entire* rank-one interface is

```text
delta U/A_N = C_c (f_c)_z + C_dz (f_d)_z + C_d f_d,

C_c  = -Y/D,
C_dz = [R a-2e/Delta]/D,
C_d  = [R a_z-2e_z/Delta+(a/D)H]/D.                           (8)
```

These three coefficients are fixed by the existing unperturbed q/E
curves and their thermal derivatives. They are not regression features
or three additional fitted mechanisms. For an actual proposed trace,
F_f and F_s determine the three arguments in (8), and hence one number.

Several consequences follow without a scaling hypothesis:

- If the same fractional contribution f is independent of z in both
  geometries, (8) is exactly zero even when F_f and F_s are nonzero.
- If f_f=f_s=f(z), the complete contribution is `-A_N Y f_z/D`.
  A thermally varying common normalizer can enter, but only by multiplying
  the already present angular value Y; its existence does not by itself
  create a new independent angular sector.
- A geometric fractional difference can enter through both its thermal
  derivative and its value. The latter includes the root displacement
  `a f_d/D`; dropping it is unjustified when individual q_g are nonzero.
- If additionally Y=0, a=0, a_z=0 and e_z=0 at the root, (8) reduces to
  `delta U/A_N=-2e (f_d)_z/(Delta D)`. These are stated conditions, not
  assumptions about the current N25 pair.

Thus “the colour channel has a nonzero torus character” is weaker than
“it has a nonzero original-U coupling.” The extra requirement is the
nonvanishing of the specific right side of (8).

## 4. What a trace or confluent four-leg proposal must actually deliver

For the rank1 part of such a proposal the sufficient finite input is
`F_g(z0), (F_g)_z(z0)` in the same two geometries, together with the
already specified baseline jets. Those values must come from the actual
transfer closure, not from an arbitrary multiple of `tr P_[2]`.

If the proposal concerns a Q derivative, F_g is the **complete finite
log-Q derivative packet** assigned to that part, after including its
measure, projector and explicit operator terms. In particular,
`ell' P_[2]+ell P_[2]'=0` for the regular endpoint remains zero; selecting
only the second term cannot manufacture an F_g. Conversely, a genuine
torus closure need not be that linear endpoint and must be evaluated
with its own prescribed functional.

A claimed activated power, logarithm or modulus relation then needs both
the colour/spectral identification and the predicted nonzero value of
(8), with its size dependence. Neither a projector pole nor a finite
trace character supplies the latter. This is an interface condition on
a named operator, not an invitation to fit additional components.

## Scientific card

- **Mechanism changed:** the rank1 normalization route from a torus trace
  to original U is explicit; a whole common constant direction is exactly
  invisible, and the surviving terms are identified.
- **Observer/sector/source:** original separately normalized pooled q/E
  and U; topology-rank1 part of a specified finite trace/source packet.
  Colour-[2] and topology rank2 are kept distinct.
- **Dependency:** finite algebra from the same action and root definition;
  no new samples, histogram scoring or independent evidence vote.
- **Not shown:** a nonzero full lattice colour-[2] amplitude, its CFT
  identity, or a nonzero value of (8) for that unidentified full sector.
- **Next discriminant:** an actual transfer insertion fixes F_f,F_s;
  evaluate (8) once. A zero value stops that proposed transmission route
  even if its unnormalized torus trace is nonzero.
