# Vacuum KdV I3 at c=0: degeneracy and fixed-point sector selection

The existing #231 oracle evaluates `K4=D2D0` on Pinson--Arguin sector numerators.  This note extracts two exact consequences that do not require fitting its Gaussian sums.

## 1. The unconditioned c=0 vacuum charge degenerates

In the normalized `Q=1` random-cluster theory,

```text
Z_total=1.
```

At `c=0`, the first KdV charge is

```text
<I3>=D2D0 Z=(D-E2/6)DZ,
```

so

```text
I3[Z_total]=0.                                      (1)
```

The same is true for the exact second-charge operator quoted in Issue #231,

```text
<I3^2>=[D^4+(1/18)E4D^2-(11/1080)E6D]Z=0.          (2)
```

Equations (1)--(2) are not an absence of the vacuum KdV current.  They say that, at `c=0`, its normalized total-torus expectation has no scalar component.  The nonzero response is a redistribution among restricted sectors.

This is a selection rule for the matching-even interpretation:

> a universal nonzero even-H4 amplitude can be called vacuum `I3` only after a topological/defect/sector instrument is declared.  An unconditioned scalar mean cannot be the normalized `c=0` vacuum charge without an additional contact or nonuniversal term.

## 2. Vector-valued modular covariance fixes the rho direction

Let `v(tau)` be the chiral `K4` response vector of the three registered primitive sectors

```text
[(1,0),(0,1),(1,1)].
```

It is a weight-4 vector-valued modular form.  The hexagonal point

```text
rho=exp(i pi/3)
```

is fixed by

```text
gamma(tau)=(tau-1)/tau,
c tau+d=tau.
```

At a fixed point, weight `w` covariance gives

```text
R_gamma v(rho)=rho^(-w)v(rho).                     (3)
```

For `w=4`, `rho^-4=omega=exp(2pi i/3)`.  In the registered forward sector cycle, the allowed vector is therefore

```text
v_KdV(rho)=A(1,omega,omega^2).                     (4)
```

The included Gaussian-series oracle verifies (4) to more than 68 decimal digits, but (4) follows from the modular fixed-point equation, not from the numerical value of `A`.

The exact C3 projectors make the purity explicit:

```text
P_omega v=v,
P_1 v=P_omega2 v=0.
```

## 3. Reflection-even response and exact nulls

A real aligned coupling adds the conjugate chiral response.  From (4),

```text
v+vbar=A(2,-1,-1).
```

Hence the amplitude-free primitive-sector prediction is

```text
R_10:R_01:R_11=2:-1:-1,                           (5)
Q_rho=0,
S_rho=0,
C_rho=3A != 0 is allowed.
```

The two zero-parameter complex residuals are

```text
K4_01-omega K4_10=0,
K4_11-omega^2 K4_10=0.
```

They should be scored before estimating any lattice amplitude.

## 4. Why this differs from the two live spin-4 competitors

For a scalar/trivial-sector thermal-Q4 one-point, (3) would demand

```text
1=rho^-4=omega,
```

which is impossible.  Its value is forced to zero; this is the ordinary `E4(rho)=0` node.  The vacuum KdV sector vector evades the zero only because it occupies the nontrivial `omega` character.

The `x=17/4` four-leg primary is in Potts representation `[2]`.  An unmarked singlet one-point does not see it linearly.  Thus a nonzero response with the exact sector ratios (4)--(5) is neither a hidden scalar thermal Q4 amplitude nor an unprojected four-leg amplitude.

There is also a built-in weight-8 adversary.  At rho, weight 8 requires the conjugate character

```text
rho^-8=omega^2.
```

The reflection-even `2:-1:-1` projection alone loses that chirality and cannot distinguish weights 4 and 8.  The complex sector response does: the order of `omega,omega^2` reverses.  This is the minimum reason to retain the chiral primitive-sector numerator rather than only `C`.

## 5. Frozen score

At one rho geometry, retain the complex response in the registered sector order and score

```text
K_01/K_10=omega,
K_11/K_10=omega^2.
```

The real fallback score is (5), but it does not reject the weight-8 conjugate.  At the square point `tau=i`, the weight-4 stabilizer instead requires the `+1` eigenspace of the sector swap, giving the independent exact equality `K_10=K_01`.

## Boundary

The modular and `c=0` statements are exact.  Identifying a matching-even lattice row with this response remains conditional on using the same sector instrument, frame, and real/chiral coupling convention.  The result intentionally forbids calling an unconditioned scalar H4 mean `I3` by exponent alone.
