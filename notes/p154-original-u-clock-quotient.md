# Original U: clock quotient, three contrast jets, and the finite-step remainder

This algebra companion supplies the exact interface for the P154 canonical
birth-source calculation. The common-clock zero prediction was already in
Issue154; it is **not a new mechanism discovery**. No profile, source alpha,
Monte Carlo block, covariance or power estimate is computed here. It does
not alter the prospective no-launch decision or supply missing statistical
resolution. The existing P267 full-profile finite-warp/J-moment theory is
also not rederived as a new result.

## The precise projection of the original U

At fixed N, let Q_t=mean(q_t), Y_t=P4(E_t), A=N^(13/8)/2, and let p_* be
the pooled root Q_0(p_*)=0. Means and P4 have fixed geometry weights. Assume
D=Q_0'>0 locally and sufficient smoothness. With x=Q_t(p), define
R_t(x)=Y_t(Q_t^{-1}(x)). Then the **original** observable is exactly

```
U_t = A R_t'(0).
```

Primes on Q,Y and source profiles below mean p derivatives; derivatives
of R mean x derivatives. For j_Q=mean(J_q), j_Y=P4(J_E), write

```
r(p) = Y'(p)/Q'(p),
Zeta(p) = j_Y(p) - r(p) j_Q(p),
partial_t R_t(x)|0 = Zeta(Q_0^{-1}(x)),
v = partial_t U_t|0 = A/D * Zeta'(p_*).
```

Expanding the last equation gives the existing moving-root expression,
including the derivative of the ratio r, not a frozen-r approximation:

```
v = A/D [j_Y' - r j_Q' - r' j_Q],
r' = (Y'' D - Y' Q'')/D^2.
```

A shared scalar coordinate change q_g,t=q_g,0 composed with phi_t and
E_g,t=E_g,0 composed with the **same** phi_t makes R_t independent of t.
Infinitesimally j_Q=wD,j_Y=wY', so Zeta vanishes pointwise for arbitrary
p-dependent w. This is scalar transport of the probabilities, not a
density pullback with an additional Jacobian. Conversely, Zeta=0 locally
identifies a common coordinate tangent for the **projected pair (Q,Y)**;
it does not reconstruct all individual geometry profiles.

The strict information hierarchy is

```
common warp of every geometry/birth
    => Zeta identically0 => Zeta'(p_*)=0 <=> v=0.
```

A resolved fresh nonzero v therefore rejects the common-coordinate source
tangent, as already recognized in Issue154. A zero v does not establish
that tangent, let alone a finite shared warp. A derivative at t=0 also
does not decide whether a nonlinear source path leaves and later returns
to the same finite-warp orbit.

## Exact four-alpha quotient: three contrast jets, not three theories

Let tau_first=-1,tau_completion=+1 and eta_g=+1,-1 for the two geometries.
Where F'_{j,g} is nonzero, define the physical effective clock
alpha_{j,g}=N J_{j,g}/F'_{j,g}. Its pointwise Hadamard decomposition is

```
alpha_jg = a00 + eta_g a10 + tau_j a01 + eta_g tau_j a11,
J_q,g = [(a00+eta_g a10)q_g' + (a01+eta_g a11)E_g']/N,
J_E,g = [(a00+eta_g a10)E_g' + (a01+eta_g a11)q_g']/N.
```

This identity verifies the signs for first/completion. Where a denominator
vanishes, retain the source J rather than inventing a finite alpha ratio.
The whole a00(p) term, including a00', cancels from Zeta and v. Define
the following three **baseline-only** functions with the same fixed weights:

```
H10 = P4(eta E') - r mean(eta q'),
H01 = P4(q')     - r mean(E'),
H11 = P4(eta q') - r mean(eta E').
```

Then the exact original-U transmission is

```
v = A/(N D) sum_{m=10,01,11} [a_m' H_m + a_m H_m'] at p_*.
```

Thus the gains multiplying the contrast values and derivatives are
A H_m'/(N D) and A H_m/(N D), respectively. Existing baseline second
p-jets suffice to calculate them without fitting or rereading source
alphas. These are three simultaneous contrast coordinates, not mutually
exclusive models. Equality of all four alphas **at the root alone** is
insufficient: their three contrast derivatives can still transmit to v.
One scalar v cannot identify the six contrast-jet values separately.

## Constant microscopic event mark: exact finite-step fingerprint

For any one birth CDF F(p), let b(k)=Pr(K=k), k=1..N, and
S(p)=sum_k b(k) BinomialPMF(N,k;p). A constant event mark a gives J=aS.
The elementary Bernstein identity is

```
N B_{k,N} + (1-p) B'_{k,N} = N B_{k-1,N-1},
[N+(1-p)partial_p]J = a F',  J(0)=0.
```

Consequently constant microscopic a does **not automatically** mean an
exact common canonical clock: the effective alpha=N J/F' can depend on
the birth and geometry. No universal nonzero U response follows either.

Let L denote the linear original-U source functional A/D*Zeta', applied
to the whole birth-profile vector, and put c=1-p. The following are exact:

```
v = -(1/N) L(c J'),
J = (a/N)F' - a/[N(N-1)] c F'' + c^2/[N(N-1)] J'',
v = -a/[N(N-1)] L(c F'') + 1/[N(N-1)] L(c^2 J'').
```

The last two lines use N>1 and common a for every birth/geometry. They
expose the finite-step correction without silently discarding its
remainder. In a controlled local derivative expansion this yields the
proposed leading `-a/N^2 L((1-p)F'')`, with the stated minus sign.

For a general smooth scalar c(p), the quotient gives, exactly,

```
L(c F'') = A [(2c Q'' + c' D) R'' + c D^2 R'''] at x=0.
```

Indeed j_Q=cQ'', j_Y=cY'', so Zeta=c D² R'' and differentiation gives
this expression. In particular the leading constant-mark fingerprint is

```
-a A/N^2 * { [2(1-p_*)Q'' - D]R'' + (1-p_*)D^2 R''' }.
```

It probes curvature and its derivative in the physically fixed Q clock,
not the rigid-translation mode. Its sign is not universal and the two
terms can cancel. The expansion is not a finite-N identity: for example
F=p has J=a[1-(1-p)^N]/N, which retains a boundary term although F''=0.
The exact identities with J'' above remain valid without any asymptotic
error assumption. No power improvement is inferred from this algebra.

## Explicit nonclock responses invisible to U

These counterexamples respect valid ordered-birth CDFs rather than merely
postulating arbitrary perturbations of a graph. Put x=2p-1 and

```
q_+(p)=q_-(p)=x,
E_±,t(p)=(1+x^2)/2 ± t h(x),
F1_±=(x+2-E_±)/2,  F2_±=(x+E_±)/2.
```

For either h0=(1-x²)² or h2=x²(1-x²)² and |t|<1/8, these are increasing
CDFs with endpoints0,1 and F1>=F2. One sufficient derivative bound is
|h'|<=8(1-|x|); the endpoint values and E<=1 give the remaining conditions.
Use P4(E)=(E_+-E_-)/2 for this example; any fixed nonzero contrast
normalization preserves its zero/nonzero conclusions.

Here R_t(x)=t h(x). With h0 the intrinsic-root value changes by t while
U_t=A t h0'(0)=0. With h2 both the value and slope at the root stay fixed,
but R_t''(0)=2t: a genuine transverse curvature change is invisible to U
for **every** t. Since the strictly monotone q profiles do not change,
their only possible common coordinate map is the identity, which cannot
produce either E deformation. These examples establish the strictness
of the last implication in the hierarchy without proposing new production
observables or treating old full-profile theory as a new discovery.

For the prospective decision this means: a nonzero v kills common-clock
explanation of that source tangent; a tightly small v only limits this
root-local angular slope transmission. It does not establish full clock
closure or exclude transverse value/curvature changes. An unresolved v
at inadequate power establishes neither statement.
