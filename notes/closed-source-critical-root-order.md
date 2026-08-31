# Positive closed-source coupling: unique finite root and its order bounds

Pure theory from `bc17b81d`; no new Monte Carlo, enumeration, numerical
response or validation run enters this note. The geometry is an honest
square torus whose full occupied graph has ambient rank2. For the
checkerboard family, periods preserve parity, so A and B each have M=N/2
sites and every NN edge joins opposite colors. Pooled means use fixed
positive weights summing to1, after normalizing each geometry separately.

## 1. Attractive form of the named source

Let beta1 be the occupied graph cycle dimension, r its ambient winding
rank, beta_null=beta1-r, and q=r-1. The closed action has the exact form

```text
J=beta1+beta_null=2 beta1-r,
S_star=C+F+Bv=J-3K+2N+1.
```

The geometric input is that J is supermodular on occupied-site sets and
`0<=Delta_v J<=6`. A short dimension argument explains these properties.
Both the full graph-cycle space and its zero-winding subspace are fixed
linear spaces on the ambient edge set. Restriction to induced occupied
edges gives the corresponding spaces for a vertex subset X. For X,Y,
their intersection is the space for X intersect Y, while their sum lies
inside the space for X union Y. The dimension formula therefore makes
each of beta1 and beta_null supermodular and increasing.

More locally, inserting v with e occupied contact edges in c touched
components creates `Delta beta1=e-c`. If h is its ambient-rank gain,
`0<=h<=e-c`, hence

```text
Delta J=2(e-c)-h.
```

For e=0 this is zero. Otherwise c>=1 and e<=4, giving the upper bound6.
This is a graph-cycle/zero-winding statement, not a claim that every new
cycle is contractible when the previous ambient rank is already positive.

For finite t>=0 the law can be written

```text
mu(omega) proportional to exp[z_A K_A+z_B K_B+t J(omega)],
z_A=logit(p_A)-3t,  z_B=logit(p_B)-3t.
```

Its log weight is supermodular: linear fields are modular and t J is
supermodular. Thus the finite law satisfies FKG positive association.
Conditioning any collection of A sites to fixed occupied/vacant values
preserves this property: restricting a supermodular function to that
coordinate subcube is still supermodular. At p_A=1 the law is precisely
the all-A-occupied conditional law, with strictly positive weights on
every B configuration whenever 0<p_B<1 and t is finite.

## 2. Strict covariance needs a pivotal configuration, not just FKG

Here is a useful finite strictness lemma. Let mu be an associated law
with full support on its free-coordinate cube. Let f be increasing and
have a strictly positive increment at coordinate v for at least one
configuration of the other free coordinates. Put

```text
f0(eta)=f(eta,v=0),  Delta_v f=f(eta,1)-f(eta,0),
pi_v=E[n_v],  f=f0+n_v Delta_v f.
```

The function f0 is increasing and independent of v. Therefore

```text
Cov(f,n_v)=Cov(f0,n_v)+(1-pi_v)E[n_v Delta_v f]>0.
```

The first term is nonnegative by association. The second is strictly
positive because 0<pi_v<1 and the pivotal configuration has positive
weight. This supplies the strict inequality that FKG alone does not.

Apply the lemma to q. With all A occupied, B empty has rank0 while B full
has rank2. Along any ordering that adds all B sites, some B insertion is
pivotal; a direct rank0-to2 jump is allowed. For s<1 this configuration
with all A occupied has positive probability in the full-support law.
For s=1 it belongs to the full-support B conditional law. Hence, for
every geometry and every interior p at finite t>=0,

```text
Cov(q,K_B)>0,  Cov(q,K_A)>=0.
```

The same lemma gives Cov(q,K)>0 for a homogeneous finite law without
requiring a checkerboard: q differs between the empty and full graphs.

## 3. The entire saturation path has one simple pooled root

Now fix `p_A=s+(1-s)p, p_B=p`, with s in[0,1], t finite and nonnegative.
For s<1 ordinary differentiation of the normalized weights gives

```text
partial_p E[q]
  =(1-s) Cov(q,K_A)/[p_A(1-p_A)]
    +Cov(q,K_B)/[p(1-p)].
```

Using `1-p_A=(1-s)(1-p)`, its nonsingular extension to the endpoint is

```text
partial_p E[q]
  =Cov(q,K_A)/[p_A(1-p)] + Cov(q,K_B)/[p(1-p)] > 0.
```

At s=1, K_A=M is constant, so the first term is exactly zero. The strict
B term survives. Pooling with positive geometry weights yields D=Q_p>0.

At p=0 all occupied sites lie in A and no NN edge is occupied, so every
configuration in the limiting support has q=-1. At p=1 the graph is full
and q=1. A finite exp(t S_star) cannot change either support conclusion.
Thus Q is continuous on[0,1], with Q(0)=-1 and Q(1)=1, and strictly
increasing on(0,1). For **every** finite t>=0 and s in[0,1], exactly one
root p0(s,t) exists, it lies in(0,1), and its D is strictly positive.

Finite partition functions and the implicit-function theorem make this
root locally real analytic, with the one-sided interpretation at s=0,1
or t=0. No turning point or zero-D root can occur in the stated domain.
Consequently the established source-preserving endpoint U map no longer
needs a provisional simple-root qualification for finite positive t.
This does not assert that U itself is nonzero.

Also Q is nondecreasing in s: increasing s increases only the A field.
For s<1, `Q_s=pool Cov(q,K_A)/[p_A(1-s)]>=0`; continuity and analyticity
give the endpoint one-sided order as well. Hence `p0,s=-Q_s/D<=0`.
This order of the critical p root is distinct from any claim about U_s.

## 4. Homogeneous chemical-potential root and a finite-t envelope

In the homogeneous law set

```text
z=logit(p)-3t,  mu_(z,t) proportional to exp[zK+tJ],
C_K=pool Cov_(z,t)(q,K)>0,
C_J=pool Cov_(z,t)(q,J).
```

These are averages of **within-geometry** covariances, not a covariance
after mixing geometries. The pooled root Q(z0(t),t)=0 satisfies exactly

```text
z0'(t)=-C_J/C_K.
```

J is increasing. Since every one-site increment of J is at most6,
`6K-J` is also increasing. Applying FKG to each with increasing q gives

```text
0<=C_J<=6 C_K,
-6<=z0'(t)<=0.
```

Thus the natural activity root is nonincreasing with positive coupling.
For the original Bernoulli parameter, however,

```text
ell0(t)=logit(p0(t))=z0(t)+3t,
ell0'(t)=3-C_J/C_K,
-3<=ell0'(t)<=3.
```

The negative density term in S_star is why monotonicity in z does not
fix the sign of the p-root drift. Integrating gives, for 0<=t0<=t1<infinity,

```text
logistic(ell0(t0)-3(t1-t0)) <= p0(t1)
                           <= logistic(ell0(t0)+3(t1-t0)).
```

In particular, the zero-coupling critical root anchors a coefficient-free
finite-t bracket with logit radius3t. Equivalently the natural root lies
between `z0(0)-6t` and `z0(0)`. These bounds do not require a fitted source
amplitude or an asymptotic scaling assumption. The fixed point of the
finite count dictionary is not being promoted to a homogeneous RG flow.

## What the order argument does not decide

q is increasing, but E=q^2 takes values1,0,1 as rank increases from0 to2;
it is not monotone. P4 has signed weights, and U is a ratio of thermal
derivatives evaluated at a moving root. Therefore none of the inequalities
above supplies a sign for U_t, U_s, their mixed gain residual, or an E
response. Those remain separate source/observer questions. The proof is
restricted to finite t>=0 and the declared honest finite graphs; it makes
no negative-coupling, infinite-volume or continuum-field claim.
