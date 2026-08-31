# The closed cycle gas forces finite-volume global-U turnover

**Prediction before the finite-coupling readout.** At fixed finite volume,
the positive closed-source law becomes an empty/full two-state mixture
on its pooled critical root as t tends to infinity. Its rank-one sector
vanishes, and the original angular observable satisfies

```text
U(t)=O(exp(-3t)) -> 0.
```

This is an upper bound, not an identified leading decay exponent. For
the existing N25 pair, the already established positive U(0) and U_t(0)
therefore force at least one maximum at a finite positive coupling.
No value at the newly frozen m=2,4,8,16 points is used in this derivation,
which starts from `b70dc4bd`. There is no new computation or enumeration.

## 1. An elementary four-edge cut bound

The following applies to any connected simple four-regular abelian Cayley
graph `Cay(G,{+a,-a,+b,-b})`, hence to honest square quotients with these
properties. Simplicity and four-regularity mean `ord(a),ord(b)>=3` and
`a!=+b,-b`. The a and b edge families are disjoint unions of cycles.

For a nonempty proper vertex subset X, each directional boundary count
is even: membership changes an even number of times around a cycle.
Connectedness excludes total boundary0. Suppose the total boundary were2.
One direction must then have boundary0, say a. This makes X invariant
under translation by a. The set of starting vertices of crossing b edges,

```text
E_b={x: 1_X(x)!=1_X(x+b)},
```

is also a-translation invariant, by commutativity. Translation by a has
orbits of length ord(a)>=3, so |E_b| is divisible by ord(a) and cannot be2.
Every undirected b edge is counted once in E_b, since ord(b)>2. This is
a contradiction. Thus every nonempty proper X has at least four boundary
edges. No catalog or enumeration of subsets is needed.

In particular the axis N25 torus has directional cycles of length5;
the tilted `(4,3)` Gaussian quotient has both directional orders25.
Both therefore satisfy the bound. The general argument is only used
where the source's original honest edge/face and winding identities also
apply; it does not justify silently deleting quotient incidence aliases.

## 2. A nonnegative defect action separates the two extreme states

Use the fixed source and notation

```text
J=2 beta1-r,  S_star=J-3K+2N+1,
g=2K-J=Bmix-2 C_B+r=2N+1-K-S_star.
```

Here `Bmix` counts occupied-vacant NN boundary edges. The expression for
g follows from `beta1=Bocc-K+C_B` and `4K=2Bocc+Bmix`.
The empty and full configurations each have g=0.

For any nonempty proper occupied configuration, apply the cut bound to
each of its C_B occupied components. There is no occupied edge between
different components, so their boundaries sum exactly to Bmix. Therefore

```text
Bmix>=4 C_B,
g>=2 C_B+r>=2.
```

In particular every rank-one configuration has g>=3; every proper rank-two
configuration has g>=4. Full occupation is the exceptional rank-two state
with g=0. Thus the only zero-defect configurations are empty and full.
The bound need not be attained in a particular rank sector or quotient.

## 3. Exact weights, two-state limit and the pooled root

Let `y=p/(1-p)` be the original homogeneous Bernoulli activity, and set

```text
m=exp(t),  h=y/m=exp(logit(p)-t),  lambda=1/m.
```

Dropping only configuration-independent factors, the exact weight is

```text
y^K exp(t S_star) proportional to h^K m^(-g).
```

For each geometry alpha the reduced partition function and q numerator
are consequently finite polynomials in h and lambda,

```text
Z_alpha=1+h^N+sum_(proper omega) h^K lambda^g,
H_q,alpha=-1+h^N+sum_(proper omega) (r-1)h^K lambda^g.
```

Uniformly for h in a fixed compact neighborhood of1, these expressions
and their h derivatives give

```text
Z_alpha=1+h^N+O(lambda^2),
Q_alpha=(h^N-1)/(1+h^N)+O(lambda^2).
```

Normalize each geometry separately before taking the original pooled
Q. Its limiting root is h=1 and its limiting slope there is N/2. The
analytic implicit-function theorem at `(h,lambda)=(1,0)` gives

```text
h0(lambda)=1+O(lambda^2),
logit(p0(t))=t+O(exp(-2t)),
Q_h(h0(t),t) -> N/2.
```

The [positive-coupling root theorem](closed-source-critical-root-order.md)
identifies this branch with the unique original pooled root for all
sufficiently large t and supplies its continuation through every finite
t>=0. There is no competing critical-root branch to select by a fit.

At fixed h the limiting measure has empty/full probabilities
`1/(1+h^N)` and `h^N/(1+h^N)`. Along the pooled root each geometry tends
to the equal mixture of those two configurations, with total-variation
error O(exp(-2t)). Although the *reference parameter* p0 tends to1,
the actual tilted mean occupied fraction tends to1/2, not1.

## 4. Rank-one suppression makes the original U vanish

Let P1_alpha denote the normalized probability of ambient rank1. Its
numerator contains only terms with g>=3. Thus near the critical branch

```text
P1_alpha=O(lambda^3),  partial_h P1_alpha=O(lambda^3).
```

The same derivative bound follows directly from the finite polynomials
and their nonzero partition denominators; it does not follow merely by
differentiating an unspecified big-O remainder.

For the fixed rank readouts, `E=q^2=1-P1`. The original P4 projector has
constant geometry weights summing to zero. Hence

```text
Y=P4(E)=-P4(P1),  Y_h=O(lambda^3).
```

At fixed t, h is a common invertible thermal coordinate in both
geometries. Its Jacobian cancels in the original slope ratio, giving

```text
U(t)=[N^(13/8)/2] Y_p/Q_p
    =[N^(13/8)/2] Y_h/Q_h
    =O(exp(-3t)) -> 0.
```

This controls the same pooled-root U, not an unnormalized E amplitude.
All constants may depend on the fixed volume and geometry pair. A larger
minimal rank-one defect action or angular cancellation can make the decay
faster; neither the exact exponent nor its coefficient or sign is claimed.

## 5. A positive-coupling maximum is now compulsory for N25

The existing exact N25 result has

```text
U(0)   = +0.8804661569633677,
U_t(0) = +0.12616536341416915.
```

Only these already established signs are needed. Analyticity and the
positive derivative imply U(delta)>U(0)>0 for some delta>0. Since U tends
to0, choose a finite T beyond which every value is below U(delta).
Continuity then puts a global maximum on[0,T], and it cannot occur at0
or T. Consequently at least one finite t_star>0 satisfies

```text
U(t_star)>U(0)>0,  U_t(t_star)=0.
```

There must also be a subsequent interval of decreasing U. This rejects
monotonic increase of the original U throughout positive coupling, not
merely a linear-response extrapolation. It does not locate the maximum,
prove it unique, or require it to fall between any two of the four frozen
coupling points. Those values remain an independent, already specified
finite-coupling readout; the prediction here is not adjusted to them.

The mechanism is finite-volume concentration onto two extremal rank
sectors: strong attraction removes the intermediate rank-one probability
that carries Y. The order of limits is **fixed N, then t to infinity**.
This is neither a thermodynamic phase transition nor a universal exponent,
and it does not reopen the stopped F4-only experiment.
