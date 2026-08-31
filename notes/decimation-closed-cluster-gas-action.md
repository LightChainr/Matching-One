# A decimation-closed square-site cluster-gas action

This is a named finite microscopic model, not an inference from the new F4
production. No F4 values, residuals, fitted coefficients, simulations or
enumerations enter this note. The geometric saturation dictionary extends
the component/fugacity map in
[`square-checkerboard-endpoint-homology.md`](square-checkerboard-endpoint-homology.md).
Here we derive its smallest closed action and its thermal-coordinate gauge.

## 1. Three statistics and a three-dimensional exact closure

On an honest square-cell torus of N sites, write

```
C = C_B,NN + C_W,matching,
F = number of fully occupied unit faces,
Bv = number of vacant-vacant NN edge orbits.
```

Occupied/vacant cluster counts exclude empty components. The NN edge set
has2N edges; faces and their incidence convention stay fixed. Saturate one
checkerboard color, divide periods by1+i, and complement to the ordinary
square child. The exact pullback dictionary, in child statistics, is

```
T C=C+F,       T F=Bv,       T Bv=0.
```

The C correction counts isolated saturated sites. A filled parent face
has its two nonsaturated corners occupied, hence corresponds after child
complement to one vacant-vacant child NN edge. There are2M parent faces and
2M child NN edges for child area M. No parent vacant-vacant NN edge survives
when every vertex of one color is occupied.

For S=aC+bF+cBv the coefficient update is

```
(a,b,c) -> (a,a,b).
```

Fix a=1. The fixed equations force b=1,c=1. Therefore the unique fixed
positive combination in this span is

```
S*=C+F+Bv,        T S*=S*.
```

In fact T²(aC+bF+cBv)=aS*. In particular the bare cluster count has orbit
`C -> C+F -> C+F+Bv -> C+F+Bv`. The span is minimal: differences between
these first three members recover F and Bv. These functions are genuinely
linearly independent; all-vacant, all-occupied and a singleton have triples
`(1,0,2N)`, `(1,N,0)`, `(2,0,2N-4)`, with nonzero determinant.

The two-step orbit is an algebraic identity. To realize **two successive
finite checkerboard saturations**, periods must be divisible by `(1+i)²`,
and the child incidence convention must remain valid. A parent with odd/odd
Gaussian periods and N=2 mod4 has an odd-area first child and does not admit
a second checkerboard operation. One legal step still fixes S* exactly.

## 2. Explicit cluster, bond and cycle-gas forms

Let K be the occupied-site count, Bocc the occupied-occupied NN edge count,
and q=r-1 the ambient-rank observer. Incidence counting gives

```
Bocc+Bmix+Bv=2N,      2Bocc+Bmix=4K,
Bv=2N-4K+Bocc.
```

The finite torus Euler identity is

```
C_B-C_W = K-Bocc+F+q.
```

Consequently all the following forms are identical configuration by
configuration, including the topological correction:

```
S* = C+F+Bocc-4K+2N
   = 2C_B+2Bocc-5K-q+2N
   = 2 beta1-3K-q+2N,
beta1=Bocc-K+C_B.
```

Beta1 is the occupied **graph** cycle-space dimension, not ambient winding
rank. If beta_null=beta1-r is the dimension of its zero-ambient-image
subspace, another useful form is

```
S*=2 beta_null+r-3K+2N+1.
```

At fixed K the action rewards graph cycles, with the displayed ambient-rank
correction. It counts more than filled elementary faces: a loop surrounding
a vacant patch can contribute even when no face is fully occupied. Dropping
the term -q would change the relative weights of the rank sectors and is
not an innocuous normalization of the global q/E problem.

## 3. A finite positive-weight statistical-mechanical family

For 0<p<1 and real t define

```
mu*_(p,t)(omega) = Z*_(p,t)^(-1)
                  p^K(1-p)^(N-K) exp[t S*(omega)].
```

This is a well-defined finite cluster-gas measure with strictly positive
weights. With y=p/(1-p), and dropping only configuration-independent
factors, its Boltzmann weight can be written either as

```
(y exp[-4t])^K exp[t C+t F+t Bocc]
```

or as

```
(y exp[-5t])^K (exp[2t])^C_B (exp[2t])^Bocc exp[-t q]
 = (y exp[-3t])^K (exp[2t])^beta1 exp[-t q].
```

Thus the proposed action has explicit site activity, cluster fugacity,
occupied-bond weight and topological-sector weight. No fitted coefficient
or assumption about a continuum operator supplies this dictionary.

For inhomogeneous independent base occupations p_A=1,p_B=p, multiply by
the same exp(tS*) and normalize. The endpoint configuration bijection and
T S*=S* give the **full tilted measure** identity

```
pushforward mu*_(parent;1,p,t) = mu*_(child;1-p,t),
Z*_(parent;1,p,t) = Z*_(child;1-p,t).
```

The second equality uses normalized Bernoulli base weights. For unnormalized
occupation-odds polynomials there is instead the familiar factor y^(N/2).
This extends the bare Bernoulli endpoint to a specified interacting family:
the same t is retained under the hard-saturation map.

## 4. Removing -4K changes only the common thermal chart

Let S_local=C+F+Bocc. Since S*=S_local-4K+2N, their normalized finite tilts
are related by

```
mu*_(p,t) = mu_local_(p_eff,t),
logit(p_eff)=logit(p)-4t.
```

The constant2Nt drops out of normalization. At fixed t this is one common
smooth increasing p-coordinate change for every orientation and both q,E.
After finding the corresponding pooled roots, its Jacobian cancels in
`U=N^(13/8)/2 * Y_p/Q_p`. Therefore the two finite families have **identical
root/slope-normalized U(t)** wherever that root and nonzero slope exist.
Their p roots differ; their literal microscopic statistics are not equal.
S_local is a thermal-gauge equivalent action, not a second fixed coefficient
vector in the three-statistic dictionary above.

For genuine homogeneous equilibrium tilts at t=0, linearity also gives
`V_S*=V_C+V_F+V_Bv=V_C+V_F+V_Bocc`, since the K tilt is a common logit shift.
This identity requires the same source convention, root and covariance.
It does not license substituting a lag1 rank-centered C response or treating
correlated terms as independent observations.

## 5. Nonthermal measure, unmeasured global transmission

Take an untwisted torus with enough room to embed the following four-site
patterns without wrapping, e.g. side length at least8. Their vacant matching
graph stays connected and q=-1:

| Occupied pattern (K=4) | C | F | Bocc | S* |
|---|---:|---:|---:|---:|
| four mutually NN-separated sites |5|0|0|2N-11|
| four-site straight path |2|0|3|2N-11|
| filled2x2 square |2|1|4|2N-9|

C, F and Bv therefore are not individually density/constant aliases, and
S* itself is not one: the square/path relative probability is multiplied
by exp(2t) at the same K. No change of homogeneous Bernoulli p can do that.
The cycle-gas form explains why distinct forests at the same K have equal
S* while adding a contractible cycle changes it by2.

This proves a nonthermal change of the **measure**, not an observed change
of the particular projected global U. A fixed-face experiment alone has
not measured the complete S* tilt or its finite-t response. No new F4
number is used here, no failed source is rescued, and no prescribed
transmission result is inferred from coefficient closure. The closure is
under a deterministic, fully saturated endpoint map; it is not an ordinary
integration-out RG fixed point, a continuum field identity, an asymptotic
amplitude, or a universality theorem.
