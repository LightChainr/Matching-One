# Exact collapse of the bivariate rank source

Status: terminal-state consequence of the typed incidence spine and the
Issue #269 digital-Alexander theorem.

The incidence-spine state sum is

```text
Phi(p;x,y)=E_p[x^r_b y^r_w].
```

Issue #269 proves configurationwise

```text
r_b+r_w=2,                 q=(r_b-r_w)/2,
r_b=1+q,                   r_w=1-q,
q in {-1,0,+1}.
```

Therefore every monomial factors before any averaging:

```text
x^r_b y^r_w = x y (x/y)^q,
```

and hence

```text
Phi(p;x,y)=x y Z_rel(p,x/y),                            (1)

Z_rel(p,Q)=E_p[Q^q]
          =P_-(p)Q^-1+P_0(p)+P_+(p)Q.                 (2)
```

Equivalently, the entire terminal support is

```text
Phi=P_- y^2+P_0 x y+P_+ x^2.                           (3)
```

The local transfer calculation can still require a large typed frontier
partition.  Equations (1)--(3) say that its **terminal rank-source output**
has only three coefficients.

## Diagonal and relative logarithmic sources

Put `x=exp(a)` and `y=exp(b)`.  Then

```text
Phi(exp(a),exp(b))=exp(a+b) G(a-b),
G(s)=E[exp(s q)].                                      (4)
```

The diagonal source is deterministic because `r_b+r_w=2`:

```text
(x d_x+y d_y) Phi=2 Phi,
(x d_x+y d_y)^m Phi=2^m Phi.                           (5)
```

All stochastic information lies in the relative Cartan source

```text
D_rel=(x d_x-y d_y)/2,
D_rel Phi at x=y=1 = E[q]=M(p).                        (6)
```

There is no independent total-rank susceptibility or higher diagonal-rank
response.

## Strict rank-one covariance

Since `(r_b,r_w)=(1,1)+q(1,-1)`,

```text
Cov(r_b,r_w)
 =Var(q) [[ 1,-1],
          [-1, 1]].                                    (7)
```

For an honest nonempty torus and `0<p<1`, both the empty mask (`q=-1`) and
the full mask (`q=+1`) have positive Bernoulli weight.  Thus `Var(q)>0` and
the matrix in (7) has rank exactly one.  At `p=0` or `p=1`, `q` is
deterministic and the covariance rank drops to zero.

More generally, every connected rank cumulant of order `m>=2` is

```text
kappa(r_i1,...,r_im)
 =kappa_m(q) v_i1 ... v_im,       v=(1,-1),            (8)
```

so every rank-cumulant tensor has the same one-dimensional relative-source
direction.

## Three-state closure inherited from Issue #54

The support polynomial `q^3=q` gives

```text
G_sss=G_s,
G=P_0+P_+ exp(s)+P_- exp(-s).                          (9)
```

Every positive odd raw source derivative equals `G_s`, and every positive
even derivative equals `G_ss`.  For `F=log G`,

```text
F_sss=F_s-3 F_s F_ss-F_s^3.                            (10)
```

At `s=0`, `mu=F_s=E[q]` and `v=F_ss=Var(q)` reconstruct the terminal vector:

```text
P_+=(v+mu^2+mu)/2,
P_-=(v+mu^2-mu)/2,
P_0=1-v-mu^2.                                          (11)
```

Thus higher terminal rank cumulants contain no sector information beyond
the first two relative-source derivatives.  This is the rank-source closure
from Issue #54; it does not close the typed connectivity frontier or thermal
derivative direction.

## Axis `L=3` exact certificate

`scripts/p144_relative_source_collapse.py` consumes the 512-mask oracle from
`p144_typed_incidence_spine.py`.  It checks that the only rank pairs are

```text
(0,2), (1,1), (2,0),
```

and verifies (1)--(11) exactly at `p=2/5`.  The three sector probabilities
are retained as rational numbers, the covariance determinant is exactly
zero while `Var(q)>0`, and the matching expectation from the relative source
equals the derivative specialization of the original bivariate state sum.
