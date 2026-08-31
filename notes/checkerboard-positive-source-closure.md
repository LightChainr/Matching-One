# A positive source closed by checkerboard endpoint decimation

**Exact result.** For the ordinary square/matching pair, define

```text
C = occupied NN components + vacant matching components,
F = number of fully occupied unit faces,
Bv = number of NN bonds with both endpoints vacant,
S_star = C+F+Bv.
```

At the fully occupied-A checkerboard endpoint, scaling the B sublattice by
`1/(1+i)` and complementing it gives an ordinary square child with

```text
T C=C'+F',     T F=Bv',     T Bv=0,     T S_star=S_star'.
```

Consequently the same finite positive weight exp(t S_star) transports
exactly, at the same t. This is a different closed source from F alone;
it neither reruns nor changes the stopping decision of the completed F4
experiment. This note starts at d507a873 and uses no revealed numerical
response, fitting, sampling or new enumeration.

## 1. Geometry, winding and the three counting identities

Let Λ be a square-torus period lattice contained in `(1+i)Z[i]`. Occupy
every even site A, retain O⊂B on the odd sites, and set

```text
Λ'=Λ/(1+i),   M=N/2,
W={(w−1)/(1+i):w in O},   U=(Z[i]/Λ')\W.
```

U is the occupied configuration on the ordinary NN child. Period vectors
map by `(x,y)→((x+y)/2,(y−x)/2)`. Work in the nonalias regime with each child
3×3 stencil injective; for Gaussian tori M>8 suffices. All bonds and faces
are their geometric orbit incidences, not silently deduplicated aliases.

The [ambient-homology endpoint proof](square-checkerboard-endpoint-homology.md)
replaces every occupied parent path `w→a→w'` by a child matching edge between
the images of w,w'. Its inverse expansion has the same lifted endpoint
displacement. Thus the entire ambient integer winding image, not merely
ordinary connectivity, maps from parent NN(A∪O) to child matching(W).
The vacant parent matching graph maps directly to child NN(U). Digital
Alexander duality therefore gives, configurationwise,

```text
rank_parent=2−rank_child(U),   q_parent=−q_child(U),   E_parent=E_child(U).
```

Every nonisolated occupied parent component contains an O site and
corresponds to a child matching(W) component. An isolated A site has all
four neighboring B sites absent; these four sites form one full U face.
The A-site/child-face correspondence is bijective. Hence

```text
C_parent=C_matching(W)+C_NN(U)+F_child(U)=C'+F'.
```

For F, each parent unit face contains two already occupied A corners and
two opposite B corners. It is full exactly when both B corners lie in O,
equivalently when their two child sites are vacant. Their images are one
child NN edge. Explicitly, a parent face with lower-left corner z even
maps its odd corners to `z/(1+i)` and `z/(1+i)+i`; a face with z odd maps
them to `(z−1)/(1+i)` and `(z−1)/(1+i)+1`. These are respectively all
vertical and all horizontal child edges, once each. Parent2M faces and
child2M bonds are therefore in bijection, proving F_parent=Bv_child.
Finally every parent NN bond touches a filled A site, so Bv_parent=0.
No loop-count or Betti-number equality is assumed in these component and
incidence identities.

## 2. Repetition and its parity requirement

On coefficients, one endpoint map sends
`a C+b F+c Bv` to `a C'+a F'+b Bv'`. Thus the positive sum S_star is fixed,
and along an admissible sequence of endpoint maps

```text
T^2 C=C''+F''+Bv''=S_star'',    T^3 C=S_star'''.
```

This repetition requires actual checkerboard parity at **every** level.
j steps require `Λ⊂(1+i)^j Z[i]`. For a Gaussian generator alpha this is
divisibility by `(1+i)^j`; the final child area N/2^j>8 is a sufficient
nonalias bound for the whole chain. A one-step legal torus need not admit
a second step: N130→65 and N170→85 stop after one; N260→130→65 and
N340→170→85 admit two. The identities T²C and T³C are not permission to
invent a checkerboard on an odd-area child. Nor are they a claim that all
unconstrained configurations undergo these repeated endpoint restrictions.

## 3. Exact finite fugacity and normalized observer law

For p_A=1,p_B=p with 0<p<1 and fixed finite real t, define the endpoint partition
functional using the Bernoulli reference weight on the M random B sites:

```text
Z_parent[O](p,t)
 = sum_(O_B subset B) p^|O_B|(1−p)^(M−|O_B|)
     exp[t S_star,parent(A union O_B)] O_parent(A union O_B).
```

The configuration bijection sends this Bernoulli factor exactly to the
child factor at r=1−p. The proven equality of S_star gives

```text
Z_parent[O](p,t)=Z_child[O_tilde](1−p,t),
Z_parent[1](p,t)=Z_child[1](1−p,t)>0.
```

Thus the normalized observer law has the same identity. There is no fitted
coupling, Jacobian or additive source correction. All weights are positive
for real finite t; O_tilde is the configurationwise mapped observable.
In particular, for each direct mapped geometry,

```text
Q_parent^end(p,t)=−Q_child(1−p,t),
E_parent^end(p,t)= E_child(1−p,t).
```

Here Q is the matching mean. For any corresponding simple pooled-root
branch, `p0_parent^end(t)=1−p0_child(t)`. This does not assert uniqueness or
nonzero slope for every t. The partition identity is exact irrespective
of those root qualifications.

## 4. Fixed-t thermal derivative and original U

Keep the same ordered pair of mapped geometries. The child angle is
theta'=theta−pi/4, so `delta_cos4_child=−delta_cos4_parent`. With
`Y=P4(E)` this yields

```text
Q_parent^end(p,t)=−Q_child(1−p,t),
Y_parent^end(p,t)=−Y_child(1−p,t),
Q_parent,p^end=Q_child,p,    Y_parent,p^end=Y_child,p.
```

The first minus in Y is the angular-projector sign; differentiating the
complemented probability supplies the second. Derivatives hold t fixed.
They vary only the random B probability at the parent endpoint—as in
`p_A=s+(1−s)p,p_B=p` at s=1—not both sublattice probabilities.

At the corresponding roots, the frozen original normalization therefore is

```text
U_N^end(t) = [N^(13/8)/2] Y_parent,p^end/Q_parent,p^end
          = 2^(13/8) U_M(t),      M=N/2.
```

Reordering a child pair swaps both projector numerator and angular
denominator. Reflection preserves the scalar q/E/S_star law and cos4;
either canonicalization is harmless only when performed consistently.
The identity holds on all corresponding simple-root branches of this
specified finite positive family, not just at t=0.

This is an exact endpoint/source dictionary, not an integration of a
fluctuating sublattice away from its filled endpoint, a full RG flow,
a Jordan structure or a continuum-field identification. It supplies no
sign prediction for a new response and does not rescue or extend the
stopped F4-only block.
