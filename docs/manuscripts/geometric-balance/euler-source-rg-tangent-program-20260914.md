# From the exact Euler source to an RG tangent map

Date: 2026-09-14

Status: programme built around an exact finite sourced-family involution.  It partially answers the structural concern of #61 without claiming a local CFT/OPE automorphism.

## 1. An exact two-parameter sourced family exists on the doubled site model

For the black NN / white matching pair define the exact topological source

```text
X=r4-1=k4-k8-K+E-F0.
```

The sourced finite weight is

```text
P_{p,h}(omega)
 proportional to
 p^K (1-p)^(N-K) exp[hX].
```

On the doubled model family, complement/matching gives exactly

```text
(G4,p,h,omega)
 <->
(G8,1-p,-h,omega^c).
```

Thus, unlike an assumed local CFT parity, this sourced finite-family exchange is an actual involution.

## 2. Center the thermal coordinate

Let `p_c^4` and `p_c^8=1-p_c^4` be the paired critical points.  Use a local thermal coordinate `t` chosen so complement sends

```text
t -> -t.
```

At linear order one can take `t=p-p_c` up to a positive metric factor, provided the two members use paired normalizations.

Then on the exact two-dimensional sourced family

```text
(t,h) -> (-t,-h).                                      (2.1)
```

So the doubled matching action on the span of these two finite tangent directions is explicit.

This is a small but genuine piece of the RG-tangent map requested by #61.

## 3. What this does and does not establish

It establishes the microscopic transformation of two **source coordinates**:

```text
thermal complement t,
topological Euler source h.
```

It does not establish that a local scaling field with a given dimension/spin is an eigenvector of matching exchange.  Under RG, the microscopic `h` direction may mix with every continuum field allowed by the same geometric/source quantum numbers.

The correct question is therefore to measure/derive the image of the `h` tangent under coarse graining.

## 4. Add local residualized sources to enlarge the tangent space

Choose a finite local source basis `g_a` after normalization and thermal residualization.  Bernoulli chaos gives a convenient exact microscopic grading:

```text
degree 2 even control,
degree 3 odd-under-complement control,
selected D4 spatial representations.
```

The full finite source coordinates are

```text
u=(t,h,g_2,g_3,...).
```

The exact pair exchange acts linearly on these microscopic source functions at the baseline (with the known chaos signs after graph/color pairing).

The RG question becomes:

```text
what is the large-scale linear map from this microscopic source basis
into the low-dimensional scaling-response space?
```

No OPE algebra is needed for this first stage.

## 5. Observable response matrix

For a declared set of finite observables / safe free energies `O_i`, compute

```text
R_ia(L)=partial_(g_a) O_i
```

with normalized physical source weights.

Useful rows include

```text
Euler/rank charge response,
fixed-b rank shape response,
safe-sector free-energy difference,
angular H4/H0 projected root response.
```

Useful columns include

```text
thermal t,
topological source h,
degree-2 residualized motif,
degree-3 residualized motif,
controlled anisotropy source.
```

The singular vectors / stable column relations of `R(L)` as `L` grows are empirical RG tangent directions.  Only after a stable block appears should one attach continuum field labels.

## 6. The topological source is especially informative because it is bounded

Although its Euler representation contains extensive cluster/local terms, configurationwise

```text
X in {-1,0,+1}.
```

Thus `h` couples directly to the exact root-information bottleneck rather than to a generic bulk density.

The mixed response

```text
partial_g partial_h log Z |_(h=g=0)
 = Cov(X,H_g)
```

is exactly the numerator controlling root motion under source `g`.

At the root, with thermal score `S_p`,

```text
T_g=-Cov(X,H_g)/Cov(X,S_p).
```

Therefore the entire #802 root-tangent programme can be viewed as measuring mixed susceptibilities with the exact Euler source.

## 7. A clean interpretation of “first noncommon correction”

Instead of saying a continuum field is “matching odd,” define the observable fact:

```text
a perturbation/source has a nonzero mixed susceptibility with X
(or the corresponding safe-phase Euler/topological free-energy difference).
```

Then ask which angular/radial scaling block carries that mixed susceptibility.

For the current square model the leading block appears H4-like with root exponent four.  A later angular-scalar block may have exponent seven.  These are statements about the response of the exact `h`-sourced observable, not parity assignments of abstract fields.

## 8. Relation to generic-Q / bond-FK source

The bond/FK Krushkal/Euler source provides an analogous exact topological tangent, with its own local cluster-gas dictionary.  At Q=1 the continuum limits of the bond and site realizations may belong to the same topological source class, but that is a universality/interface statement, not an exact lattice equality.

This suggests a future high-value cross-model test: compare normalized mixed `h-g` response ratios between square-site and a rigorously controlled bond/FK realization after matching thermal metrics and modulus.

## 9. Acceptance levels for #61-style parity language

Use three levels:

```text
Level 1: exact microscopic source exchange
         (established for t,h and Bernoulli-chaos source functions);

Level 2: RG tangent block / scaling-response exchange
         (to be inferred/proved from response matrix or transfer algebra);

Level 3: local OPE/interchiral automorphism
         (not required for finite-size selection and currently unproved).
```

Most current Matching-One conclusions only need Level 2, not Level 3.

## 10. Claim boundary

Exact:

- finite sourced-family involution `(t,h)->(-t,-h)` on the doubled site pair;
- mixed susceptibility `Cov(X,H)`;
- residualized root response formula.

Programme:

- construct the empirical/theoretical RG tangent map from a small source-response matrix;
- identify continuum blocks only after stable angular/radial/source structure emerges.

This is a constructive replacement for prematurely assigning a scalar matching parity to continuum fields.