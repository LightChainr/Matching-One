# Square-site Euler localization of the rank source: the Mertens--Ziff integrand as `r-1`

Date: 2026-09-14

Status: exact finite configuration identity, but not a novelty claim.  The configuration-level formula is Mertens--Ziff 2016, Eq. (11), derived from Euler's Gem for finite matching lattices.  The useful step here is to identify its three-valued right-hand side exactly with the ambient-homology rank charge `X=r_4-1` used throughout Matching One, and to exponentiate the identity into a source dictionary.

Literature boundary: Mertens & Ziff, *Percolation in Finite Matching Lattices*, arXiv:1603.07289v2, §II.

## 1. The exact configuration identity

Consider an honest square torus with a black site configuration `omega`.

Define

```text
k4(omega)      : number of black NN clusters,
k8(omega^c)    : number of white matching (NN+NNN) clusters,
K               : number of black sites,
E               : number of black NN edges (both endpoints black),
F0              : number of elementary square faces whose four corners are black,
r4              : rank im[H1(black NN complex)->H1(T^2)] in {0,1,2}.
```

Mertens--Ziff Eq. (11) states configurationwise

```text
k4-k8-(K-E+F0)
 = +1  if black has a cross-wrapping cluster,
 = -1  if complementary matching white has a cross-wrapping cluster,
 =  0  otherwise.
```

Their classification also records that single/spiral wrapping occurs in paired black/white components and contributes zero to this difference.

But these three cases are exactly

```text
r4=2 : black cross -> +1,
r4=0 : white cross / black no ambient homology -> -1,
r4=1 : one-dimensional ambient homology -> 0.
```

Therefore

```text
boxed:
r4(omega)-1
 = k4(omega)-k8(omega^c)-K+E-F0.                    (1.1)
```

This is the square-site analogue of the bond/FK Euler/Krushkal localization, with the important difference that the matching white graph and the elementary black-face term are intrinsic to the site construction.

Independent implementation check during this audit: lifted NN homology, black NN clusters, white G8 clusters and `(K,E,F0)` were recomputed for every configuration at honest `L=3` (512 states) and `L=4` (65,536 states); (1.1) had zero violations.  This is a regression, not the proof; the published Euler argument is the all-size source.

## 2. The familiar matching function is literally `E[r-1]`

Take expectation under Bernoulli(p).  Since

```text
E K/L^2  = p,
E E/L^2  = 2 p^2,
E F0/L^2 = p^4,
```

we obtain

```text
E[r4-1]
 = N_L(p)-Nhat_L(1-p)-L^2[p-2p^2+p^4].              (2.1)
```

The left side is

```text
P2-P0,
```

because `r-1` takes values `-1,0,+1`.  Thus the Matching-One topological balance observable and the finite Sykes--Essam/Mertens--Ziff matching function are not merely asymptotically related or two different estimators:

```text
boxed:
M_L(p)=P2-P0=E[r-1]
       = cluster-count difference - local Euler polynomial.    (2.2)
```

This is exactly the dictionary that should be used whenever cluster-number and rank-language branches meet.

## 3. The entire bounded topological source localizes, not only its derivative

Because (1.1) is configurationwise, for every real/complex source `h`,

```text
exp[h(r-1)]
 = exp[h k4]
   exp[-h k8]
   exp[-h K]
   exp[+h E]
   exp[-h F0].                                      (3.1)
```

Hence the exact finite source partition function

```text
Z_L(p,h)=E_p exp[h(r-1)]
        =P0 e^-h + P1 + P2 e^h
```

can also be read as a **two-colour cluster gas with local Euler interactions**:

```text
black NN cluster fugacity        : e^h,
white matching cluster fugacity  : e^-h,
black site local factor          : e^-h per occupied site,
black NN-pair factor             : e^+h per occupied NN edge,
black plaquette factor           : e^-h per all-black elementary face.
```

Multiplying the Bernoulli weight explicitly,

```text
p^K q^(N-K) e^{h(r-1)}
 = q^N
   [(p/q)e^-h]^K
   e^{hE-hF0}
   e^{h k4-h k8}.                                   (3.2)
```

So the source is not a mysterious closure label.  It is an exact coupled black/white cluster fugacity deformation plus finite-range site/edge/plaquette factors.

The cluster fugacities are nonlocal in a spin Hamiltonian sense, but are local in a connectivity transfer: they can be paid when a component retires/closes, exactly as in random-cluster transfer methods.

## 4. Exact complement action on the sourced family

The digital-Alexander/matching identity gives

```text
r8(omega^c)=2-r4(omega),
```

so

```text
X_hat=r8-1=-(r4-1)=-X.
```

Thus on the doubled square-site matching pair the rank source transforms exactly as

```text
(p,h,G4) <-> (1-p,-h,G8).                            (4.1)
```

This is a genuine involution on the **doubled sourced finite model family**.

It still does not prove that every local continuum scaling field has a scalar matching parity.  The source `h` is a global/topological deformation represented by cluster fugacities.  But it supplies one exact tangent direction whose pair-exchange action is known without CFT assumptions.

## 5. Relation to the bond/FK Krushkal localization

For bond/FK states the separate exact identity is

```text
r-1 = k(A)-k(A*)+|A|-|V|,
```

which turns the source into primal/dual cluster fugacities and an edge factor.

The site formula (1.1) is structurally parallel:

```text
bond/FK : cluster imbalance + edge Euler term;
site    : black/white matching cluster imbalance + site-edge-face Euler term.
```

Therefore the generic research idea “topological charge is an Euler imbalance between two complementary cluster gases” applies to both settings, although the microscopic local terms differ.

This corrects an earlier overly sharp boundary that treated Euler localization as a bond/FK-only structural option.

## 6. A new transfer/source route for the square-site master object

The exact source object already used in the master scaling proposal is

```text
Z_L(p;h,z)=P0 e^-h + P1 H_p(z)+P2 e^h.
```

Equation (3.2) supplies a concrete square-site implementation for the `h` direction:

1. retain black NN and complementary white matching component closure counts;
2. add local weights for black site, NN edge and full plaquette motifs;
3. multiply retiring black/white components by `e^h/e^-h` respectively.

The rank-one neutral source `z` remains a separate count of parallel essential components.  Thus a two-source transfer no longer needs to treat `h` as an external after-the-fact rank label.

A practical implementation should first reproduce the ordinary `h=0` safe/rank transfer and finite small-torus source polynomial before any massive/continuum interpretation.

## 7. Relation to source normalization and #802

This Euler source makes the normalization issue explicit.  If the local factors in (3.2) are inserted as unnormalized row weights, the physical free energy must include the row partition normalizer before comparing source derivatives.

Moreover the local source contains degree-0/1/2/4 Bernoulli-chaos pieces plus cluster-count terms.  Profiling the thermal direction therefore cannot be done by inspecting one raw motif derivative.  The normalization-safe/source-quotient rules from the #802 audit still apply.

The advantage is conceptual: the particular combination that equals `r-1` is fixed configurationwise, so all nuisance pieces are tied together by an exact identity rather than chosen ad hoc.

## 8. New continuum question

The useful continuum question is no longer simply

```text
which local CFT field is matching-odd?
```

but rather

> how does the exact finite Euler-imbalance source `h` decompose under RG into thermal/nuisance directions, topological torus sectors, and irrelevant angular/scalar corrections?

The leading root correction is a correction to the zero of

```text
partial_h log Z_L(p,h)|_(h=0).
```

This is a source-defined object.  A continuum candidate is relevant only if it has nonzero overlap with this particular sourced response after the common/thermal pieces have been removed.

This formulation avoids assigning an OPE parity before the source-to-RG map is built.

## 9. Claim boundary

Exact / published:

- Mertens--Ziff configuration Euler identity;
- its three topological cases;
- identification of those cases with `r-1` under the rank convention;
- exponential source rewrite (3.1)--(3.2);
- complement action `h->-h` on the doubled sourced family.

New synthesis / programme:

- use the Euler-localized source directly inside the square-site transfer;
- study its RG decomposition as the primary source-defined route to the noncommon correction;
- combine it with the neutral source `z` in the master topological generating object.

This does not establish local continuum operator parity or a new percolation theorem; it turns a known finite matching identity into the exact microscopic source dictionary needed by the current research programme.