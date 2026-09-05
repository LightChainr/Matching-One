# Jordan identifiability: a closure obstruction and a C4-protected escape

Status: proved conditional model statements with exact rational controls;
application to Matching One is a proposed discriminator, not a data result.
Baseline: `4d70c1787ff97dbb98cb5e96022f947bb8fad97e`.
Related: #218, #249, #250, #370; PR #273.

## 1. Why this is the next question, not another scalar fit

PR #273 rejects its frozen analytic q=2 model and leaves its Jordan alternative
close to the acceptance boundary. That is a comparison of declared models,
not an exclusion of every ordinary two-mode model with arbitrarily small
spectral splitting.

The later P250 radius-five reveal, commit `11130ae`, has already selected
`Alexander_R2_conjugation` over the other four frozen maps at alpha 0.01:
R2 has p=0.013379 and identity has p=0.000728. The selected object is a
truncated annihilator-line morphism, not a closed five-state translation
algebra. This note neither reopens that vote nor treats rank-five compatibility
as established flatness. Its application begins *after* a higher-degree flat
extension and a correctly typed symmetry action have been obtained.

Issue #370 already allows a spectral-separation condition. The first result
below explains why such a condition, or a different exact obstruction, is
mathematically necessary for a robust Jordan-versus-ordinary certificate.
The second result supplies one explicit exact symmetry obstruction.

## 2. A common ordinary semigroup approaches the Jordan semigroup

Put

```text
N = [[0,1],[0,0]],
G_epsilon = [[0,1],[epsilon^2,0]],  epsilon > 0.
```

Then `N^2=0`, `N!=0`, whereas `G_epsilon` has two distinct real eigenvalues
`+epsilon,-epsilon` and satisfies `G_epsilon^2=epsilon^2 I`. Hence

```text
T_epsilon(t) = exp(t G_epsilon)
 = cosh(epsilon t) I + sinh(epsilon t)/epsilon * G_epsilon
 -> I+tN = T_0(t).
```

Entry by entry, the discrepancy is O(epsilon^2), uniformly for t in a bounded
interval. In particular the upper-right discrepancy is
`sinh(epsilon t)/epsilon-t = epsilon^2 t^3/6+O(epsilon^4)`;
the lower-left entry is `epsilon sinh(epsilon t)=epsilon^2 t+O(epsilon^4)`.
For every epsilon the full family obeys exactly

```text
T_epsilon(t+s)=T_epsilon(t) T_epsilon(s),
[T_epsilon(t),T_epsilon(s)]=0.
```

Multiplying by a fixed leading character and setting `t=log Norm(m)` gives
ordinary powers with nearby exponents approaching the same logarithmic
Gaussian-cover law. This is one common diagonalizable realization across
all multipliers, not independently refitted scalar curves.

**Finite-observation closure statement.** Fix finitely many times/words and
fixed source/readout vectors. Their Jordan prediction vector is in the closure
of the prediction vectors of this ordinary two-mode family. Consequently a
positive-radius confidence set containing that Jordan point in its relative
interior also contains ordinary predictions, provided those approximants obey
all other hard constraints and stay in the same exact covariance-support
subspace. This follows simply from continuity and the definition of relative
interior. Boundary-only intersections require separate treatment.

Therefore there is no positive-distance certificate excluding the entire
unrestricted ordinary class while retaining an interior Jordan point. A
numerical optimizer's failure to find the near-colliding alternatives cannot
supply the missing separation. In particular, commutation and composition
alone do not resolve this counterexample.

This does **not** obstruct exact noiseless identities, or an ordinary class
with an independently justified spectral gap, normalized eigenbasis condition
bound, physical spectrum restriction, positivity condition, or symmetry
constraint that actually excludes the approximating family. It does not
invalidate the specific frozen q=2 rejection in PR #273.

The common eigenbasis can be chosen as

```text
S_epsilon = [[1,1],[epsilon,-epsilon]],
condition_inf(S_epsilon) = 1+1/epsilon,  0<epsilon<=1.
```

The matrices remain bounded while their eigenbasis degenerates. A bound on
matrix entries alone is not a bound on this conditioning. Such bounds must
also have a declared normalization/gauge and physical justification.

At fixed dynamic range and with nonzero quadratic sensitivity, mean separation
is O(epsilon^2). Under the usual n^(-1/2) sampling error model, resolving this
particular local family costs order epsilon^(-4) samples. This is a conditional
resolution scaling, not a universal percolation sample-complexity theorem.

## 3. Independently verified rational finite-word controls

`scripts/jordan_nonseparation.py` stores six ordinary realizations with
`epsilon=10^-1,...,10^-6`. The synthetic generators are

```text
A_epsilon(lambda) = I+lambda G_epsilon,
lambda_a=1, lambda_b=3.
```

They share the exact eigenbasis S_epsilon and have distinct positive real
eigenvalues. All 63 words of length at most five are verified in each of the
six realizations. The verifier reads the stored rational matrices without
calling the builder or a numerical optimizer. It checks 378 word-error bounds.

A telescoping product proves, for |epsilon|<=1,

```text
||A_w(epsilon)-A_w(0)||_inf <= epsilon^2 C_w,
C_w = sum_j |lambda_j| product_(k!=j) (1+|lambda_k|).
```

An additional independent binomial calculation gives the b^5 error
`360 epsilon^2 + 648 epsilon^4` on the declared epsilon range. At epsilon
`10^-6`, the maximum stored context error is less than `10^-9`, despite every
approximating generator being ordinary and simultaneously diagonalizable.

The numbers 1 and 3 are rational synthetic controls, **not** log 2 and log 5.
The exact physical radial-semigroup closure statement is the analytic proof
in Section 2, not an interpretation of these rational constants.

## 4. Why passing from roots to a spectral scheme does not alone solve it

The quotient-algebra formulation already proposed in #250 is the right
basis-independent language. But without extra constraints, nonreducedness is
not automatically noise-separable from reducedness.

For a smooth complex surface, the Hilbert scheme of length-n subschemes is
irreducible and contains the open locus of n distinct points. Thus any finite
cyclic quotient of C[X,Y], or of the Laurent algebra on (C*)^2, can be approached
by reduced quotients in the unrestricted Hilbert scheme. In a local frame of
the tautological bundle, multiplication matrices and finite moments vary
continuously. This is the higher-rank analogue of the explicit two-mode limit.

This uses the classical smooth-surface theorem of Fogarty; a current account
is the introduction of Craw--Yamagishi, arXiv:2607.08913. Their additional
results about canonical singularities are not needed here. Reducedness of a
*Hilbert scheme as a moduli space* is also not reducedness of each subscheme
that it parametrizes.

Crucially, a smoothing need not preserve a specified finite-group character.
An equivariant Hilbert scheme/fixed-character locus is not the Hilbert scheme
of the quotient surface. This distinction makes the following escape possible.

## 5. C4-character obstruction for a cyclic length-five Laurent quotient

### Exact hypotheses

Let U,V be commuting invertible complex matrices acting on a five-dimensional
space. Suppose b is cyclic for their Laurent algebra and R obeys

```text
R^4=I,
R U R^-1=V,
R V R^-1=U^-1,
R b=b.
```

Cyclicity identifies the space with a length-five quotient
`A=C[U^+-1,V^+-1]/I`; the invariant source fixes its algebra-unit convention.
The quotient must be genuinely closed/flat, not only a rank-five truncated
moment window. No assumption about these hypotheses having been measured in
P250 is made in this note.

### Theorem

If A is reduced, its rotation character must be

```text
(tr I, tr R, tr R^2, tr R^3) = (5,1,1,1),
(m_1,m_i,m_-1,m_-i) = (2,1,1,1).
```

**Proof.** A reduced A is the function algebra of five distinct joint spectral
points in (C*)^2. The declared rotation acts on that support by
`(lambda,mu) -> (mu,lambda^-1)` (using the inverse convention gives the same
orbit counts). There are exactly two fixed points, `(1,1)` and `(-1,-1)`.
There is exactly one length-two orbit, `{(1,-1),(-1,1)}`. Every remaining orbit
has length four. A reduced invariant five-point set must therefore consist of
one fixed point and one four-point orbit: with at most two fixed points and
one pair, no other decomposition of five exists.

Because the cyclic source is the invariant unit, R acts by permutation on the
function algebra. A fixed point contributes the trivial character and a
four-orbit contributes the regular C4 character. The stated trace and
multiplicity vectors follow. QED.

A different verified character therefore excludes **all** reduced length-five
quotients satisfying these exact hypotheses, not just one fitted root model.
Finite-group character multiplicities are locally constant in a flat
equivariant family, so the exclusion cannot be evaded by making eigenvalues
arbitrarily close while preserving the specified character.

Matching the allowed character is necessary, not sufficient, for reducedness.
If the cyclic source has a nontrivial scalar character, its known character
must first be divided out, or the multiplicities shifted accordingly. A
multi-source/noncyclic representation requires a different theorem.

Use the global **Laurent** rotation law. Replacing inverse translation by
minus a translation matrix is not valid. The linear law X->Y, Y->-X used
below is justified only for logarithms in the explicit unipotent control.

## 6. A nonreduced control with the forbidden character

Consider the exact five-dimensional algebra

```text
A_* = C[x,y] / ((x,y)^3, x^2+y^2),
basis = (1,x,y,x^2,xy).
```

Let X,Y be multiplication by x,y; let R send x to y and y to -x, fixing 1.
The ideal is invariant, X and Y commute, all cubic monomials vanish, and
X^2 is nonzero. Define

```text
U=exp(X)=I+X+X^2/2,
V=exp(Y)=I+Y+Y^2/2.
```

These are invertible and satisfy exactly the global Laurent rotation laws in
Section 5. The source 1 is cyclic and invariant. Direct calculation gives

```text
(tr I,tr R,tr R^2,tr R^3)=(5,-1,1,-1),
(m_1,m_i,m_-1,m_-i)=(1,1,2,1).
```

This is a concrete symmetry-protected nonreduced quotient: unrestricted
smoothing is possible, but C4-equivariant smoothing to five distinct spectral
points with the same unit/source character is impossible.

`scripts/c4_five_state_obstruction.py` records the exact rational matrices,
a separate ordinary five-point positive control, and 20 exact construction
checks. Its separate verifier rechecks the stored translations, rotations,
cyclicity, nilpotency, trace characters and orbit arithmetic without calling
the builder. Tests deliberately corrupt matrices, characters and support to
ensure rejection.

This control has nilpotency order three. It is not asserted to be the actual
P250 state, a rank-two Jordan pair, a CFT module, or a new percolation field.

## 7. The next discriminator for Matching One

First complete the higher-degree flat/border-basis test in the already chosen
P250 R2 gauge. Then establish whether the source is cyclic and whether a
rotation acts internally on a five-dimensional quotient with the exact
Laurent relations above. The present two-charge, two-hand data may require
another source character or an enlarged module; that is a legitimate outcome,
not a reason to force the theorem onto the data.

If the hypotheses do hold, inspect the four exact symmetry projectors

```text
P_j=(1/4) sum_(k=0)^3 i^(-jk) R^k,  j=0,1,2,3.
```

Their ranks are the discrete character multiplicities. An allowed-character
failure is a finite algebraic obstruction to ordinary reduced modes; an
allowed-character pass leaves ordinary and confluent alternatives unresolved.
A noisy fitted R must not simply be projected onto C4 and then reported as an
exact experimental symmetry certificate. Carry uncertainty through a typed
feasibility test, or derive the action from exact microscopic symmetry.

For #370, compile the closure test before generic SOS: state which hypotheses
make the ordinary class a separated set. With unrestricted near-collisions,
report a resolution/spectral-gap frontier. With verified exact symmetry,
try the character obstruction before optimizing over root locations.

The shift is from "which fit looks logarithmic?" to "what prevents every
ordinary realization from approximating the declared observation?".

## 8. Reproduction and claim boundary

```bash
python -m compileall -q scripts tests
python scripts/jordan_nonseparation.py --verify results/jordan-nonseparation/latest.json
python scripts/c4_five_state_obstruction.py --verify results/jordan-nonseparation/c4-five-state.json
python -m unittest discover -s tests -p 'test_jordan_nonseparation.py' -v
python -m unittest discover -s tests -p 'test_c4_five_state_obstruction.py' -v
```

The new isolated bundle passed 28 tests on CPython 3.13. No Monte Carlo was
run, no production data were refitted, and the complete repository suite was
not run. The exact proofs and rational controls concern the declared model
classes only. Their application to Matching One remains C0 until the listed
flatness, source and symmetry hypotheses are established. This work does not
identify x=21/4, prove a logarithmic continuum overlap, or determine p_c.

## Sources and existing work

- P250 radius-five result: `notes/p250_projective_leg_radius5_morphism_result.md`
  at commit `11130ae`; issue #250 comment `5467890317` (2026-08-30 09:28:59 UTC).
- #218: emerging-Jordan diagnostics; #249: context-Hankel realizations;
  #250: quotient/annihilator formulation; #370: proof-carrying model elimination.
- Liu--Jacobsen--Saleur, *Emerging Jordan blocks in the two-dimensional Potts
  and loop models at generic Q*, arXiv:2403.19830. Existing physical motivation,
  not a claim that a logarithmic scalar fit determines a Jordan module.
- Craw--Yamagishi, *Hilbert schemes of points on canonical surfaces*,
  arXiv:2607.08913v1 (2026-07-09), Introduction: smooth-surface Fogarty theorem.
  The explicit C4 orbit proof above is self-contained and does not identify
  an equivariant Hilbert scheme with the Hilbert scheme of a quotient surface.
- Camia--Feng, *The percolation energy field and its logarithmic partner*,
  arXiv:2508.16047v2 (2026-06-01): actual triangular-site lattice logarithmic
  fields, a distinct physical identification task from the present controls.
