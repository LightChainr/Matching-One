# The actual S4 [2,2] seam trace is confined to topology rank1

**Exact finite result.** At the already specified m=2, Q=4 closed-source
law, three ordinary colour-permutation seam partitions give the canonical
S4 `[2,2]` character projection. On an honest square-cell torus its q and
E numerators vanish configuration by configuration. Its partition term
need not vanish: the projection has support in topology rank1 only.
Consequently its possible effect on original U is a normalization
transmission, not a direct nonzero q/E numerator.

This is an actual finite torus closure, not the regular all-ones endpoint
or an arbitrary multiple of a categorical character. The conclusion uses
the existing integral saturation and disjointness properties of occupied
components. It makes no CFT-field identification. The notation `[2,2]`
denotes an S4 irrep; topology rank2 is a different concept.

## 1. Fix the seam and all three positive local sums

Use the period basis `(a,b),(-b,a)` of each of the two fixed geometries.
Across the first period apply a colour permutation pi in S4; across the
second period apply the identity. The local vacant/active interaction,
activity and ambient-rank multiplier `2^-r` are unchanged. Call the
resulting occupation sum Z_pi. It is nonnegative for each pi and includes
all occupations, not just the seam-crossing patterns.

For a fixed occupied NN component C with winding image H_C, its allowed
colour count is

```text
f_pi(C)=#Fix{pi^u : (u,v) in H_C}.                            (1)
```

Contractible components have four colours. With reduced odds activity
`a0=y/2^5`, the full expression is

```text
Z_pi=sum_A a0^K 4^B 2^-r product_(occupied components C) f_pi(C).
                                                                    (2)
```

This is the local colour gas with one spatial seam and the same stipulated
topological factor as before. Equivalently expand a finite transfer
history and close it with the global colour action of pi, retaining the
rank weight in the closure. For q/E numerators insert their actual rank
values in (2). The identity seam recovers the original partition.

## 2. The central projector has three nonzero class coefficients

On the unordered distinct-pair carrier, the `[2,2]` character is the
number of fixed unordered pairs minus the number of fixed colours. This
subtracts the singlet and standard summands of that carrier. Directly
counting gives

| Class | Size | Fixed colours | Fixed pairs | chi_[2,2] |
|---|---:|---:|---:|---:|
| identity | 1 | 4 | 6 | 2 |
| transposition | 6 | 2 | 2 | 0 |
| double transposition | 3 | 0 | 2 | 2 |
| three-cycle | 8 | 1 | 0 | -1 |
| four-cycle | 6 | 0 | 0 | 0 |

The centre of the group algebra therefore gives the fixed projection

```text
Z_[22] = (1/6) Z_id + (1/2) Z_(12)(34) - (2/3) Z_(123).        (3)
```

Indeed the projector is `(dim/|S4|) sum_pi chi(pi^-1) U(pi)` with
dim=2 and |S4|=24. The same linear functional acts on the rank numerators
and on their thermal derivatives. The coefficient of Z_id in (3) is
1/6, not the character ratio of one selected permutation. No transposition
measurement alone would reconstruct this full projection.

Equation (3) is a signed trace functional of positive seam sums. Its
individual configuration contribution is not a probability; no assertion
of positivity for every projected connection pattern is needed.

## 3. Rank0 and rank2 disappear for different exact reasons

In topology rank0 every occupied component has H_C=0. Every Z_pi term
then has the same colour multiplicity `4^C_B`. The coefficients in (3)
sum to zero. Thus the character-filtered contribution vanishes on each
such configuration.

In topology rank2, there is exactly one occupied component carrying both
torus directions, and its winding image is all Z² by the repository's
saturation theorem. All other occupied components are contractible:
two disjoint components cannot carry winding cycles of nonzero algebraic
intersection. The distinguished component has respectively4,0,1 colours
under the three seams in (3). Its coefficient is

```text
(1/6)*4 + (1/2)*0 - (2/3)*1 = 0.                             (4)
```

Representation-theoretically that component closes in the point-colour
representation V, which contains only a singlet and the standard irrep,
not `[2,2]`. All other component factors are common. This proves the
second configurationwise zero.

Since q=E=0 on rank1, all configurations now obey

```text
N_q,[22]=N_E,[22]=0,                                         (5)
```

and their complete thermal derivatives vanish as well. The argument
uses the original topology and rank weight, not endpoint annihilation.
It applies to the two fixed honest square quotients and to any quotient
in the same proved saturation scope.

## 4. The surviving rank1 weight is computable without colour enumeration

Let b2(A) indicate whether at least one occupied component has a winding
whose first deck coordinate is nonzero modulo2. Let n3(A) count components
with such a coordinate nonzero modulo3. For a double transposition,
an odd winding leaves no allowed colour, while an even one leaves four.
For a three-cycle, a nonzero winding modulo3 leaves one colour instead
of four. Relative to the original identity-seam occupation weight, (3)
therefore becomes

```text
beta(A)=1/6 +(1/2) 1_(b2=0) -(2/3) 4^(-n3),
Z_[22]=sum_A w_star(A) beta(A).                               (6)
```

The deck coordinates are integral in the specified period basis, not
Cartesian displacements reduced modulo2/3. Per-component constraints are
necessary: ambient rank alone does not determine n3.

For two separated essential rows crossing this seam once, b2=1 and n3=2,
so beta=1/8. This agrees with the independent physical
[two-cluster trace](closed-source-two-winding-cluster-trace.md): two of
the16 ordered colour states belong to `[2,2]`. Thus (5) does not imply
Z_[22]=0.

## 5. The remaining finite question has one prescribed score

Use the unique infinitesimal insertion determined by (6),
`w_epsilon=w_star(1+epsilon beta)`. It is well-defined near epsilon=0,
since beta is bounded. It varies the coefficient of this exact trace
packet without fitting a new microscopic combination. Equation (5) gives

```text
delta Z_g=Z_[22],g,       delta N_q,g=delta N_E,g=0.            (7)
```

Set f_g=Z_[22],g/Z_g and use the already derived
[rank1 transmission functional](closed-source-rank-one-trace-transmission.md).
The complete response is

```text
delta U/A_N = C_c (f_c)_z + C_dz (f_d)_z + C_d f_d.            (8)
```

Every coefficient and input is now defined by the actual lattice closure.
The unknown is a number, not an unspecified projector-to-observable map.
The frozen [Q4 calculation](../analysis/p337_s4_trace_transmission_contract.json)
evaluates (8) once at the saved original m2 pooled root on the N25 pair.
It adds only the missing per-component modulo2/3 closure information.

A nonzero answer would establish finite trace-to-U normalization
transmission despite the exact direct-numerator zero (5). A zero answer
would retain the trace but stop this transmission explanation at the fixed
point. Neither outcome assigns a Q=1 activated amplitude, a continuum
four-leg field, or a scaling exponent. Those require the further generic-Q
and spectral interface; no colour number or seam is changed to rescue the
fixed calculation.

**Completed:** the one frozen score now gives
`V_beta=+5.440121494634842e-6`, with a strictly positive exact rational
enclosure. The [result](closed-source-s4-trace-transmission-result.md)
therefore establishes this finite normalization route; it does not supply
a Q1 derivative or a continuum field assignment.

## Scientific card

- **Changed mechanism space:** the canonical Q4 colour projection is
  neither a generic endpoint nor a free trace ansatz. Its direct q/E
  channel vanishes exactly, leaving a precisely specified normalization
  route whose original-U coupling can be decided.
- **Observer/sector/source/geometry:** original pooled-root U; S4 `[2,2]`
  central insertion along first deck seam; Sstar at m2; N25(5,0)/(4,3).
- **Dependency:** exact topology/colour algebra plus a named insertion
  in the same finite configuration population. No independent data vote.
- **Boundary:** no identification of this entire finite isotypic trace
  with one particular CFT state; no real-Q derivative inferred from three
  fixed-Q seam evaluations.
