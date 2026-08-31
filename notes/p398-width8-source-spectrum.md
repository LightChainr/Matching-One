# P398: two protected readout rays, but 93 propagation directions per ray

Parent: `dbd408154b4215ca41fbf26c0fd962997074f05d`.
This is a deterministic width-eight calculation in the same positive
continuous square-bond frontier process. It uses 1,430 noncrossing states,
not row-by-row bond enumeration or Monte Carlo.

**Outcome:** the two old readout rays survive exactly by Kreweras duality,
but their single-mass kernels do not. Each ray generates a 93-dimensional
propagation subspace. The formerly fast ray now contains the smallest mass.

## Same physical readouts, same cyclic character

At widths four and eight use

\[
G=\sum_j(J_j-I)+\sum_j(D_j-I),\qquad
A=\sum_j i^j\mathbf1\{j\sim j+1\},\quad
L=\sum_j i^j\mathbf1\{j\text{ is singleton}\}.
\]

These are exactly the previous AP and landing Gram readouts, not replacement
observables. Joining a configuration to the AP reference gives block count
`b-1+1{j~j+1}`. Joining it to the one-singleton landing reference gives
`1+1{j singleton}`. Their constant terms vanish in the Fourier sum.
The one-site cyclic character remains i: momentum index is k=1 at width four
and k=2 at width eight. Width and circumference change; lattice momentum
pi/2 and the definitions of the two local statistics do not.

An orbit basis gives an **exact complex dimension 186** for the width-eight
character-i sector. The block Krylov ranks from `(A,L)` are
`2,4,6,...,186`. These are not floating-point rank calls: integer/Gaussian
integer coefficients are reduced modulo the prime 65537 with `i=256`.
A full-rank minor modulo that prime is nonzero over Q(i), and 186 is also
the exact orbit-count upper bound. Hence the source-generated dimension is
exactly 186. This says nothing about the number of continuum fields.

## Why closure fails: two explicit current-configuration statistics

Let

\[
T_2=\sum_j i^j\mathbf1\{|\text{block}(j)|=2\}.
\]

If j and j+1 belong to different blocks, let n_j count the adjacent edges
whose two endpoints are in precisely those two blocks. Define

\[
R=\sum_j i^j\mathbf1\{j\not\sim j+1\}(n_j-1).
\]

The exact generator identities are

\[
\boxed{GL=-3L+T_2,\qquad GA=-3A+R.}
\]

T2 counts the extra way another site's detachment can make j a singleton.
R counts extra neighboring contacts through which a join can reconnect the
two endpoint clusters. Thus the new directions are size-two membership and
boundary-contact multiplicity, not abstract fitting parameters.

At width four these collapse to
`T2=-(1+i)A-2L` and `R=-2A-(1-i)L`. At width eight,
`rank(A,L,R,T2)=4` exactly, from the first Krylov increment.
An immediate counterexample is

```text
partition 00000000: A=L=0, (-GA,-GL)=(0,0)
partition 00000011: A=L=0, (-GA,-GL)=(0,1+i).
```

No fixed two-by-two generator on the recorded A/L values can represent both.
Even the optimal stationary-L2 projection leaves leakage covariance
approximately `[[.60031444,.12469100(1+i)],
[.12469100(1-i),.60031444]]`.

## What survives: a duality-protected two-ray decomposition

Represent a noncrossing partition by its increasing block-cycle permutation
p, and define its Kreweras complement from the cycles of `p^-1 c`, with
`c=(0 1 ... w-1)`. On the complete finite state sets the script checks exactly

```text
K^2 = rotation by -1;   KG = GK;
A(K configuration) = -i L(configuration);
L(K configuration) = A(configuration).
```

Therefore the unchanged combinations

\[
\psi_-=(A-e^{-i\pi/4}L)/\sqrt2,\qquad
\psi_+=(A+e^{-i\pi/4}L)/\sqrt2
\]

have different Kreweras characters `-exp(-i*pi/4)` and `+exp(-i*pi/4)`.
The unique stationary measure is K invariant, so their cross-correlation
vanishes at every distance. This is **two protected readout rays**, not two
transfer states. K has trace zero on the 186-dimensional sector, and K squared
is -i there, giving exact ray-sector dimensions 93+93. Separate modular
Krylov ranks attain 93 from each respective single source.

The named-channel mixing invariant `I_c=1` thus survives. The width-four
mass relation `25 I_m^2(2-I_c)=2` does not: the projected two-channel kernel
contains many masses and is not a semigroup.

## Actual low spectrum and the most useful source projection

The numerical stationary/spectral solve uses float64 in the exact orbit
block. Its spectrum is real to 4.6e-15; the complete spectrum and residues
are retained in the JSON. No exact radical formulas or general
diagonalizability theorem are claimed.

| Width | Fixed ray | Lowest visible mass | Residue of that mass | Variance at s=0 |
|---:|---|---:|---:|---:|
| 4 | psi-minus | 3.585786438 | 1.665264893 | 1.665264893 |
| 4 | psi-plus | 6.414213562 | .0490208215 | .0490208215 |
| 8 | psi-minus | 2.819658633 | 1.923917979 | 3.295348526 |
| 8 | psi-plus | 1.955750138 | .1243186013 | .2641691831 |

Further masses in psi-minus begin `4.051285678,4.711931602,4.806883620`;
psi-plus begins `4.738819472,5.150224573,5.153880704` after its lowest mode.
The width-eight **psi-plus** projection cancels the much larger psi-minus
response and isolates the new lowest mass. At s=1 its lowest exponential
contributes .94768 of the complete psi-plus covariance. This is a descriptive
ratio of signed terms, not a probability or explained-variance estimate:
some higher spectral residues are negative, e.g. -.0578702. Positive
transition probabilities do not imply a positive spectral measure for these
readouts in a nonreversible process.

The two lowest masses have ratio **1.441727436**, below the `3/2` lower bound
of the entire width-four kappa family. A common metric rescaling or a change
of kappa cannot reproduce this pair.

For `U(s)=C0^-1 C(s)`, the linked statistic is .0320504 at s=.5,
.2013605 at s=1 and .4593400 at s=2, instead of 2. The relative Frobenius
defect `||U(1)-U(.5)^2||/||U(1)||` is .277046.
There is even a nonzero distance **s=.265657320** where the two normalized
ray correlations coincide and U(s) is scalar. This is a finite-lag crossing
of sums of ordinary modes, **not a collision of generator eigenvalues or a
Jordan point**.

## Scientific boundary and reproduction

The robust continuation is duality-protected rays plus hidden propagation
within each ray. Exact source-generated dimensions and numerical masses
describe this finite positive frontier model only. They do not count
site-Matching continuum fields or identify a thermal/Jordan family.

```sh
/Users/lc/python-envs/research-py311/bin/python scripts/p398_width8_source_spectrum.py
/Users/lc/python-envs/research-py311/bin/python -m unittest discover -s tests -p 'test_p398_width8_source_spectrum.py'
```

Two focused checks cover the old-readout identity and new finite-rank/
duality statements; no repository-wide suite or new sampling is involved.
