# A single hairpin has positive flux but vanishing physical coefficient

## Corrected result

Complement pairing turns the two orientations of a rooted single
hairpin into a discrete divergence.  On the full height line that
divergence vanishes.  On the noncrossing Dirichlet half-line it leaves
a positive boundary current.

That current is **not** the coefficient `kappa_hp` defined in
`1163bebe`.  The stripped physical kernel contains the longitudinal
placement average `(1/L) sum_j`.  In capillary time `s=j/m`,

\[
 {1\over L}\sum_j\longrightarrow {1\over c}\int ds,
 \qquad c={L\over m}.                              \tag{1}
\]

Let

\[
 \rho(c)=1-\left({I_1(2c)\over I_0(2c)}\right)^2.
\]

The unaveraged endpoint-relative Ward flux is

\[
 \eta_{flux}(c)=4{I_2(2c)\over I_1(2c)}>0.         \tag{2}
\]

After the physical placement average, and then the bulk normalization
by `I_0(2c)^2` required in `1163bebe`, the exact coefficient is

\[
 \boxed{\kappa_{hp}(c)=
 {4\over c}{I_2(2c)\over I_1(2c)}
 \left[1-\left({I_1(2c)\over I_0(2c)}\right)^2\right]>0.}         \tag{3}
\]

Consequently

\[
 \boxed{\kappa_{hp}(\infty)=0,\qquad
        c^2\kappa_{hp}(c)\longrightarrow2.}         \tag{4}
\]

Thus the positive Dirichlet flux is compatible with the physical
single-hairpin null.  The `alpha=c beta` crossover term proposed in
`1163bebe` receives no nonzero limiting constant from this packet.
The value `4` belongs to the unaveraged endpoint-relative flux, not to
the repository's bulk-normalized `kappa_hp`.

## 1. Relative-gap Dirichlet transfer

Let `S_+|d>=|d+1>` on `ell^2(N)` and let

\[
 S_-|1\rangle=0,\qquad S_-|d\rangle=|d-1\rangle\quad(d\ge2).
\]

The directed relative-gap generator is

\[
 H_D=S_++S_- .                                     \tag{5}
\]

Its heat kernel is the image kernel

\[
 K_c(d,e)=\langle d|e^{cH_D}|e\rangle
 =I_{d-e}(2c)-I_{d+e}(2c).                         \tag{6}
\]

The subtraction is exactly the noncrossing Dirichlet condition at
`d=0`.  The centre-of-mass bridge factor is common to the unperturbed
and hairpin terms and cancels from endpoint-relative ratios.

## 2. Complement pairing is a boundary commutator

Erase a rooted horizontal reversal from one boundary.  Sliding its
attachment by one directed step gives the complement packet with the
opposite signed occupation mark.  In relative-gap coordinates their
current is

\[
                         B=S_+-S_- .                \tag{7}
\]

The signed difference is the discrete divergence

\[
                         V_{hp}=[H_D,B].            \tag{8}
\]

On the full line the shifts commute.  On the Dirichlet half-line,

\[
 S_-S_+=1,\qquad S_+S_-=1-P_1,qquad
 P_1=|1\rangle\langle1|,
\]

and hence

\[
                         \boxed{V_{hp}=2P_1.}       \tag{9}
\]

Complement pairing therefore kills the bulk hairpin but localizes a
strictly positive current at the noncrossing wall.  The two-cloud root
is needed: at `h=1+m^-2`, black exterior and white interior packet
weights are equal, so no separate area-tilt operator accompanies (9).

## 3. Ward telescoping before placement normalization

Insert (9) over capillary time.  Duhamel's formula gives

\[
 \begin{aligned}
 \mathcal H_{raw}(c)
 &=\int_0^c\langle1|e^{(c-s)H_D}V_{hp}e^{sH_D}|1\rangle ds\\
 &=\langle1|e^{cH_D}B-Be^{cH_D}|1\rangle\\
 &=2K_c(1,2).                                      \tag{10}
 \end{aligned}
\]

Equivalently,

\[
 \mathcal H_{raw}(c)=
 2\int_0^cK_{c-s}(1,1)K_s(1,1)ds,                 \tag{11}
\]

which makes positivity explicit.  Relative to the endpoint kernel,

\[
 \eta_{flux}(c)=2{K_c(1,2)\over K_c(1,1)}.         \tag{12}
\]

Substitute (6) and use

\[
 I_0(2c)-I_2(2c)={I_1(2c)\over c},\qquad
 I_1(2c)-I_3(2c)={2I_2(2c)\over c}
\]

to obtain (2).  This is the step at which an unnormalized calculation
finds the nonzero limit `4`.

## 4. Physical placement and bulk normalization

The repository definition is

\[
 \widehat H_c=\beta^{-1}H_c^{phys},\qquad
 \beta={L\over m^2},
\]

with the rooted packet averaged over `L` longitudinal positions.  One
column is capillary time `1/m`, so (1) changes (12) to

\[
 {\langle1|\widehat H_c|1\rangle\over K_c(1,1)}
 ={1\over c}\eta_{flux}(c).                       \tag{13}
\]

The endpoint determinant satisfies

\[
 {Z_{endpoint}\over I_0(2c)^2}=\rho(c).
\]

Therefore

\[
 I_0(2c)^{-2}\langle1|\widehat H_c|1\rangle
 =\rho(c){\eta_{flux}(c)\over c},                 \tag{14}
\]

which is exactly (3).  Since

\[
 \rho(c)={1\over2c}+O(c^{-2}),\qquad
 {I_2(2c)\over I_1(2c)}=1+O(c^{-1}),
\]

equation (4) follows.  Restoring the physical fugacity, the actual
endpoint-relative perturbation is

\[
 \beta{\eta_{flux}(c)\over c}
 ={4\beta\over c}\{1+O(c^{-1})\},                 \tag{15}
\]

not `4 beta`.

## 5. What the trace Ward identity says

For a finite separation cutoff with Dirichlet boundaries at both ends,
`tr[H_D,B]=0`: the lower current `+2P_1` is canceled by the upper cutoff
current.  Sending the upper boundary to infinity leaves the positive
lower flux (9).  This explains why the raw endpoint current is nonzero,
but it does not supply the physical `1/c` placement normalization.

A one-boundary overhang far from the other interface is a common
Toeplitz renormalization and carries no `P_1` term.  Equations (3)--(4)
concern only the gap-sensitive, complement-paired, rooted single
hairpin.

## Scientific boundary

Any fixed orientation multiplicity multiplies (3) but still leaves
`kappa_hp(infinity)=0`.  An enlarged collar state can evade this result
only if its stripped placement sum grows faster than the single-current
Duhamel term.  This note does not determine such a multi-state packet,
the total fixed-alpha correction, a continuum field, or a fixed-m
limit.
