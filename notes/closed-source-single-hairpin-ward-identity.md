# A single hairpin is a positive Dirichlet boundary current

## Result

The natural telescoping Ward identity does **not** imply
`kappa_hp(infinity)=0`.  Complement pairing turns the two orientations
of a rooted single hairpin into a discrete divergence.  On the full
height line that divergence vanishes, but the noncrossing relative-gap
transfer lives on the Dirichlet half-line `d=1,2,...`.  The missing
`d=0` state leaves a positive boundary current.

After factoring the microscopic positive hairpin fugacity and its
longitudinal embedding count, the normalized coefficient is

\[
 \boxed{\kappa_{hp}(c)=
 2\,{I_1(2c)-I_3(2c)\over I_0(2c)-I_2(2c)}
 =4\,{I_2(2c)\over I_1(2c)}>0,}                   \tag{1}
\]

and therefore

\[
                         \boxed{\kappa_{hp}(\infty)=4.}            \tag{2}
\]

Thus a single endpoint-sensitive hairpin is not killed by signed
complement pairing.  Its sign is positive in the finite-gap transfer.
The zero Ward identity survives only after tracing over all separations,
or for a hairpin packet that never sees the noncrossing boundary.

This resolves one ambiguity left by
[`closed-source-axis-growing-capillary-window.md`](closed-source-axis-growing-capillary-window.md):
at `alpha=L^2/m^3=O(1)`, the width-one hairpin can carry a genuine
first-order coefficient.  Equation (2) does not supply its microscopic
fugacity; that finite digital-packet count remains a separate local
input.

## 1. Relative-gap Dirichlet transfer

Let `S_+|d>=|d+1>` on `ell^2(N)` and let

\[
 S_-|1\rangle=0,\qquad S_-|d\rangle=|d-1\rangle\quad(d\ge2).
\]

The directed relative-gap generator is

\[
 H_D=S_++S_- .                                     \tag{3}
\]

Its heat kernel is the image kernel

\[
 K_c(d,e)=\langle d|e^{cH_D}|e\rangle
 =I_{d-e}(2c)-I_{d+e}(2c).                         \tag{4}
\]

The subtraction is exactly the noncrossing Dirichlet condition at
`d=0`.  The centre-of-mass bridge factor is common to the unperturbed
and hairpin terms and cancels from their ratio, so (4) is the minimal
finite-gap matrix required here.

## 2. Complement pairing is a commutator

Erase a rooted horizontal reversal from the lower boundary.  Sliding
its attachment one directed step to the right gives the corresponding
upper/complement packet with the opposite signed occupation mark.  In
relative-gap coordinates the two attachments are represented by the
antisymmetric current

\[
                         B=S_+-S_- .                \tag{5}
\]

Their signed difference is the discrete divergence

\[
                         V_{hp}=[H_D,B].            \tag{6}
\]

On the full line the two shifts commute and (6) is zero.  On the
Dirichlet half-line,

\[
 S_-S_+=1,\qquad S_+S_-=1-P_1,qquad
 P_1=|1\rangle\langle1|,
\]

so the exact Ward identity is instead

\[
                         \boxed{V_{hp}=2P_1.}       \tag{7}
\]

This is the minimal obstruction to `kappa_hp=0`.  It is not a failure
of complement pairing: complement pairing is what makes the bulk a
commutator.  Noncrossing converts that commutator into boundary flux.

The two-cloud root is needed here.  At `h=1+m^-2`, black exterior and
white interior packet weights are equal after complement, so no area
tilt is left in (5).  Without that equality an additional diagonal
collar operator would accompany (7), and one could not call its
coefficient the pure single-hairpin response.

## 3. Telescoping leaves a positive endpoint term

Insert (7) at every possible longitudinal time.  Duhamel's formula and
(6) give

\[
 \begin{aligned}
 \mathcal H(c)
 &=\int_0^c\langle1|e^{(c-s)H_D}V_{hp}e^{sH_D}|1\rangle ds\\
 &=\langle1|e^{cH_D}B-Be^{cH_D}|1\rangle\\
 &=2K_c(1,2).                                      \tag{8}
 \end{aligned}
\]

Equivalently, using (7),

\[
 \mathcal H(c)=2\int_0^cK_{c-s}(1,1)K_s(1,1)ds.   \tag{9}
\]

Both forms make the sign strict: every term in (9) is positive.  Divide
by the unperturbed endpoint kernel `K_c(1,1)` to obtain

\[
 \kappa_{hp}(c)=2{K_c(1,2)\over K_c(1,1)}.         \tag{10}
\]

Substituting (4) and using

\[
 I_0(2c)-I_2(2c)={I_1(2c)\over c},\qquad
 I_1(2c)-I_3(2c)={2I_2(2c)\over c}
\]

proves (1).  Since `I_2(2c)/I_1(2c)->1`, equation (2) follows.

## 4. What the zero trace does and does not say

For a finite cutoff `1<=d<=D` with Dirichlet boundaries at both ends,

\[
                         \operatorname {tr}[H_D,B]=0.              \tag{11}
\]

The lower current `+2P_1` is canceled by the upper cutoff current.
Sending `D` to infinity before taking an endpoint matrix element removes
the latter; it does not cancel (7).  Therefore summing the hairpin over
all stripe widths can telescope to zero while the width-one thermal
response remains positive.  The original global-U mechanism reads the
endpoint, not the unweighted trace, so (11) is not the relevant null.

Similarly, a one-boundary overhang far from the other interface is a
translation-invariant Toeplitz renormalization and carries no `P_1`
term.  Equation (2) concerns only the gap-sensitive hairpin whose erased
packet reaches the width-one Dirichlet wall.

## Scientific boundary

The value `4` uses the unit normalization in (6): one complement-paired
rooted hairpin corresponds to one discrete current difference.  If the
microscopic digital packet has fugacity `w_hp`, orientation multiplicity
or an additional collar state, its contribution is `4w_hp` times that
multiplicity, or the analogous matrix element in the enlarged finite-gap
kernel.  Signed pairing alone fixes the positive Dirichlet factor; it
does not determine those microscopic prefactors.

No claim is made about a continuum field, fixed-`m` limit, or the total
fixed-`alpha` correction after other endpoint packets are added.
