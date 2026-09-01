# Complete leading axis capillary transfer

## Result

Put `r=m^-1` on the axis `L x L` quotient and take

```text
L,m -> infinity,        Lr -> c in [0,infinity),
Lr^2 ->0.                                                   (1)
```

The second condition is displayed because it is the bulk/interface
separation used below; it follows automatically when c is finite.  After
the dilute singleton pressure is divided out, the **complete leading
one-carrier rank-one interface partition** is

```text
G(c)=I0(2c)^2.                                               (2)
```

This includes both periodic noncrossing boundaries, all their simple NN
overhangs, narrow stripes and the alternating-face resolution.  The last
three do not create extra factors in (2): their total relative correction
is `o(1)`.  Thus the earlier two-row `cosh(c)` is a strict positive
subfamily, while `I0(2c)` is the complete leading factor per boundary.

The area-refined transfer is equally explicit.  If u is the effective
centered width field after removing the rank-zero singleton pressure, then
uniformly for c and u in compact sets,

```text
G(c,u)
 =I0(2c)^2 * sinh(u/2)/(u/2).                                (3)
```

The quotient at `u=0` is understood as one.  The Bessel factor is
independent of u at this order because an interface bridge changes the
area by `O(L)`, whereas changing the macroscopic stripe width changes it
by `O(L^2)`.

Equation (2) completes the *unmarked* capillary partition problem.  It
does **not** by itself fix the sign of original U.  In fact (3) is even in
u, so its centered first thermal derivative vanishes at the two-cloud
root.  The sign lives in the first relative correction: boundary-excluded
singleton sites, narrow-width collisions, and the island/hole insertion.
Those are precisely the terms which cancel to two extra powers at fixed L.

Keeping the endpoint interaction and **both** dilute gases produces a
sharper signed prediction.  Besides black singletons outside the carrier,
a white singleton hole is allowed in every occupied row not adjacent to a
boundary.  At the two-cloud root their pressures agree exactly.  The
continuous-time, strictly noncrossing result is therefore

```text
Phi_2gas(c)=[I0(2c)^2-I1(2c)^2]/Delta.             (S1)
```

Since `I0(x)>I1(x)>0` for finite positive x,

```text
Phi_2gas(c)>0,                                      (S2)
Ustar/A_N=-L^2m^(-(2L+1))Phi_2gas(c)[1+o(1)]<0.    (S3)
```

Thus the complete two-gas noncrossing transfer has no finite-c zero.  A
black-gas-only truncation would instead produce a spurious zero at
`c=3.9308883020798478...`; it is excluded by an actual rank-one topology,
not by a preference for the negative sign.  Formula (S3) is a theorem for
the leading two-gas noncrossing transfer.  Its promotion to the unrestricted
square-lattice source uses the overhang and connected-polymer estimates
stated below.

No sampling, finite-c fit or finite-volume enumeration is used.

## 1. The exact finite directed transfer

Fix horizontal homology and describe one boundary by its height after
each forward column.  Between successive columns a directed boundary may
make an arbitrary monotone vertical run.  The exact height kernel is

```text
T_r(a,b)=r^|a-b|,
t_r(z)=sum_(d in Z) r^|d| z^d
      =(1-r^2)/[(1-rz)(1-r/z)].                              (4)
```

Its closed partition, before imposing the finite vertical period, is

```text
B_L(r)=[z^0] t_r(z)^L.                                      (5)
```

At `r=c/L`,

```text
L log t_(c/L)(exp(i theta))=2c cos(theta)+O_C(L^-1),          (6)
```

and therefore

```text
B_L(c/L) -> I0(2c)
 =sum_(j>=0)c^(2j)/(j!)^2.                                  (7)
```

The combinatorics in (7) are transparent.  A leading bridge has j upward
and j downward unit steps at asymptotically independent column positions.
The factor `1/(j!)^2` counts the two unordered position sets.  The old
two-row family forces the signs to alternate and retains only
`sum c^(2j)/(2j)!=cosh(c)`.

Runs of vertical length at least two have aggregate activity
`O(Lr^2)=o(1)`.  Keeping them in (4) is useful because (5) is an exact
directed transfer before taking the limit.

## 2. Horizontal reversals vanish in the complete simple-curve sum

It remains to show that (7) is not only a directed submodel.  Let an
oriented NN curve in homology `(L,0)` have E,W,U,D steps.  Closure gives

```text
E-W=L,      U=D=n,
length=L+2(w+n),       w=W.                                  (8)
```

The source weight relative to a straight boundary is `r^(2(w+n))`.
Forget simplicity and cyclic identifications; a based step word with
fixed `(w,n)` is then bounded by

```text
C_L(w,n)
 =(L+2w+2n)!/[(L+w)! w! n! n!].                             (9)
```

For `r=c/L`, every west step removes one power of positional entropy:

```text
C_L(w,n) r^(2w+2n)
 =O_C[L^-w c^(2w+2n)/(w! n!^2)]                             (10)
```

after summation on any fixed compact c interval.  Splitting
`w+n<=L/4` and its factorial tail makes (10) uniform; the tail is smaller
than every fixed inverse power of L.  Hence

```text
sum_(w>=1,n>=0) C_L(w,n)r^(2w+2n)=O_C(L^-1).                 (11)
```

Self-avoiding curves are a subset of these relaxed words.  Equation (11)
therefore proves that every curve with a horizontal reversal is negligible
relative to (7).  It also explains why a finite overhang cluster does not
survive at fixed c: its compensating forward step costs two edges but does
not supply a second independent column position.

Vertical wrapping and self-contact inside the directed class require
total vertical variation at least L.  The same factorial tail makes their
weight `o(1)`.  Thus (7) is the complete one-boundary leading partition,
not only a transfer ansatz.

## 3. Two noncrossing periodic boundaries factor at leading order

Choose the two boundary heights at one reference column, with cyclic gap
`d in {1,...,L-1}`.  At zero excess this is exactly the usual stripe-width
zero mode.  Conditional on total vertical variation S of the two bridges,
the curves can collide only when

```text
d<=S   or   L-d<=S.                                          (12)
```

There are at most `2S` such starting gaps.  Under the weights (5), all
moments of S are bounded uniformly for c in a compact interval.  Averaging
over the `L-1` gaps gives

```text
Pr_weighted(collision)=O_C(L^-1).                            (13)
```

The same estimate covers a shared unsmoothed face centre.  The local
occupied-corner resolution either keeps the two curves separated or puts
the configuration in the collision error.  Narrow stripes are exactly
the `O(S)` exceptional gaps in (12), not a separate extensive family.

Combining (7), (11), and (13), the normalized exact annulus sum for one
fixed orientation satisfies

```text
1/[L(L-1)] sum_(side-decorated simple annuli A)
 r^[|partial_1 A|+|partial_2 A|-2L]
   -> I0(2c)^2.                                             (14)
```

The two lattice orientations and the common translation factor are the
same as for rigid stripes and cancel in this normalization.  This proves
(2).

Multiple essential occupied components are not part of the same leading
interface: the component winding bound charges another `2L-2`, which is
smaller than any capillary entropy factor after multiplication by
`r^(2L-2)`.  Contractible nonsingleton components have total relative
activity `O(L^2r^4)=o(1)`.  Contractible singletons are not discarded;
they form the bulk pressure treated in Section 6.

## 4. Area-refined transfer

For an annulus A let `K_A` be its occupied area.  At the reference column
write its width as w.  If S is the total vertical variation, then

```text
K_A=Lw+delta K_A,       |delta K_A|<=LS.                      (15)
```

Use the centered area field `exp[u(K_A-N/2)/N]`.  Since S has bounded
exponential moments on compact c sets,

```text
exp[u delta K_A/N]=1+O_(C,U)(S/L).                            (16)
```

The interface factor is consequently independent of u at leading order.
The width average is a Riemann sum:

```text
1/(L-1) sum_(w=1)^(L-1) exp[u(w/L-1/2)]
 -> integral_0^1 exp[u(x-1/2)]dx
  =sinh(u/2)/(u/2).                                         (17)
```

Equations (14), (16), and (17) prove (3).  In particular,

```text
partial_u G(c,u)|_(u=0)=0,
partial_u^2 log G(c,u)|_(u=0)=1/12.                          (18)
```

The `1/12` is the macroscopic stripe-width variance.  The bridge-area
variance derived in the directed calculation is a genuine next-order
`O(L^-2)` correction to `(K/N)` and is not a leading first derivative.

## 5. The endpoint hard-core kernel and the two-gas signed closure

The `O(1/L)` narrow-width term in the unmarked partition cannot be thrown
away in original U.  Its fraction is `O(1/L)`, but its centered occupation
leverage is `O(L^2)`.  The continuous-time unit-jump limit makes this term
explicit.

The one-particle height kernel is

```text
p_c(a,b)=I_(b-a)(2c).                                       (19)
```

Start two ordered boundaries at heights `0<d` and require the same
endpoints after one horizontal period.  Karlin--McGregor gives the strict
nonintersection partition

```text
J_d(c)=det [[I0(2c), I_d(2c)],
            [I_d(2c), I0(2c)]]
      =I0(2c)^2-I_d(2c)^2.                                  (20)
```

For a bulk width, `d->infinity` and `J_d->I0^2`.  The two endpoint widths
`d=1,L-1` instead carry the survival ratio

```text
rho(c)=J_1(c)/I0(2c)^2
      =1-[I1(2c)/I0(2c)]^2.                                 (21)
```

This is the exact endpoint correction in the named continuous-time
noncrossing model.  It is not obtained by multiplying every width by the
bulk `I0^2` factor.

The second dilute gas is forced by the actual topology.  Let a straight
rank-one carrier have width w.  A vacant site in one of its `w-2` interior
rows has four occupied NN neighbours.  Removing it adds four mixed bonds
but changes neither the occupied component count nor its ambient rank:

```text
Delta Bmix=4,       Delta C_B=Delta r=0,
Delta g=4.                                                   (22)
```

It is therefore a legal isolated white hole of exact relative activity
`a^2/h=m^-4h^-1`.  The matching convention causes no alias: the excluded
two boundary rows keep the hole outside the matching-white exterior, and
the identity `g=Bmix-2C_B+r` already includes the white-component change.
Adjacent or matching-connected hole pairs begin at total activity
`O(N/m^6)=o(N^-1)` here, so they cannot alter the marked O(1) limit.

For a stripe of width w define the available interior counts

```text
H_w=L max(w-2,0),
M_w=L max(L-w-2,0).                                         (23)
```

The complete singleton factor is

```text
W_w(h)=h^(Lw)(1+a^2/h)^H_w(1+ah)^M_w.                       (24)
```

For bulk widths `2<=w<=L-2`, put

```text
C=h+a^2,       A=1+ah.
```

Then

```text
W_w(h)=h^(2L) C^[L(w-2)] A^[L(L-w-2)].                      (25)
```

At the ideal two-cloud root `h=1+a`,

```text
C=A=1+a+a^2.                                                 (26)
```

Thus the width tilt vanishes exactly, not merely asymptotically.  Moving
a boundary trades a C-site for an A-site, so the bridge-area tilt vanishes
at the same time.  The earlier black-gas-only terms

```text
-c^4/12 + c^2[c I1(2c)/(3I0(2c))]                            (27)
```

are precisely canceled by the white-hole gas.  Omitting it would give the
spurious bracket

```text
1-[I1/I0]^2-c^2/6+(2c/3)I1/I0,                              (28)
```

whose numerical zero is `3.9308883020798478...`.

Differentiate (24) and subtract the separately normalized two-phase mean.
For `2<=w<=L-2`, the A and C terms form a discrete difference in w; pairing
w with `L-w` cancels the bulk sum at (26).  The max functions in (23) stop
that difference at `w=1` and `w=L-1`.  This is the transfer version of the
fixed-L island-versus-hole cancellation: the surviving mark is an endpoint
insertion, not a uniform mark on all widths.  Equation (20) changes its
rigid transfer factor from one to `J_1(c)`.  Hence

```text
Phi_2gas(c)=J_1(c)/Delta
 =[I0(2c)^2-I1(2c)^2]/Delta,                                (29)
```

which proves (S1)-(S3) inside the two-gas noncrossing transfer.

Positivity is exact: the integral representation makes
`0<I1(x)<I0(x)` for every finite `x>0`.  At large c,

```text
J_1(c)=exp(4c)/(8 pi c^2)[1+O(c^-1)],                        (30)
```

so capillary roughening enhances the magnitude but cannot reverse the
sign in this model.

The remaining full-lattice boundary is now narrower.  A horizontal
overhang has `o(1)` total mass in (2); its relative `O(1/L)` correction to
the endpoint kernel multiplies the already scaled endpoint insertion and
is `o(1)` in (29).  Corner-dependent losses change H and M by `O(S)` for
S capillary jumps, hence their log weight is `O(aS)=O(L^-2)` and also
vanishes after the O(1) marked scaling.  Adjacent black sites which form a
forest already have product activity and are included in A.  The first
black correction is a cycle, of aggregate size `O(N/m^6)`; the first
interacting white-hole pair has the same bound.  Their marked effects are
`O(N^2/m^6)=o(1)` under `L/m->c`.  These estimates promote (29) to the
leading one-carrier square-lattice interface law, subject only
to the already proved suppression of extra essential components.

## 6. Restoring both dilute gases and the common root

The black exterior gas alone gives the misleading effective field
`R=h/(1+ah)`.  The correct carrier has the two pressures A and C in
(25).  Define instead

```text
u=N log(C/A).                                                (31)
```

The ideal two-cloud root is `h=1+a`, where u is exactly zero by (26).  The
connected-pressure correction moves this root by `O(a^3)`, so

```text
u_root=O(Na^3)=O(L^-4),
N u_root=O(N^2a^3)=O(L^-2).                                 (32)
```

The second line is the actual centered-occupation leverage.  It vanishes.
Therefore neither the macroscopic width variance nor the bridge-area
variance contributes to the leading signed answer once both physical gases
are retained.

For bulk widths the complete uncentered one-carrier generating function is

```text
Z1,axis
 =2L(L-1)r^(2L-1)h^(2L)C^(-2L)A^(N-2L)
   I0(2c)^2 (e^u-1)/u [1+o(1)].                             (33)
```

The continuous value at u=0 is one.  The two endpoint widths replace the
bulk Bessel factor by `J_1`, and their marked difference yields (29).
Equations (29), (32), and the positive two-phase denominator give the
original-U law (S3).  The tilted companion is exponentially later by the
independent shortest-winding theorem, so it cannot alter this sign.

This also pinpoints the error in a black-only capillary reduction.  It
retains the u displacement and bridge-area variance while dropping the
white holes that set `C=A`; the resulting zero is a truncation artifact.
The cancellation is microscopic and exact, not an appeal to continuum
symmetry.

## 7. The full-lattice error scale

The signed endpoint term is O(1) after the fixed
`L^2m^(-(2L+1))` normalization.  Every omitted class has a smaller marked
scale:

```text
horizontal reversal in endpoint kernel     O(L^-1),
corner correction to H_w,M_w               O(m^-2),
black cycle / interacting white holes      O(N^2/m^6),
second essential occupied component        O(m^(-2L+2)).     (34)
```

All four vanish when `L/m->c<infinity`.  The first follows from the lost
positional-entropy bound (11), now applied with the endpoint insertion.
The second uses the tight capillary jump count.  The third is
`O(c^6/L^2)`, and the fourth beats every Bessel entropy factor.  These are
the corrections that had to be checked after the unmarked theorem; simply
calling the endpoint fraction `O(1/L)` would not have been enough because
its occupation leverage is macroscopic.

Consequently the leading one-carrier result (S3) is stable in the full
axis rank-one sector under the existing two-cloud contour gate.  It proves
that original U remains negative for every fixed finite c.  It does not
address `c->infinity` jointly with L, where the uniform constants in the
factorial and endpoint bounds must be reopened.

## Scientific card

- **Complete theorem:** every leading simple axis annulus, including
  overhangs and noncrossing interaction, has normalized partition
  `I0(2c)^2`; the full area transform is (3).
- **Mechanism changed:** `cosh(c)` is replaced by a two-boundary Bessel
  transfer, while the dilute pressure restores a centered width zero mode.
- **Signed mechanism:** bulk widths cancel after both black exterior and
  white interior gases are included; the surviving mark is the hard-core
  endpoint kernel `I0(2c)^2-I1(2c)^2`.
- **Sign:** under the existing two-cloud contour gate, original U remains
  negative for every fixed finite c.  The black-only zero at
  `3.930888...` is an omitted-white-hole artifact.
- **Remaining boundary:** a regime with `c` itself diverging; the present
  estimates are uniform only on compact c intervals.
- **Scope:** finite-lattice strong-source joint asymptotics.  No continuum
  interface field, fixed-m limit or all-rank-one finite-c enumeration is
  asserted.
