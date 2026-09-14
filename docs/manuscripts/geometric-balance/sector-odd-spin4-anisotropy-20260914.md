# Sector-odd rotational anisotropy as the source of the 4,6,... pseudo-critical exponents

2026-09-14.  CONJECTURAL CFT/lattice-symmetry interpretation, replacing the earlier scalar Kac-(4,2) guess.

The starting facts are unusually restrictive:

1. square-site and kagome-bond percolation are expected to have the same continuum percolation CFT;
2. Jacobsen's semi-infinite-cylinder eigenvalue criterion has leading pseudo-critical exponent `4` on the square lattice;
3. on the kagome lattice the amplitude of that `4` correction vanishes and the leading exponent is `6`;
4. Jacobsen explicitly suggests the change from four-fold to three-/six-fold lattice rotation as the reason for the missing kagome amplitude;
5. on the square-site safe transfer, the critical open/closed sector mismatch scales consistently with `w^-17/4`, while the thermal derivative scales as `w^-1/4`.

A scalar correction does not naturally explain item 3.  A spinful lattice-anisotropy correction does.

## 1. Thermal family and sector parity

The percolation thermal field has

\[
x_t=5/4,
\qquad h_t=\bar h_t=5/8.                                      \tag{1.1}
\]

Under primal/dual exchange the thermal perturbation changes sign: moving one side of criticality for the primal model corresponds to the opposite side for the matching/complementary description.  It is therefore natural for the **sector-odd** correction spectrum to live in the thermal conformal family rather than the identity family.

The relevant lattice observable is

\[
\Theta_w(p_c)
=I^0_{4,w}(p_c)-I^0_{8,w}(1-p_c),                             \tag{1.2}
\]

the difference of the two magnetic/topological excitation energies.  Sector-even irrelevant fields can be large in each term and cancel from (1.2).

## 2. Square lattice: first allowed anisotropic thermal descendant

A level-`s` chiral descendant of the thermal primary has weights

\[
(h_t+s,\bar h_t)
\quad\hbox{or}\quad
(h_t,\bar h_t+s),                                             \tag{2.1}
\]

and therefore

\[
x=x_t+s,
\qquad
\text{spin}=\pm s.                                           \tag{2.2}
\]

A reflection-invariant lattice perturbation uses the real combination of the two spins.

On a square lattice, `C4` invariance permits spin `s` only when

\[
s\equiv0\pmod4                                               \tag{2.3}
\]

among nontrivial rotational anisotropies.  The first possibility is therefore `s=4`:

\[
\boxed{
x_{odd,\square}=x_t+4=21/4,
\qquad \text{spin}=4.}                                       \tag{2.4}
\]

Its contribution to a **per-row** cylinder excitation energy scales as

\[
\Theta_w(p_c)\sim A_4 w^{1-x_{odd}}
=A_4 w^{-17/4}.                                               \tag{2.5}
\]

The thermal response scales as

\[
\Theta'_w(p_c)\sim B w^{-1/4}.                               \tag{2.6}
\]

Hence the zero of the sector difference moves by

\[
\boxed{
p_w^{ch}-p_c
\sim-\frac{A_4}{B}w^{-4}.}                                   \tag{2.7}
\]

This reproduces the square-site pseudo-critical exponent four without requiring a scalar field of dimension `21/4`.

## 3. Kagome / triangular rotational selection

The kagome and triangular geometries have six-fold spatial symmetry in the bulk; Jacobsen phrases the relevant change relative to the square basis as three-fold rotational symmetry.  In either description, a spin-four perturbation is not invariant:

\[
e^{i4(2\pi/3)}\ne1,
\qquad
e^{i4(\pi/3)}\ne1.                                        \tag{3.1}
\]

So the square spin-four amplitude must vanish on a symmetry-preserving kagome/triangular realization:

\[
\boxed{A_4^{\rm kagome}=0.}                                  \tag{3.2}
\]

With inversion/reflection excluding odd-spin anisotropies, the first allowed thermal-family anisotropy is spin six:

\[
\boxed{
x_{odd,\hex}=x_t+6=29/4,
\qquad \text{spin}=6.}                                       \tag{3.3}
\]

It gives

\[
\Theta_w(p_c)\sim A_6 w^{1-29/4}=A_6 w^{-25/4},              \tag{3.4}
\]

and after division by the same thermal response,

\[
\boxed{p_w-p_c\sim w^{-6}.}                                  \tag{3.5}
\]

This is exactly the leading exponent Jacobsen finds for kagome bond percolation after setting the square-like `A_1` term to zero.

## 4. Why the ordinary square spin-four identity-family field does not already give the root shift

Square critical lattice models generically possess lower-dimensional spin-four / identity-family anisotropy corrections.  A familiar candidate has total dimension around four and produces ordinary corrections in each finite-size energy.

The safe-transfer data indeed show a much larger common correction in the average magnetic excitation energy, consistent with a per-row `w^-3` term (`x=4`).

But the matching root is controlled by the **difference**

\[
I^0_4-I^0_8.                                                   \tag{4.1}
\]

The two sectors represent the same continuum magnetic primary.  A sector-even identity-family anisotropy can therefore have equal amplitudes and disappear from (4.1).  What survives must be both

1. allowed by lattice rotation;
2. odd under primal/dual topological-sector exchange.

A spin-four descendant of the thermal family satisfies exactly that structural request.

## 5. A more precise operator statement to test

The conjecture should not be phrased as the literal operator `partial^4 epsilon`.  Total derivatives can integrate to zero, and Virasoro descendants mix.  The correct target is:

> **Sector-odd spin-four conjecture.**  The lowest-dimension `C4`-invariant, reflection-even, primal/dual-odd irrelevant scaling field with a nonzero matrix-element difference between the two magnetic cylinder sectors is a spin `±4` quasi-primary/descendant in the thermal Virasoro module, of total dimension `x_t+4=21/4`.

For six-fold lattices the corresponding first allowed field is the spin `±6` thermal-module descendant of dimension `x_t+6=29/4`.

This formulation leaves room for logarithmic mixing inside the `c=0` theory while making the symmetry and dimension claims falsifiable.

## 6. Immediate falsification tests

### 6.1 Square orientation test

A genuine spin-four amplitude changes phase/sign under rotation of a weak anisotropic perturbation according to `e^{i4 theta}`.  Introduce a controlled rectangular/diagonal anisotropy while keeping the continuum thermal coordinate fixed and test the angular harmonic of the sector mismatch.

### 6.2 Six-fold lattice test

For a C6-symmetric percolation realization, the direct sector mismatch at criticality should have no `w^-17/4` term.  Its first rotational thermal-family contribution should instead scale as `w^-25/4`, yielding a pseudo-critical `w^-6` shift.

Jacobsen's kagome data already supply strong qualitative support for this selection rule; a direct fit of the **sector free-energy difference**, rather than only the root, would be the sharper test.

### 6.3 Symmetry breaking on kagome/triangular

Add a small perturbation that reduces C6/C3 to C2 or C1 while preserving critical tuning.  If the missing exponent four is a spin-four anisotropy channel, an amplitude linear in the symmetry-breaking coupling should reappear.

### 6.4 Self-matching triangular site control

Critical triangular-site percolation is exactly self-matching at `p=1/2`.  In a period convention that preserves the matching identification, the primal/dual charge-sector difference should vanish much more strongly than on the non-self-matching square lattice.  This provides a clean control separating duality-odd amplitudes from rotationally allowed but sector-even corrections.

## 7. Relation to Jacobsen's correction ladder

On the square lattice, further spin-four thermal-family descendants and/or higher derivative levels can generate exponents `6,8,...` after the leading spin-four contribution.  On a six-fold lattice the spin-four amplitude is removed, making exponent six the first visible term.

The precise origin of every later even power is not fixed by this note.  The robust claim is narrower:

\[
\boxed{\text{rotation selection naturally distinguishes square }4
\text{ from kagome }6.}                                      \tag{7.1}
\]

That distinction is not explained by a scalar `(4,2)` assignment.

## 8. Literature boundary

Jacobsen, arXiv:1507.03027, reports `Delta_1=4` for square-site percolation, `Delta_2=6` as the leading kagome-bond correction, and explicitly proposes rotational symmetry as a possible reason the kagome `A_1` amplitude vanishes.  He does not identify a specific irrelevant CFT field.

The general symmetry principle that square-lattice anisotropy permits spin multiples of four whereas triangular/hexagonal symmetry permits multiples of six is standard lattice-CFT finite-size-scaling logic.  The specific identification with a **primal/dual-odd thermal-family spin-four descendant of dimension `21/4`** is a new conjecture here and requires operator/sector verification.

## 9. Claim boundary

The symmetry selection `C4` allows spin four while `C3/C6` forbids it is exact group theory.  The observed square/kagome correction exponents are literature facts.  The assignment of the sector-odd mismatch to thermal-family spin `4` / `6` descendants is conjectural, albeit more structurally compatible with the lattice comparison than the withdrawn scalar Kac-(4,2) guess.
