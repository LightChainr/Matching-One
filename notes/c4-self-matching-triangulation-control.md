# C4 self-matching triangulation: a site-percolation parity control

Status: proposed exact control for the two-spin4/matching-parity hypothesis.

## Construction

Start from the square lattice. In every square face add exactly one diagonal, with orientation alternating by checkerboard parity of the face. Equivalently, choose each diagonal so that it connects the two vertices with even `x+y` parity.

Every original square is split into two triangles. Therefore the embedded graph is a periodic planar **triangulation**. Vertices of one checkerboard parity have degree 8 and the other parity degree 4; the pattern is invariant under a 90-degree rotation (face parity and diagonal orientation both flip).

This is a `C4` microscopic site-percolation lattice, unlike the usual uniformly oriented triangular lattice embedding.

## Why it is self-matching

For a triangulation, every face boundary is already a clique, so the site matching operation adds no edges:

\[
\widehat G=G.
\]

For an amenable periodic matching pair the critical probabilities obey the matching relation; hence a self-matching triangulation has

\[
\boxed{p_c=1/2}.
\]

This gives an exact site-percolation control in the same matching formalism as the square target.

## Finite-size central matching difference

With identical finite quotient and wrapping convention,

\[
M_L(1/2)=R_G(1/2)-R_{\widehat G}(1/2)=0
\]

identically because `G=Ghat`.

Thus the matching-odd critical-center sector is killed exactly, not merely asymptotically.

At the same time an ordinary wrapping probability on this `C4` lattice can still have orientation-dependent matching-even spin-4 finite-size corrections. This cleanly separates

- existence of a square-lattice spin-4 anisotropy;
- matching parity of that anisotropy.

## Gaussian-torus implementation

The checkerboard diagonal rule descends consistently to a Gaussian quotient when both period vectors preserve vertex parity. For periods

\[
(a,b),(-b,a),
\]

this requires `a,b` to have the same parity. Primitive pairs therefore use both odd, giving `N=a^2+b^2 = 2 mod 4`.

Useful already-validated orientation sizes include

- `N=130`: `(11,3)` vs `(9,7)`;
- `N=170`: `(13,1)` vs `(11,7)`;
- later `N=290`: `(17,1)` vs `(13,11)`;
- later `N=370`: `(19,3)` vs `(17,9)`.

For odd `a,b`, the cyclic label

\[
j=ax+by\pmod N
\]

also satisfies `j mod 2 = x+y mod 2`, so the high-degree/even checkerboard sublattice is directly identifiable in the existing `Z/NZ` representation.

Primal NN steps remain `+/-a,+/-b`. Diagonal steps `+/-(a+b), +/-(a-b)` are included according to the checkerboard/site-parity rule so that each square carries exactly one diagonal.

## Primary predictions

At `p=1/2`:

1. the matching difference `M` is exactly zero for each orientation (subject only to a correct finite matching implementation);
2. the orientation projector of the matching-even/single-lattice observable may remain nonzero;
3. if the identity-family spin-4 interpretation is correct,

\[
P_4[R](1/2)\sim A_I^{control}N^{-1}
\]

(up to `c=0` logarithmic/subleading effects);
4. there can be no analogue of the target's nonzero central matching-odd `N^-13/8` amplitude unless the observable/finite matching convention breaks the asserted self-matching identity.

The amplitude `A_I^{control}` is not predicted to equal the square/matching-pair amplitude; only exponent, spin and parity are the control targets.

## Why this is stronger than triangular-site or square-bond alone

- ordinary triangular site percolation is self-matching but microscopic `C6` symmetry suppresses a generic spin-4 lattice scalar;
- square bond percolation has exact self-duality and `C4`, but uses bond-dual rather than site-matching topology;
- this checkerboard triangulation is simultaneously **site**, **self-matching**, and **C4**.

It is therefore the closest exact control for the proposed matching-parity projector.

## Required exact tests

Before Monte Carlo:

1. verify every finite quotient face is triangular and the matching-edge generator adds nothing;
2. exhaustively verify `M(1/2)=0` on the smallest compatible quotients;
3. verify the checkerboard pattern is invariant under the intended `C4` operation up to quotient automorphism;
4. verify orientation/cyclic relabeling leaves the graph degree multiset and local face structure unchanged.

## Production test

Use fresh paired random fields at N=130 and N=170. Measure the individual wrapping orientation differences at `p=1/2`, with batch covariance. Do not spend effort estimating a central matching difference that should be algebraically zero; use it as a runtime invariant/assertion.

Then compare the angular data against a frozen `N^-1 cos(4 theta)` model and a zero-anisotropy model. Extend to N=290 only after the exact checks and first two sizes pass.
