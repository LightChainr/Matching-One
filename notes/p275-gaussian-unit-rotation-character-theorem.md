# Exact Gaussian unit-rotation theorem for the primitive C3 character

Status: finite exact mechanism; discovered after the held-out N145 H0 gate.

Let the square-bond torus be the Gaussian ideal quotient

\[
G_g=\mathbb Z[i]/(g),
\]

with nearest-neighbour bonds in directions `+/-1,+/-i`.  Use the period basis
`(g,ig)` and the three unoriented primitive homology lines

\[
\ell_0=(1,0),\qquad \ell_1=(0,1),\qquad \ell_2=(1,-1).
\]

## Theorem

At every finite norm, every bond probability and every Gaussian generator,

\[
P_g(\ell_0)=P_g(\ell_1).
\]

Consequently, after subtracting the common `tau=i` continuum baseline, the
complex character used by the N65/N145 gates obeys

\[
z_g=(P_0-\pi_0)+\omega^2(P_1-\pi_1)+\omega(P_2-\pi_2)
    =\omega\big[(P_2-\pi_2)-(P_0-\pi_0)\big]
    \in \omega\mathbb R,
\]

where `omega=exp(2 pi i/3)`.  Hence any two nonzero Gaussian-ideal responses
are related by a signed real gain.  Their physical relative rotation cannot
identify H4, H8 or any other local spin from this character.

## Proof

Multiplication by the Gaussian unit `i` is a well-defined automorphism of the
quotient because `i(g)=(g)`.  It maps horizontal bonds to vertical bonds and
vertical bonds to oppositely oriented horizontal bonds, so it preserves both
the finite graph and the iid bond measure.  On homology coefficients,

\[
i(mg+nig)=-ng+mig,
\]

therefore `(m,n)->(-n,m)`.  This exchanges the unoriented lines `ell_0` and
`ell_1`.  Their event probabilities are exactly equal.  The intrinsic square
torus baseline has the same equality `pi_0=pi_1`; inserting the two equalities
and `1+omega^2=-omega` gives the displayed character line.

## Relation to the two productions

The N145 observation

```text
X = 0.999926 +/- 0.007664,
p_H0 = 0.968628
```

is a Monte Carlo realization of the exact identity, not evidence for a new
spin-zero local field.  The earlier N65 two-model gate selected H8 because its
chosen angle placed the H8 line almost on the exact real-gain line while H0
was omitted from the frozen candidate set.

Likewise, a negative norm-2 cover ratio cannot by itself select H4: the exact
theorem fixes the projective character line but permits its real amplitude to
change sign with the quotient.  An H4 interpretation needs an independently
proved positive radial gain.

## Stop and scope

No further Gaussian-ideal physical-angle production with this same character
can distinguish local spin; the answer is fixed algebraically.  A future spin
experiment must leave at least one premise of the theorem: use a non-Gaussian
period lattice, a homology orbit not collapsed by the unit action, or a
different observer with an explicit map to the global channel.  This theorem
does not constrain the global `A_top`/quotient-prism H4 result.
