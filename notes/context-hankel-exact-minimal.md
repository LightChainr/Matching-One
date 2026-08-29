# Exact minimal context-Hankel realization for Gaussian covers

This note isolates the smallest experiment that can distinguish a two-state
unmarked endpoint realization from a third charged/marked cover state.  It is
sector-resolved and does not repeat #255's Gaussian-versus-annulus mixed
rectangle.

## Endpoint blindness is exact

Let

\[
a=1+i,\qquad d=2i.
\]

Then \(a^2=d\) as Gaussian integers.  Both paths have endpoint norm four and
Smith type \((2,2)\).  Therefore every unmarked endpoint observable that is
defined only from the final graph satisfies

\[
f(a a)=f(d).
\]

No amount of precision on such a row can reveal whether the endpoint was
presented as two degree-two covers or one degree-four cover.  Treating `aa`
and `d` as two independent endpoint observations would be pseudoreplication.

The committed exact oracle uses the constant and norm characters as a
two-dimensional product-sector witness.  On the four columns
\((\epsilon,a,aa,d)\),

\[
H_{\rm end}=\begin{pmatrix}
1&1&1&1\\
1&2&4&4
\end{pmatrix},\qquad \operatorname{rank}H_{\rm end}=2.
\]

The repeated last two columns display the blindness directly.  This table is
an exact realization witness, not a fit of the physical thermal amplitudes.

## A symmetry-allowed context adds one sector, not one scalar field

The selection rule of #244 gives a block decomposition.  Unmarked endpoint
rows live in the trivial deck sector.  An opposite-character marked matrix
element

\[
R_\chi=\frac{\langle O_{\bar\chi}S_\chi\rangle}{p(1-p)}
\]

lives in a nontrivial charged sector but is invariant as a complete
source/readout pair.  Cross-sector Hankel entries vanish exactly.  The first
enriched minor therefore has the form

\[
H_{\rm sectors}=\begin{pmatrix}
1&1&0\\
1&2&0\\
0&0&1
\end{pmatrix},\qquad \det H_{\rm sectors}=1.
\]

Its minimal dimension is exactly three: two trivial endpoint states plus one
charged state.  If the same pattern occurs in the lattice data, the third
direction is not a third scalar correction field.  It is a one-dimensional
deck/marked sector coupled to the trivial sector only by declared insertions.

An equivalent three-state witness makes the composition issue explicit:

\[
A_a=\operatorname{diag}(1,2,-1),\qquad
A_d=\operatorname{diag}(1,4,0).
\]

On the endpoint subspace, \(A_d=A_a^2\).  On the charged coordinate,

\[
A_d-A_a^2=\operatorname{diag}(0,0,-1)
\]

has rank one.  Endpoint covectors annihilate this defect.  This is an exact
identifiability counterexample: perfect product composition on every unmarked
row does not constrain the charged context.  The displayed matrices define a
model class, not a claim that the lattice defect already equals `-1`.

## Frozen risky prediction

For the declared Matching One experiment family:

```text
unmarked thermal endpoint sector: rank 2
first allowed charged deck sector: rank 1
combined unstructured table:      rank 3
```

The prediction is falsified if held-out endpoint contexts reject rank two, if
the allowed charged matrix element vanishes after exact null controls, or if
the charged block itself needs rank greater than one.

## Current archive is not enough

The repository has unmarked N65/N130/N260 endpoint material and norm-four
pilots, but it does not have one common-replica block containing the
intermediate `Z2` charged context for `aa` and the corresponding flagged
character of the direct `Z2 x Z2` cover.  Consequently no empirical Hankel
rank is reported here.  Missing entries are not completed from a latent fit.

The minimum acquisition is one N65 parent with the exact common fiber labels,
four declared contexts \((\epsilon,a,aa,d)\), two frozen endpoint outputs,
and one allowed charged source/readout pair.  Rebuild the whole partial
Hankel inside delete-one blocks.  Score endpoint rank two first, then the
charged nonzero gate, then the sector ranks `(2,1)` versus combined rank
three.  The unmarked `aa` and `d` columns are one endpoint fact, not two votes.
