# Suggested re-scope of existing #636: width-four annular/lifted continuation closure

This is a proposed next brief for the existing issue, not a new duplicate issue,
production authorization, or claim that a weighted engine has already been built.
It follows the delivered width-three exact analysis. Ordinary CPU only; no GPU,
no width scan, no new Monte Carlo, no free-exponent fit.

## One decision

Can a specified local frontier state, carrying appropriate connectivity and
horizontal topological information, preserve the square-site ambient rank under
arbitrary continuation and periodic closure at circumference four? Supply a
sufficient construction with proof, or a counterexample to a precisely stated
candidate quotient. Mere small-size agreement is not the closure proof.

## Two compulsory witnesses, already solved

Bit j denotes column j. First compare open prefixes [0,5] and [7,5], followed
by the common suffix [13,0]. Their final 4x4 ranks are 0 and 1. Frontier occupancy
alone is insufficient.

The stronger pair is [13,5] and [7,5]. Both have five occupied sites, frontier
mask 5={0,2}, ordinary frontier partition {{0,2}}, and no existing transverse
winding. Their lifted 0-to-2 connecting paths have displacements -2 and +2.
Append [13,0] and close through the empty row: ranks 0 and 1, D=-1 and 0.
A state that tracks only the ordinary partition, occupation count and already
completed wrapping flag STILL fails. Both prefix weights are p^5(1-p)^3, so
weighting cannot compensate for identifying the two histories.

## Deliverable

Define the state, local row update, equivalence relation and periodic closure
functional; explicitly retain enough annular or lift information to distinguish
these histories. Prove future equivalence for every pair the representation does
identify. Explain how primitive winding, independent windings and rank saturation
are preserved without silently truncating necessary integer offsets.

Then implement only a bounded width-four reference calculation. Check the existing
4x4 Bernstein polynomial and a small non-square torus by an independent lifted
traversal. Preserve parallel periodic edges. Use integer/rational coefficients;
no asymptotic exponent or p_c claim. If a finite state bound has not been proved,
label it and do not start enumerating more widths.

If a pTL formulation from Jacobsen is reused, give the actual map between its
states/closure and the black-site rank observable; matching a root numerically
is not an intertwiner. Output a full finite probability closure before discussing
leading eigenvalues; amplitudes and all subleading terms must remain explicit.

## Stop

Stop at a demonstrated sufficient state and closure, or one precise obstruction.
Do not claim minimal state count unless a distinguishability argument proves it.
Do not re-prove the width-three formula or pay for another all-honest-torus
five-cell census. A broader engine needs a separate information-gain decision.
