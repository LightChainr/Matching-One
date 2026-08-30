# Observable descriptor map audit

Status: exact registry-audit slice of Issue 146.

The canonical descriptor type has a finite valid state space. Scalar values require raw normalization
and no orientation order; orientation contrasts require one of two signed orders and may be raw or
angular-normalized. Combining these constraints with the current enums gives exactly 200 descriptors and
40,000 ordered pairs.

The executable audit applies the registered mapper to every pair. A nonzero signed angular-factor fixture
`(+2,-2)` is supplied solely to exercise raw/normalized conversion and order reversal. It is not a claim
about any physical Gaussian geometry.

Every registered affine map is checked in both directions. The reverse map must compose to the exact
identity. Every two-step composable path is also checked against the direct registered map. Unsupported
topology or quantity changes remain hard failures and are counted rather than skipped.

This establishes internal completeness and coherence of the current finite registry. It does not inspect
every prediction/scorer in the repository, add missing mathematical channel identities, or show that all
channel-bearing artifacts have adopted descriptors. Those broader Issue 146 audit tasks remain open.
