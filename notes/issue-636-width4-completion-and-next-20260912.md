## Width-four sufficiency task completed — proposed update to this existing issue

The latest #636 comment's requested state is now constructed and checked.
Do not commission another width-four representation or census for that question.
The additive submission bundle contains the proofs, code and complete finite table;
it is not yet committed by this message.

For the fixed width-four NN square-site row language, retain pinned seam and
current boundary, component partition, path gains modulo the accumulated GLOBAL
cycle span, and the latter span even after its component leaves the boundary.
Replacing old interiors by gained stars preserves all future cycle spans.
A common integer shear of current-to-pinned path gains changes final windings
by (hx,hy)->(hx+c*hy,hy) and preserves rank, but NOT directional/primitive labels.

Complete closure: 1448 representatives, all 23168 row transitions independently
reconstructed. Deterministic continuation refinement: 3 -> 164 -> 509 -> 509.
Every class is reachable at a physical length >=2; every pair is distinguished
by some suffix of at most two rows. This is a rank-language minimum, not a
linear state count or continuum claim.

The explicit normalized operator is M_(4,m)(p)=b(p)^T K(p)^(m-1)c for all m>=2,
where c=r-1 and every row mask receives its original Bernoulli weight.
At p=1/2 the scalar integer sequence 16^m M has minimal recurrence order 15,
with six trace factors and weights (-2,1,-1,2,1,1). The identity has a 509-
moment exact operator certificate and a nonzero physical-tail 15x15 Hankel
determinant; it is not a finite-sequence fit promoted to a theorem.

Executed controls: 74640 full physical-graph configurations, every occupation
coefficient, all finite transitions, 7 local tests. No Monte Carlo/GPU. Full
repository CI for the NEW patch has not run. #707's own CI 34689086221 is now
successful and is a different check.

### One remaining bounded task worth doing here

From the supplied K(p), compute/prove the parameter-dependent observable
reduction over Q(p), separately for P0, P2 and P2-P0. Give the actual visible
spectral denominator and closure weights; specialize exactly to the p=1/2
certificate. A fraction-free Krylov relation or a rigorously degree-bounded
interpolation followed by exact identity verification is acceptable. A handful
of numerical p specializations alone is not.

Question to decide: do the two leading observable weights stay equal, and
which other visible modes survive the closure? This is not supplied by the
stochastic matrix's Perron eigenvalue 1 or by the mere existence of K(p).
Keep pTL intertwining and width-uniform/fixed-aspect asymptotics distinct.

Use ordinary CPU and the existing bounded resource envelope. Do not extend
widths, launch production, replace proofs by another histogram, or open a
duplicate issue. #275's original-U candidate map remains a separate problem.
