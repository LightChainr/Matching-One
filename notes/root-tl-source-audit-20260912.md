# Source and claim audit: root consistency / P398 TL identification

2026-09-12. Status PRIMARY_TEXT_READ below means the stated sections and equations
were actually opened. It does not mean every proof in each paper was independently
reproved. No broad novelty search was performed.

## Repository inputs actually read

- #716, `8f60bd28e72eb2e00399bdb8f5528190dfb0f045`, thin-torus result and source modes.
- #715, `308f66c7c6e34af7feb4fb689f27dfccea3eefea`,
  `notes/lit-double-pulse-observability-20260912.md`, especially Q3a.
  Review comment `5646126309` repeats the negative P398/TL identification.
- #709, `83e011a7ceeab23b9ae5eab741f41ea15c983459`,
  `scripts/p398_double_pulse_exact.py`, plus the already supplied unmodified archive.
- `main@eb89e9422791d9e3c3a78f0e65d56912b815a7bd`,
  `scripts/planar_state_operations.py`: exact definitions of detach and cyclic join.
- #711--#714 retrieval tickets; #276/#613 probability contract. No replacement
  retrieval tasks were created and no issue/PR lifecycle was changed.

The reused #709 script copy has Git blob
`e211c5cf9e3677dd7f95c27deec5402723361f1b`, SHA256
`975bd781a27a4e2f03fe84f15669b1cd8ac0b9030a6c86dd3cc378df6825bb28`.
The standalone archive check verifies that copy before comparing G/H; it is not
an assertion that a complete repository checkout or complete CI was executed.

## Primary probability sources

**Duminil-Copin--Tassion**, *A new proof of the sharpness of the phase transition
for Bernoulli percolation and the Ising model*, arXiv:1502.03050v3:
https://arxiv.org/html/1502.03050v3
PRIMARY_TEXT_READ: Theorem 1.1, §1.2, specifically the paragraph Site percolation.
The finite-range exponential-decay conclusion and explicit site adaptation are
used. The printed bond-square pc=1/2 discussion is NOT used for square-site pc.
The present finite-torus probability-ratio bound is derived in our note; it is
not attributed as a theorem appearing in this source.

**Grimmett--Li**, *Percolation critical probabilities of matching lattice-pairs*,
Random Structures & Algorithms 65 (2024) 832--856, DOI 10.1002/rsa.21226:
https://onlinelibrary.wiley.com/doi/10.1002/rsa.21226
PRIMARY_TEXT_READ: introductory Eqs. (1.1),(1.3), amenable pc=pu passage,
and theorem context. The introduction carries its own provenance for the matching
relation. The square NN graph and NN+NNN matching graph satisfy the relevant
one-ended/transitive/amenable planar-primal setting. No nonamenable equality of
pc and pu is imported.

## Primary process/representation sources

**Pearce--Rittenberg--de Gier--Nienhuis**, *Temperley-Lieb Stochastic Processes*,
J. Phys. A 35 (2002) L661--L668, arXiv:math-ph/0209017v2:
https://arxiv.org/html/math-ph/0209017v2
PRIMARY_TEXT_READ: §2 Eqs. (2.1),(2.4),(2.15)--(2.17), IC/DC distinction;
§3 stationary-conjecture context. The representation required here is periodic IC
at even endpoint number L=2w, dimension Catalan(w), not the larger DC space.
Its statement about irreversible joining in line filtrations does not rule out
point detach in a different state encoding. We establish the needed map directly.

**Cantini--Sportiello**, *Proof of the Razumov-Stroganov conjecture*,
J. Combin. Theory A 118 (2011) 1549--1574, arXiv:1003.3376v1:
https://arxiv.org/html/1003.3376v1
PRIMARY_TEXT_READ: §2.1 ASM/FPL normalization, §2.2 Eq. (4) reconnection,
§2.4 Eqs. (22)--(24), and proof statement in §3. The static FPL boundary-pattern
count vector is an eigenvector of sum e_i with eigenvalue 2w. The previously
conjectural periodic IC stationary correspondence is proved here. We import
that theorem after identifying the generator; we do not claim to have rechecked
the entire combinatorial proof or derived the alternating-rate stationary law.

The arXiv HTMLs display a recent auto-rendered date in their body. Bibliographic
attribution uses the original arXiv version/publication dates above, not the
HTML rendering timestamp. No PDF parsing or figure interpretation was required.

## Limits of the audit

No originality certification for the fattening map, complement identities,
Hankel theory, the stationary law, or the all-aspect root consequence. The first
four are established structural ingredients. The last consequence has an explicit
proof and a precisely scoped novelty/retrieval question for the team.
Not observing a process name in a keyword search is not a no-go; the #715 mistake
is resolved by an actual conjugacy, not by finding a more similar abstract.
