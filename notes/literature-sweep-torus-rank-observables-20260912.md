# Literature sweep: torus rank observables, exact finite-size site polynomials, matching dictionary

2026-09-12. Executed at the owner's request ("broadly search the related
literature") as a bounded novelty/prior-art sweep around PR #733 and the
#711/#717 dictionary correction. It is a **retrieval note**, not a result.
Nothing here is a claim about our data, and nothing here makes any output
publication-ready.

## 0. Scope and how this was searched

Sources: arXiv API (`export.arxiv.org/api/query`, relevance-sorted,
`max_results` 8–10 per query) and general web search. Date: 2026-09-12.

Query strings actually used (arXiv `all:` / `abs:` fields):

```text
all:"critical polynomial" AND all:"percolation"
all:"wrapping probability" AND all:"percolation"
all:"percolation" AND all:"homology" AND all:"torus"
all:"site percolation" AND all:"transfer matrix"
all:"matching lattice" AND all:"percolation"
all:"matching polynomial" AND all:"percolation"
all:"Pinson" AND all:"percolation"
all:"percolation" AND all:"conditional sampling"
all:"percolation" AND all:"winding" AND all:"torus"
all:"homological percolation"
all:"percolation" AND all:"Euler characteristic" AND all:"torus"
all:"Arguin" AND all:"percolation" AND all:"torus"
all:"Scullard" AND all:"critical polynomial"
```

Plus targeted web queries on: Mertens–Ziff finite matching lattices; torus
wrapping/homology probabilities; Pinson's analytic formula; exact site
percolation probabilities on plane/cylinder/torus; Jacobsen's site-percolation
critical polynomials; homological percolation on a torus; Newman–Ziff wrapping
estimators.

**Search boundary (must be stated wherever this is cited).** No MathSciNet,
zbMATH, Scopus or Google Scholar citation-graph traversal; no non-English
literature; no books or theses; no patent search; no systematic forward
citation chasing of Pinson (1994) or Arguin (2001). Several entries below are
known only through secondary quoting in arXiv preprints, and are flagged as
such. **`NO_MATCH_FOUND_IN_THIS_SEARCH` below means no match in this sweep, not
a novelty proof.**

## 1. The one finding that matters most: our rank variable already has a name

The central observable of this repository is

```text
r(omega) = rank im[ H_1(K(omega)) -> H_1(T^2) ] in {0,1,2},   X = r - 1.
```

The **inclusion-induced map on first homology of the occupied subcomplex into
the ambient torus is exactly the object** studied, under different names, by:

- **Langlands, Pouliot, Saint-Aubin (1994)** — introduced crossing
  probabilities on compact Riemann surfaces; for percolation on a surface `S`
  they consider `phi: H_1(X_s) -> H_1(S)` and the probabilities `pi_G` that the
  image is a given subgroup `G <= H_1(S)`.
  (Quoted in arXiv:1402.0879 and in arXiv:2011.11903v4 §1 as [LPSA94].)
- **H. T. Pinson, "Critical spanning probability on a torus", J. Stat. Phys. 75,
  1167 (1994)** — analytic expressions for those subgroup probabilities as
  functions of the torus modular parameter, via an orientation argument on
  cluster-boundary curves; rigorous except for the mesh-to-zero limit, which
  uses a Nienhuis renormalization-group step.
  (Quoted throughout arXiv:0812.2925v2 §1 as [17] and arXiv:0905.3521.)
- **L.-P. Arguin, "Homology of Fortuin–Kasteleyn clusters of Potts models on the
  torus", arXiv:hep-th/0111193** — extends Pinson to FK clusters of the
  Q-state Potts model, `Q in [1,4]`, with closed forms in Jacobi theta
  functions for `Q = 1` (percolation) and `Q = 2`; `Q = 1` is our case.
- **Morin-Duchesne & Saint-Aubin, arXiv:0812.2925v2** — asymptotic behaviour of
  `pi({a,b})` for thin tori; exponents tied to weights `h_{r,s}` of the extended
  Kac table with half-integer entries as well.
- **Duncan, Kahle & Schweinhart, "Homological percolation on a torus: plaquettes
  and permutohedra", arXiv:2011.11903v4; Ann. IHP 61, 2235 (2025)** — for a
  random subcomplex `S` of `T^d` they take the natural inclusion
  `phi: S -> T^d` and call nontrivial elements of
  `im[ phi_*: H_i(S;Q) -> H_i(T^d;Q) ]` **giant cycles**. They prove a sharp
  transition (via Friedgut–Kalai), convergence of the threshold function, and
  `p_c = 1/2` in the middle dimension `i = d/2`. They explicitly cite
  [LPSA94], [Pin94] and [MDSA09] as the 2-dimensional precursors.

**Dictionary that follows immediately, and that I have not seen written down in
our notes:**

```text
our P0  <->  pi({0})                       (trivial ambient image)
our P1  <->  sum over coprime (a,b) of pi({a,b})   (rank-one primitive family)
our P2  <->  pi(Z x Z)                     (rank-two / cross topology)
```

so that `A_top = <X> = P2 - P0` is a difference of two Pinson/Arguin subgroup
probabilities, and `E_top = P2 + P0` is their sum. If that identification is
correct at the level of events, then **the scaling-limit values of our two
coordinates are already known in closed form (Jacobi theta / Dedekind eta)**,
and the finite-size object we have been fitting is the correction to them.

Two cautions that must travel with this:

1. Pinson/Arguin is a **scaling-limit** statement (plus one non-rigorous
   renormalization step). Our numbers are **exact finite-size**. Agreement is
   expected only after finite-size corrections, and the *shape* of those
   corrections is precisely what this repository has been arguing about.
2. The subgroup lattice of `H_1(T^2)` is what our `r in {0,1,2}` is recording.
   This is prior art for the *observable*, not for our exact finite-width
   computations, and not for the finite-amplitude source work in #733.

**Verdict: PRIOR_ART_FOUND for the observable definition and (in the continuum
limit) for its distribution. This is the single most important thing to tell
any referee, and it is better coming from us.**

### 1a. The thin-torus limit is our width-4 regime

Morin-Duchesne & Saint-Aubin compute the asymptotics of `pi({1,0})` along
`tau_r = 0, tau_i -> infinity`, i.e. exactly the geometry of a **fixed-width
strip whose length grows** — the same aspect-ratio regime as our width-3/width-4
closures. Their exponents come from Kac-table weights including half-integer
indices. Any claim we make about an orientation-sensitive, aspect-dependent
finite-size exponent should be checked against that asymptotic family before
being called new.

### 1b. Euler characteristic zero vs homological threshold

**Bobrowski & Skraba, Phys. Rev. E 101, 032304 (2020), arXiv:1910.10146** —
numerical evidence across four models (site percolation on cubical and
permutahedral lattices, Poisson–Boolean, Gaussian random fields; flat torus,
`d = 2,3,4`) that **the zeros of the expected Euler characteristic curve
approximate the homological-percolation thresholds**, with a discussion of the
approximation error.

This is directly adjacent to our use of the Sykes–Essam matching polynomial
(`chi(p)`) and of `M_B(p) = P2 - P0` as a threshold functional. It is a
different model class (higher-dimensional / continuum) and it is numerical, but
"EC-zero approximates the topological threshold" is now an explicit, published
heuristic with a named error. **Verdict: PARTIAL_OVERLAP; cite it when we use
an Euler-characteristic zero as a threshold estimator.**

## 2. Exact finite-size site percolation: who has computed what, to what size

This is the competitive landscape for PR #733's width-4 exact rank closure and
for #708's automaton.

| Work | Object | Reach | Relation to us |
|---|---|---|---|
| **Akhunzhanov, Eserkepov, Tarasevich, J. Phys. A 55, 204004 (2022), arXiv:2204.01517** | exact percolation-probability **polynomials**, site percolation on `L x L`: plane (crossing), cylinder (spanning), **torus (wrapping along one direction)** | `L <= 17` plane, `L <= 16` cylinder, **`L <= 12` torus**; dynamic programming + topology-based state reduction; divisibility properties proved; naive FSS gives `p_c = 0.59269` | Closest competitor. They reach full two-dimensional `L x L` tori; we reach **fixed width 4 at arbitrary length** with a rank decomposition they do not make. Their `R^{(e)}, R^{(1)}, R^{(b)}, R^{(h)}, R^{(v)}` are *direction* events; our `r` is an *ambient-rank* event. The map between them is a dictionary exercise we have not done. |
| **Mertens, J. Phys. A 55, 244002 (2022), arXiv:2109.12102** (already in `references.bib`) | exact `R_n(p)` for spanning an `n x n` square (open boundaries) | `n <= 24`, `O(lambda^n)`, `lambda ~ 2.6` | Open boundary, spanning — not torus, not rank. Prior art for the *technique* of exact enumeration + extrapolation, not for our observable. |
| **Newman & Ziff, PRL 85, 4104 (2000), arXiv:cond-mat/0005264** | microcanonical union-find MC; `p_c = 0.59274621(13)` site/square | MC, not exact | The number everyone compares against. |
| **Newman & Ziff, cond-mat/0203496 ("Convergence of threshold estimates")** | exact enumeration of `R_L(p)` for site percolation on `L x L`, crossing | Table up to `L = 7` (L=2..5 from Reynolds–Stanley–Klein; 6,7 from Ziff 1992) | Small-size exact crossing polynomials; superseded in reach but useful as an independent check of any small-L polynomial we produce. |
| **Jacobsen, J. Phys. A 47, 135001 (2014), arXiv:1401.7847** | graph polynomials `P_B(q,v)` and site polynomials `P_B(p)` via a transfer matrix in the **periodic Temperley–Lieb algebra**; "We discuss in detail the role of the symmetries and the embedding of `B`" | bases to 882 edges (bond), site polynomials with up to 243 vertices; `p_c` to ~1e-8..1e-9 | **The symmetry-and-embedding discussion is the direct precedent for our D4/lumping work.** Different algebra (pTL vs our finite automaton), different target (thresholds vs rank decomposition). |
| **Jacobsen, arXiv:1507.03027** (in `references.bib`) | eigenvalue identities in periodic TL; semi-infinite cylinders of circumference `n` | site square `n_max = 21`; `p_c = 0.59274605079210(2)` | Already read for #681/#717. |
| **Scullard, arXiv:1111.1061 ("The percolation critical polynomial as a graph invariant")** | critical polynomial as graph invariant | — | Background for the #717 dictionary question; not yet read at section level here. |
| **arXiv:2010.02887** ("Critical polynomials in the nonplanar and continuum percolation models") | non-planar / continuum extensions | — | Flagged, unread at section level. |
| **arXiv:2205.02734** ("Percolation critical probabilities of matching lattice-pairs") | matching-lattice pair thresholds | — | Flagged, unread at section level; potentially relevant to the `p_c + p_c(NN+NNN) = 1` anchor. |

**Verdict for this group:** the *exact finite-size* frontier on full `L x L`
tori is `L <= 12` (AET 2022). Our contribution is not "bigger `L`" — it is a
different decomposition (ambient rank) and a fixed-width/arbitrary-length
regime. Say it that way.

## 3. Matching lattice / critical polynomial dictionary (the #711/#717 channel)

- **Mertens & Ziff, Phys. Rev. E 94, 062152 (2016), arXiv:1603.07289** — already
  PRIMARY_TEXT_READ in `notes/finite-critical-polynomial-dictionary-20260912.md`.
  Finite-size generalisation of Sykes–Essam:
  `N_L(p) - N̂_L(1-p) - L^2 chi(p) = R_L^x(p) - R̂_L^x(1-p)`.
- **Mertens, Jensen & Ziff, Phys. Rev. E 96, 052119 (2017), arXiv:1602.00644** —
  cluster-number universality; uses the Sykes–Essam matching polynomial for
  exact lattice/matching-lattice relations. Not yet read at section level.
- **Sykes & Essam (1964)** — in `references.bib`. `chi(p) = p - 2p^2 + p^4` for
  the square/NN+NNN pair; `chi_Delta(p) = p - 3p^2 + 2p^3` for self-matching.
- **arXiv:1906.10543** ("Critical p = 1/2 in percolation on semi-infinite
  strips") — flagged, unread.

**Verdict:** the dictionary note in #733 is consistent with the primary texts,
and this sweep found **nothing that contradicts it**. It did find that the
surrounding literature is larger than the three texts we read, in particular
on the matching-lattice-pair side.

## 4. Symmetry reduction and lumping

- **May & Wierman, Combin. Probab. Comput. 14, 549 (2005)** — use the graph
  **automorphism group** to reduce the computational work of the substitution
  method; on `(3,12^2)` bond percolation they cut the bound interval width by
  62%. Methodologically the closest published relative of "reduce by symmetry
  before you enumerate" — but it is substitution-method bound improvement, not
  a transfer-matrix lumping, and it does not produce state counts.
- **Jacobsen (2014), arXiv:1401.7847** — explicit discussion of the role of
  symmetries and of the embedding of the basis in the transfer-matrix
  construction.
- No arXiv hit for `abs:"lumping" AND abs:"Markov" AND all:"percolation"`.

**Verdict:** the *idea* of symmetry reduction is standard; the specific result
in #733 — seven D4 column-probability types, common strong lumpings of size
94/303/179/262/509/303/509, each equal blockwise to the colour-preserving D4
orbit partition — has **NO_MATCH_FOUND_IN_THIS_SEARCH**. That is a weak
statement and should be re-checked against the transfer-matrix literature
(Jacobsen, Jensen, Enting, Guttmann) before being used in any novelty claim.

## 5. Site sources and conditioned/rare-sector sampling

- `all:"percolation" AND all:"conditional sampling"` returned **zero** arXiv
  hits.
- Web searches for a point/dipole site-source response of the specific kind in
  #733 (zero linear response, nonzero mixed response; four-sign finite-amplitude
  extraction by degree-<=2 multiaffinity) returned nothing on topic — the
  "dipole" hits are electromagnetic-wave localization, unrelated.
- Exact backward/conditional samplers that condition on a **topological**
  final state (our final-rank conditioning) have **NO_MATCH_FOUND_IN_THIS_SEARCH**.

**Verdict: NO_MATCH_FOUND_IN_THIS_SEARCH for both.** Given the search boundary
in §0, this is the weakest kind of evidence and must not be phrased as novelty.
It does suggest these two are the least contested parts of #733, which is an
argument for writing them up first — not for claiming them.

## 6. What I would do with this, in priority order

1. **Check the Pinson/Arguin dictionary numerically before anything else.**
   Take the `Q = 1` closed forms (`pi({0})`, `pi(Z x Z)`, `pi({a,b})` in Jacobi
   theta functions, as restated with full formulas in arXiv:0905.3521 §"crossing
   probabilities", which gives `P_{a,b}(r)` and `P_X(r)` explicitly in terms of
   `Z_{m,n}(g,r)` and the Dedekind eta function). Evaluate at our aspect ratios.
   Compare with the width-4 exact `P0/P1/P2`. Report the finite-size discrepancy
   as the object of interest. **This converts "we think there is a
   matching-odd structure" into "here is our measured deviation from the known
   continuum answer", which is a much stronger sentence.**
2. **Cross-validate against AET 2022 (arXiv:2204.01517) wherever our geometries
   overlap** (their torus polynomials reach `L = 12`; supplemental material
   contains the polynomials). Even a single shared small case is worth more than
   an internal consistency check.
3. **Read arXiv:2205.02734 and arXiv:1111.1061 at section level** before any
   further statement about the matching-`p_c` anchor or about critical-polynomial
   novelty. Both are cheap and both sit directly on the #711/#717 channel.
4. **Do not** phrase the rank observable as new. Phrase the finite-width exact
   rank decomposition, the source response and the conditioned oracle as the
   contributions, and cite Pinson/Arguin/DKS as the definition's origin.

## 7. BibTeX for the entries not yet in `references.bib`

Not added to `references.bib` by this note; copied here so the decision to add
is explicit and separate.

```bibtex
@article{pinson1994,
  author  = {Pinson, H. T.},
  title   = {Critical spanning probability on a torus},
  journal = {Journal of Statistical Physics},
  volume  = {75}, pages = {1167--1177}, year = {1994}}
  % known to us only via secondary quoting (arXiv:0812.2925, 0905.3521, 1402.0879)

@article{arguin2001,
  author       = {Arguin, L.-P.},
  title        = {Homology of {Fortuin--Kasteleyn} clusters of {Potts} models on the torus},
  journal      = {Journal of Statistical Physics}, year = {2002},
  eprint       = {hep-th/0111193}, archivePrefix = {arXiv}}

@article{mdsa2009,
  author       = {Morin-Duchesne, Alexi and Saint-Aubin, Yvan},
  title        = {Critical exponents for the homology of {Fortuin--Kasteleyn} clusters on a torus},
  eprint       = {0812.2925}, archivePrefix = {arXiv}}

@article{dks2025,
  author       = {Duncan, Paul and Kahle, Matthew and Schweinhart, Benjamin},
  title        = {Homological percolation on a torus: plaquettes and permutohedra},
  journal      = {Annales de l'Institut Henri Poincar\'e, Probabilit\'es et Statistiques},
  volume       = {61}, number = {3}, pages = {2235--2261}, year = {2025},
  eprint       = {2011.11903}, archivePrefix = {arXiv}}

@article{bobrowski2020,
  author       = {Bobrowski, Omer and Skraba, Primoz},
  title        = {Homological percolation and the {Euler} characteristic},
  journal      = {Physical Review E}, volume = {101}, pages = {032304}, year = {2020},
  eprint       = {1910.10146}, archivePrefix = {arXiv}}

@article{aet2022,
  author       = {Akhunzhanov, R. K. and Eserkepov, A. V. and Tarasevich, Yu. Yu.},
  title        = {Exact percolation probabilities for a square lattice: site percolation on a plane, cylinder, and torus},
  journal      = {Journal of Physics A: Mathematical and Theoretical},
  volume       = {55}, pages = {204004}, year = {2022},
  eprint       = {2204.01517}, archivePrefix = {arXiv}}

@article{jacobsen2014,
  author       = {Jacobsen, Jesper Lykke},
  title        = {High-precision percolation thresholds and {Potts}-model critical manifolds from graph polynomials},
  journal      = {Journal of Physics A: Mathematical and Theoretical},
  volume       = {47}, pages = {135001}, year = {2014},
  eprint       = {1401.7847}, archivePrefix = {arXiv}}

@article{newmanziff2000,
  author       = {Newman, M. E. J. and Ziff, R. M.},
  title        = {Efficient {Monte Carlo} algorithm and high-precision results for percolation},
  journal      = {Physical Review Letters}, volume = {85}, pages = {4104}, year = {2000},
  eprint       = {cond-mat/0005264}, archivePrefix = {arXiv}}

@article{newmanziff2002,
  author       = {Newman, M. E. J. and Ziff, R. M.},
  title        = {Convergence of threshold estimates for two-dimensional percolation},
  eprint       = {cond-mat/0203496}, archivePrefix = {arXiv}}

@article{scullard2011,
  author       = {Scullard, Christian R.},
  title        = {The percolation critical polynomial as a graph invariant},
  eprint       = {1111.1061}, archivePrefix = {arXiv}}

@article{wierman2005,
  author       = {May, William D. and Wierman, John C.},
  title        = {Using symmetry to improve percolation threshold bounds},
  journal      = {Combinatorics, Probability and Computing},
  volume       = {14}, pages = {549--566}, year = {2005}}
```

## 8. Standing cautions

- A missing formula in the texts read is not proof that no prior formula exists.
  §5 rests entirely on this negative and should be treated as provisional.
- Everything above is literature positioning. No number in this note has been
  recomputed against our own data; §6.1 is the first place where that must
  happen.
- This note does not authorize new acquisition, a next-width scan, a new
  transfer engine, or any STATUS promotion.
