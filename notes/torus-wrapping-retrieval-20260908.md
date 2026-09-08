# #642 torus wrapping retrieval (cite-or-gap)

Theory input for [issue #642](https://github.com/LightChainr/Matching-One/issues/642). Distinct from #637 (transfer matrices / critical polynomials). **Does not enter** `docs/STATUS.md`. Does not close #642. No claimed exact threshold is ingested through `scripts/threshold_claim_intake.py` because none is claimed.

Read this session, at first hand unless marked:

- Mertens–Ziff, *Percolation in finite matching lattices*, PRE **94**, 062152 (2016), [arXiv:1603.07289v2](https://arxiv.org/abs/1603.07289)
- Newman–Ziff, *A fast Monte Carlo algorithm for site or bond percolation*, PRE **64**, 016706 (2001), [cond-mat/0101295](https://arxiv.org/abs/cond-mat/0101295) (PDF)
- Pruessner–Moloney, *Winding clusters in percolation on the torus and the Möbius strip*, [cond-mat/0310361](https://arxiv.org/abs/cond-mat/0310361)
- Repository notes already on published wrapping ground: `notes/literature-officer-20260905-issue576-wrapping.md`, `notes/pinson-arguin-primitive-baseline.md`
- This session’s axis/diamond wrapping-type census: `results/wrapping-type-census/`, `notes/wrapping-type-census-l3l4-20260908.md`

Pinson, *Critical percolation on the torus*, J. Stat. Phys. **75**, 1167–**1177** (1994) is **not** opened as a PDF this session (Short Communication; no arXiv version; `hep-th/9309029` is a different paper). Formulae and ten-figure values below are quoted from Newman–Ziff and Pruessner–Moloney, who attribute them to Pinson. That is a read-through, not a primary-page check.

---

## Q1 — exact wrapping results, by type

**Universality class.** Two-dimensional percolation (Fortuin–Kasteleyn / Q=1 Potts, c=0). Pinson’s torus wrapping probabilities are **universal at criticality for a given aspect ratio / modular parameter**, not lattice-specific. Newman–Ziff evaluate them for the **square torus** (aspect 1) and use them as the L→∞ targets of **square-site** wrapping. Pruessner–Moloney check the same formulae numerically for **both site and bond** percolation on very large tori and report no excess site-versus-bond deviation from the formula.

**They do transfer to square site.** Newman–Ziff §III.A treat square **site** percolation as the working example and quote Pinson’s evaluations as the exact R_∞(p_c) for that geometry. Akhunzhanov–Eserkepov–Tarasevich (2022), already scored in the #576 wrapping note, give exact **square-site** Bernstein polynomials for wrapping along one specified direction through L=12 (L=10 ancillary file corrupt). Independent enumeration at L=3,4 matched those polynomials bit-for-bit (`notes/literature-officer-20260905-issue576-wrapping.md`).

**They are for one lattice copy, not for the matching lattice.** Pinson / Newman–Ziff / Akhunzhanov give wrapping of occupied clusters on one graph. Matching-lattice wrapping (vacant sites, NN+NNN) is **not** a second table in those papers. The matching half enters only through Mertens–Ziff (Q2–Q3).

**Type resolution, Newman–Ziff definitions** (cond-mat/0101295, p. 10). On an L×L torus:

- R^{(h)}: wraps the specified axis, and **may** also wrap the other
- R^{(e)}: wraps either axis, or both
- R^{(b)}: wraps **both** axes (cross **or** spiral; their Fig. 8)
- R^{(1)}: wraps one specified axis **and not** the other

Identities (their (11)–(12)):

```text
R^{(e)} = 2 R^{(h)} − R^{(b)}
R^{(1)} = R^{(h)} − R^{(b)} = (R^{(e)} − R^{(b)}) / 2
```

**Pinson evaluations at p_c, aspect 1**, quoted by Newman–Ziff as (13)–(15):

```text
R^{(e)}_∞(p_c) = 1 − [ϑ3(e^{−3π/8}) ϑ3(e^{−8π/3}) − ϑ3(e^{−3π/2}) ϑ3(e^{−2π/3})]
                 / (2 [η(e^{−2π})]^2)

R^{(1)}_∞(p_c) = [3 ϑ3(e^{−6π}) + ϑ3(e^{−2π/3}) − 4 ϑ3(e^{−8π/3})]
                 / (√6 [η(e^{−2π})]^2)
```

and, to ten figures,

```text
R^{(h)}_∞(p_c) = 0.521058290
R^{(e)}_∞(p_c) = 0.690473725
R^{(b)}_∞(p_c) = 0.351642855
R^{(1)}_∞(p_c) = 0.169415435
```

Duality remark, Newman–Ziff p. 10–11, verbatim in substance: R^{(e)}_∞(p_c) = 1 − π(Z×Z), because if there is no wrapping around either axis then there is a cross configuration on the **dual** lattice. That duality is **bond / dual**, not site / matching. Do not import it as a site identity.

**Neither-wrap at aspect 1.** 1 − R^{(e)} = 0.309526275. Pinson’s surprising identity, restated by Pruessner–Moloney: the probability of a **cross** topology equals the probability that **all** clusters are homotopic to a point, π(X, r) = π(0, r). At r=1 this is the same number as 1 − R^{(e)}.

**Winding-number formula** (Pruessner–Moloney eq. (1), attributed to Pinson):

```text
P̂((a,b), ≥1, r)
  = Σ_ℓ Z_{a 3ℓ,     b 3ℓ}(2/3; r)
    − ½ Σ_ℓ Z_{a(3ℓ+1), b(3ℓ+1)}
    − ½ Σ_ℓ Z_{a(3ℓ+2), b(3ℓ+2)}
    −   Σ_ℓ Z_{a 2ℓ,     b 2ℓ}
    +   Σ_ℓ Z_{a(2ℓ+1), b(2ℓ+1)}
```

with Z_{m,n}(g; r) as in that paper (g=2/3 in **this** normalization). This is Arguin’s π({a,b}), **not** Newman–Ziff R^{(h)}. At r=1, π({1,0}) = 0.169415435… = R^{(1)}, already frozen in `notes/pinson-arguin-primitive-baseline.md`. At r=2,4 they diverge: see Q5.

**Gap inside Q1.** No published exact table for **Sq8 / NN+NNN site** wrapping polynomials. Akhunzhanov is square NN only.

---

## Q2 — joint law of primal vs matching wrapping

**Cite, combinatorial, all finite L.** Mertens–Ziff §II, not a CFT formula. Notation trap: \(\hat R(1-p)\) is wrapping of the **vacant** matching colouring at black density \(p\), i.e. it is an equal-configuration (equal-\(p\)) statement. It is not a second independent copy of the primal lattice occupied at \(1-p\). On a torus, Euler’s formula plus the matching construction give, configuration-wise (their (9)–(11)):

```text
N_black − N_white − (V − E + F0)
  = +1  if black is cross-wrapping
    −1  if white is cross-wrapping
     0  otherwise
```

and the pairing

- no black wrap ⇒ exactly one white **cross**-wrapping cluster
- k single-wrapping black clusters ⇔ k single-wrapping white clusters
- black cross-wraps ⇔ white has **no** wrapping
- spiraling counts match: R^{s}(p) = R̂^{s}(1−p)
- one-direction only: R^{1}(p) = R̂^{1}(1−p)  (their (19))

**Not independent copies of Pinson.** The two wrapping events live on complementary colourings of the **same** configuration. At p=p_c both sides are critical (p_c(matching)=1−p_c), so each **marginal** tends to the Pinson numbers, but the **joint** is supported only on the pairing above.

**Exact identity for the difference**, Mertens–Ziff (20), the main theorem:

```text
M_L(p) := N_L(p) − N̂_L(1−p) − L² χ(p)
        = R^x_L(p) − R̂^x_L(1−p)
        for x ∈ {c, b, e, h}
```

with χ_□(p) = p − 2p² + p⁴ on square site. They state explicitly: **the only contribution to the right-hand side is the cross-wrapping probabilities**; the other wrapping types cancel by the pairing.

This **is** the repository observable. `scripts/exact_matching_polynomial.py`:

```text
D(C) = 1{black NN wraps} − 1{white NN+NNN wraps}
M(p) = E_p[D] = Σ_k a_k p^k (1−p)^{N−k}
```

“Wraps” here is Newman–Ziff **either** (any nontrivial homology). By MZ (20) that equals the **cross** difference. `scripts/matched_torus_reference.py` already names the equality as the Mertens–Ziff finite matching relation.

**What is not published.** A closed modular/CFT formula for the **joint** law (primal type, matching type) at p_c, beyond the combinatorial support and the one-dimensional difference M. Pinson does not treat the matching colouring. The equal-occupancy difference R(p)−R̂(p) (matching occupied at p, not vacant at 1−p) is a different observable and is also unpublished. Those are real gaps for a continuum joint. They are **not** a gap for the finite-L pairing.

**#640 is a type-resolved census of a published pairing**, not a rediscovery of an unknown identity. If the wrapping classifier implements MZ types, the 4×4 can only be supported on the five cells the pairing allows. The L=3,4 census (`notes/wrapping-type-census-l3l4-20260908.md`) finds exactly those cells and zero MZ-forbidden cells. Axis L=5 (2^25, tripwire vs PR #649 Bernstein pass) and diamond L=4 K2 (2^32, tripwire pass) repeat the same five cells; `both-two` mass is zero (#651 A-continues). D mass lives only on `neither×both` / `none×both-same` and `both×neither` / `both-same×none`. Complement-transpose of the 4×4 **fails** (NN ≠ NN+NNN), which is the finite-L source of `M(p)+M(1−p)≠0`.

**Do not cite #628’s 118133 bond dual_fail as physics.** That count is an implementation artifact, repaired in PR #653. Site `M(p)+M(1−p)≠0` still stands, and is the expected square-site (non-self-matching) statement.

**Scullard–Jacobsen connection**, MZ after (21): the criterion R^c_L(p) − R^0_L(p) = 0 is identical to M_L(p)=0, because R^0(p)=R̂^c(1−p). That sentence belongs to #637 as well; it is recorded here because it is how wrapping types become a threshold polynomial.

---

## Q3 — the matching function D(C) = 1{black wraps} − 1{white wraps}

**Where introduced.** Mertens–Ziff 2016. Two faces of the same object:

1. Cluster-count form (15): M_L(p) = N_L(p) − N̂_L(1−p) − L² χ(p)
2. Wrapping form (20): M_L(p) = R^x_L(p) − R̂^x_L(1−p), x ∈ {c,b,e,h}

The repository uses (2) with x = either. They prove (1)=(2) at **every finite L** on a torus.

**What is proved about the root.**

- M_L is strictly increasing, range in [−1,1] (because it equals a difference of probabilities).
- Unique root p*_L ∈ (0,1).
- lim_{L→∞} M_L(p) = −1 for p<p_c and +1 for p>p_c (their (31)). Hence p*_L → p_c. This is the same qualitative convergence the repository already has from subcritical decay + matching duality (`notes/homological-balance-root-ledger-20260906.md` §2; `docs/astra/ANSWER-610-20260907.md`).
- On **self-matching** lattices (triangular site, square bond, …), M_L(p_c)=0 for **every** L, so p*_L = p_c with zero displacement (MZ (22), (24)). The L^{−4} displacement is a square-site matching-odd residual. The repository already isolated that as a theorem in the homological-balance note.

**Rate.** **Not a theorem.** MZ: empirically p*_L − p_c ∼ L^{−w} with w≈4, citing Jacobsen 2014/2015; their own fit from M_L(p_c)∼L^{2−x} gives w=2−x−1/ν≈4.17, and a direct plot of p*_L−p_c gives slope −4.07. They say larger L is needed. Jacobsen’s transfer-matrix critical polynomials are the high-precision engine behind the “~L^{−4}” lore; that engine is #637, not a wrapping-rate theorem.

#618’s report that Q_N(u)→p_c gives **no polynomial rate from H1–H3** is compatible: MZ also do not prove a rate. The empirical L^{−4} is an extra input, not a corollary of the finite identity.

**Do not promote L^{−4} through threshold_claim_intake.** It is not an exact threshold.

---

## Q4 — finite-size corrections: wrapping vs the difference

**Wrapping probabilities themselves.** Newman–Ziff Fig. 10 and eq. (18): they **conjecture** R_L(p_c) − R_∞ ∼ L^{−2} (fits −1.95(17) site, −2.003(5) bond for R^{(1)}). Combined with the critical-window slope L^{1/ν}, the estimator defined by R_L(p)=R_∞(p_c) converges as

```text
p_L − p_c  ∼  L^{−2 − 1/ν} = L^{−11/4}
```

This is the wrapping-only rate. It is **not** the matching-function rate.

**The difference M.** MZ scaling (36)–(39): in the scaling limit M_L(p)=f(z)−f(−z), z∝(p−p_c)L^{1/ν}. Even powers cancel, so M is analytic in z even at criticality. Corrections: M_L(p_c)∼L^{2−x} with a numerical 2−x≈−3.42, hence the root shifts as L^{2−x−1/ν}≈L^{−4.17}. If one assumes w=4 exactly, 2−x=13/4, which is the L^{−13/4} the repository already writes for M_L(p_c).

**Shape after location/scale.** The repository’s ω≈0.970±0.031 (`results/p612-n725-score/latest.json`, spin0 exponent fit; unity not excluded) is the residual of the **threshold-law shape** after projecting out location and scale. MZ do **not** quote a correction exponent for that projected shape.

A literature comparison, not a new fit: Newman–Ziff’s R_L(p_c)−R_∞ ∼ L^{−2} is N^{−1} in site count N=L². A leftover N^{−1} after location/scale would be ω=1. The measured 0.97±0.03 does not exclude unity. That is a **possible** identification, not a published theorem that the matching-function shape correction is the wrapping θ=2.

**Gap.** No paper found that analyses the affine-invariant shape of M_L, or of Q_N, or quotes ω≈0.97. Q4 is cite for the **unprojected** wrapping and matching-root corrections, gap for the projected shape.

---

## Q5 — aspect ratio vs the N=580 ladder

Pinson’s wrapping probabilities depend on the modular parameter. The repository already evaluated π({1,0})(ir) at r=1,2,4 (`notes/literature-officer-20260905-issue576-wrapping.md`):

| r | π({1,0})(ir) | ratio to r=1 |
|--:|--:|--:|
| 1 | 0.169415435321 | 1 |
| 2 | 0.503035897695 | **2.969244784222** |
| 4 | 0.855969321054 | **5.052487215408** |

Newman–Ziff R^{(h)}(i)=0.521058290 is a **different** observable (specified-direction wrap, including simultaneous wrap). Do not score the ladder against 0.521 when the competitor is π({1,0}).

**Does this predict the aspect-ladder amplitude ratio?** Only if the measured object **is** a Pinson wrapping. The #576 wrapping note already recorded a **non-claim**: matching-odd slope is not identified with π({1,0}). N=290’s measured 1.880±0.177 is ~6σ from 2.969, so Pinson does not explain that run. It remains a named competitor (`pinson_pi10_ratio`) that must sit on the next freeze **before** the data, alongside:

| competitor | r=2 | r=4 |
|---|--:|--:|
| weight-4 Ê4(ri)/Ê4(i) | 2.75 | 10.99 |
| Pinson π({1,0}) ratio | 2.969 | 5.052 |
| bare aspect r | 2 | 4 |
| area r² | 4 | 16 |

The three 11/4’s that must stay apart: modular weight-4 ratio 11/4; Newman–Ziff estimator L^{−11/4}; Pinson r=2 ratio 2.969.

**Score.** Theory number exists for **wrapping homology class {1,0}**. It does **not**, on present evidence, predict the matching-odd aspect ladder. Converting the underpowered three-hypothesis result into a Pinson test would first require identifying the ladder’s observable with π({1,0}) — which #576 explicitly refused. That identification is still refused here.

If a later readout **is** wrapping-flavoured, score it in the Cardy/Pinson function space, not against E4, and use the table above.

---

## What this does to the live tickets

- **#640 / #651.** The finite-L 4×4 is the MZ pairing, resolved by wrapping type. Axis L=5 is a **classifier tripwire** (does a sixth cell appear?) more than an open topological question. both-two mass is the one structural caveat MZ do not name in wrapping-type language; it is the rank-2-versus-two-component distinction inside `both`.
- **#635.** MZ prove, at every finite L, that M is exactly the difference of two **cross-wrapping** amplitudes. Verdict A in the wrapping-form is published. Whether that is the same map as Jacobsen’s two transfer-matrix sectors is #637/#635, not a wrapping-probability gap. #646’s “degenerate both-same” labels need to be matched onto MZ cross vs spiral; the L=3,4 census is compatible with that degeneracy and does not prove it.
- **#618 / #622.** Location Q_N(u)→p_c is proved (qualitative). Polynomial rate of the root is empirical L^{−4}, not a wrapping theorem. Projected shape ω is a gap relative to this literature.
- **#577 / aspect ladder.** Pinson supplies a number; it is not, without a new identification, the ladder’s theory value.

---

## Boundaries

- No STATUS edit, no ticket closed, no threshold ingested.
- #637 still owns: connectivity transfer matrices, critical polynomials, Jacobsen eigenvalue identities, Scullard–Jacobsen as a computational method.
- #620 still owns: inverse-CDF windows, DKS, F1 connecting lemma.
- #601 still owns: symmetry quotients and selection rules.
