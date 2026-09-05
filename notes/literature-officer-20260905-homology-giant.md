# Literature officer — giant cycles / wrapping homology, 2026-09-05

Repo search for `2011.11903` / Duncan / Kahle returned **zero**. The tickets name matching-odd homology; this is the paper that defines giant cycles on the torus. A 2026 follow-up splits the observable the engine might actually be counting. Theory input; **does not enter** `docs/STATUS.md`.

X this pass: one July 2026 post on CSS codes measuring “homological percolation”; nothing on square-site wrapping or matching-odd.

Do not re-read as new: Pinson 1994, Pruessner–Moloney, Akhunzhanov polynomials, Arguin (already in the deeper note), #576 `pinson_pi10_ratio`.

---

## Giant 1-cycle is wrapping; `p_c=1/2` is not square site

**Duncan, Kahle, Schweinhart,** [arXiv:2011.11903](https://arxiv.org/abs/2011.11903) v4 (29 Sep 2023), *Ann. Inst. H. Poincaré Probab. Statist.* (2025). Plaquette percolation (high-d analogue of square **bond**) and permutohedral site percolation (analogue of triangular **site**) on `T^d`.

A cycle in the random subcomplex is **giant** when its image in `H_i(T^d)` is nontrivial. Sharp transition, for every `i,d`, from no giant cycles to giant cycles spanning the homology. In middle dimension `i=d/2`:

> `p_c = 1/2`

for both models. Finite-volume analogue of Kesten: square **bond** and triangular **site**.

**Square site is not in that list.** Wrapping on the square-site torus still has a threshold (it is `p_c(Sq) ≈ 0.5927`), but that threshold is not `1/2` and is not a theorem of this paper. Matching (Sq / Sq8) is a different involution; Duncan–Kahle do not treat it.

For `d=2,i=1`, a giant 1-cycle **is** a wrapping path on `T^2`. That is the homological name of the #576 channel, not a new number.

---

## 2026: event A ≠ event S

**Schweinhart and Shuman,** [arXiv:2601.00793](https://arxiv.org/abs/2601.00793) (2026). *Voronoi percolation: topological stability and giant cycles.* Voronoi (Poisson), not the square lattice. Same homological events as Duncan–Kahle, now with a 2-d dictionary:

- **Event A:** an `i`-dimensional giant cycle exists. For `i=1` on `T^2`: **one** periodic path (e.g. vertical wrapping).
- **Event S:** a basis of non-homologous giant cycles. For `i=1` on `T^2`: **two** non-homologous wrappings (horizontal **and** vertical).

On the `2i`-torus, `ℙ_p(A)→0` for `p<1/2` and `ℙ_p(S)→1` for `p>1/2`.

**Why #576 has to name A vs S.** Pinson `π({1,0})` (and Newman–Ziff `R_h`) is A-type: wrapping in a specified direction, evaluated **at criticality**, not as a threshold. Both-direction wrapping `R_b` is closer to S. They are different numbers (Pinson `R_h(r=1) ≈ 0.521058290` vs `R_b` smaller). The #576 note already computed `π({1,0})` at `r=1,2,4`. If the engine wrapping channel is “any wrapping” vs “specified direction” vs “both axes”, those three are A, A-specified, and S. Comparing any of them to `E₄` `11/4` is still the wrong function space (Ridout / Diamantis–Kleban / SKZ). Comparing A to S without naming which is how a later ladder fights the published number.

Matching-odd is **not** claimed to be A or S. If a later readout *is* wrapping-flavoured, write which event, then score it against Pinson / Akhunzhanov, not against `11/4`.

---

## Opinions

1. **#576 / #567.** Freeze the wrapping observable as A-specified (`π({1,0})` / `R_h`), A-any, or S. The published numbers already exist for A-specified. Do not leave it unnamed.
2. **Matching-odd homology ticket.** Duncan–Kahle is giant cycles, not matching-odd. Matching-odd remains the Sq/Sq8 involution on the engine’s homology readout. Do not cite `p_c=1/2` as if it applied to square site.
3. **Q2.** Unchanged: wrapping-flavoured leftovers live with Pinson / Cardy–Watts / Ridout-`t`, not holomorphic `E₄`.
4. **P2.** Unchanged. Homological `p_c=1/2` is bond/triangular/Voronoi self-dual, not an algebraic form for square site.

## Not established

- matching-odd = giant 1-cycle;
- square-site wrapping threshold equals `1/2`;
- anything in the claim ledger.
