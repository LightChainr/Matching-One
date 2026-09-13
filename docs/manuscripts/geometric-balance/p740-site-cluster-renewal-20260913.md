# #740 — the site-cluster renewal: sources, the identity that is missing, and a measured obstruction

Issue #740, returned to #739. Scope: the once-per-COMPONENT cylinder intensity
`nu_w^G(p)` for independent SITE percolation on `G = NN` and `G = NN+NNN`, and the
question whether `nu_w^G(p) = A_G(p) w^{-1/2} exp[-w kappa_G(p)] (1+o(1))` uniformly on
compact subcritical `p` intervals.

Answer in one line: **the analyticity half is answerable and is answered by citation; the
`1/2` half is not, and the obstruction is named and measured here rather than argued
around.** Deliverable class (c) plus a partial (a).

Two things are NOT done and are not pretended: no `A`, no `beta` is derived, and no
counterexample to the `1/2` is exhibited — what is exhibited is the specific multiplicity
that any proof of it must control first.

---

## 1. Retrieval matrix

Each row states the model, the hypotheses, what the source actually gives, and — explicitly —
what it does not.

| source | model | hypotheses | what it gives | what it does **not** give |
|---|---|---|---|---|
| Campanino–Ioffe, *Ann. Probab.* **30** (2002) 652–682, doi 10.1214/aop/1023481005 | Bernoulli **bond** percolation on `Z^d`, `d >= 2` | `p < p_c(d)`; bond; nearest neighbour | a precise Ornstein–Zernike asymptotic for the two-point function `P_p(0 <-> x)` **in any direction `x`** and any subcritical `p` | not site; not the ANN+NNN range; not a cylinder-component density; no amplitude for anything but the two-point function |
| Campanino–Ioffe–Velenik, arXiv:math/0610100 (= mp_arc 06-275) | subcritical random-cluster measures, general `q` | Assumption (1.2): exponential decay of finite-volume **wired** connectivities in rectangles. Known to hold for `q = 1`, `q = 2` **in any dimension**, and for `q` sufficiently large; in `d = 2` it holds whenever infinite-volume connectivities decay exponentially | sharp OZ two-point asymptotics; **analyticity and strict convexity of the inverse correlation length**; an invariance principle; and structurally, a description of long clusters as **"essentially one-dimensional chains of irreducible objects"** with a random-walk representation | not a component-density theorem; the amplitude is the connectivity one, i.e. a *linear* chain between two distant points |
| D'Alimonte–Manolescu, arXiv:2510.13648v3 (23 Jun 2026) | 2D random-cluster, `1 <= q < 4` | 2D; `1 <= q < 4`; random-cluster (bond FK) | an OZ asymptotic for the two-point function holding **uniformly for `p < p_c`**; **strict convexity of the inverse correlation length** at the correlation-length scale, uniformly in `p < p_c`; the exploration is a **killed Markov renewal process** | bond FK; two-point function; no amplitude of a component density; no site statement; no matching-diagonal statement |

**Direction versus `p`.** The two analyticities are different and the ticket is right to
separate them. Campanino–Ioffe 2002 supplies regularity in the **direction** `x` — the
asymptotic holds for every direction, so the directional dependence of the decay rate is
not the obstruction. Campanino–Ioffe–Velenik and D'Alimonte–Manolescu supply regularity in
**`p`**: analyticity (CIV) and strict convexity, the latter uniformly in `p < p_c` (DM,
in 2D). For our model `q = 1`, so CIV's Assumption (1.2) is **known**, not conjectural.

Site versus bond: none of the three is stated for site percolation. CIV's skeleton/renewal
machinery is model-agnostic for finite-range independent percolation and its Assumption
(1.2) is exactly the hypothesis one would verify for the site model; but that verification
is not in any of the three, and the ticket's own note that "a negative search is not an
originality certificate" applies to the converse reading too.

---

## 2. The secondary question — the exceptional `d` set — is answered, and the answer is "discrete, not empty"

The earlier round flagged an at-most-countable exceptional set of `d` at which the
finite-median centred fluctuations need not be two independent Gumbels.

**What is now citable.** `kappa_G(p)` is analytic on the whole subcritical interval and, in
2D, strictly convex there, uniformly in `p < p_c`. That is CIV plus D'Alimonte–Manolescu
specialised to `q = 1`. So the possibility that `kappa` fails to be differentiable at some
`p`, which would have wrecked the local inversion used to move the centre by `log(w)/w`, is
closed. `a(d) = kappa_NN^{-1}(d)` and `b(d) = 1 - kappa_matching^{-1}(d)` are therefore
analytic in `d` on the relevant range.

**What that does and does not buy.** Define `h(d) = a(d) + b(d) - 1`. The exceptional `d`
are the zeros of `h`. Analyticity makes the zero set **discrete** — it cannot contain an
interval unless `h` vanishes identically, in which case `a = 1 - b` identically, which is a
strong and separable coincidence. So:

- the exceptional set is at most countable and has no accumulation point inside the valid
  interval, which is what one needs to bound it numerically; and
- it is **not** removed. Removing it means proving `a(d) + b(d) != 1` everywhere on the
  range, and that is a statement about the square-site chain, not about regularity. No
  cited source gives `A`, so none gives `a`, `b`, or `h`.

The honest status of the secondary question is therefore: **analyticity obtained, emptiness
not.** Reported as asked.

---

## 3. The identity that is missing, stated so that it can be attacked

Write the renewal object of the handoff exactly as eq. (5.1) of
`winding-intensity-and-prefactor.md`:

```
L_w = w [z^w y^0] { -log(1 - A(z,y)) }
    = w sum_{n>=1} (1/n) sum_{sum x_i = w, sum y_i = 0} prod_i a(x_i, y_i).
```

The factor `w` is horizontal translation; the `1/n` removes the marked renewal cut. For
this to compute `nu_w` rather than a different object, three separate identities are needed.
They are not three statements of one thing, and only the first has support in the literature.

**(S1) Chain decomposition.** Every complete winding component is a *closed* chain of the
CIV irreducible objects, with the object displacement law having finite mean and variance.
*Support:* CIV §1.2 and §3.1 construct precisely such a representation for the linear
(0-to-x) case, with a local limit theorem for the displacement. The cyclic case is not in
CIV; it is the natural closure of their statement, and it is the part that has to be
supplied.

**(S2) Multiplicity.** Up to the `w` translations, each component carries **exactly one**
admissible marking. *Status: false as it stands, measured below.*

**(S3) Boundary weight.** The external vacant boundary weight `(1-p)^{|dC|}` is absorbed
into the object weights `a(x,y)` as a product over objects. *Status: no source supplies
this.* `|dC|` is not a sum over chain objects: a void site can border two objects, and the
outside of a winding component is itself one connected region whose weight is not
distributed over the chain. This is the point the ticket names and it survives every source
checked.

**Why the two-point OZ amplitude cannot be substituted.** The OZ amplitude of CIV/CI/DM is
the `w^{-1/2}` obtained from the *directional* second derivative of the OZ surface along the
line `0 -> x`, in the geometry of the diamond/tube decomposition. The cylinder prefactor of
(5.1) is a **transverse closure probability on a periodic strip**: it is
`(2 pi D w)^{-1/2}` with `D = sigma^2/mu` from the *transverse* displacement law, and the
periodicity of the transverse direction is what selects coefficient `y^0`. These agree only
if an extra identity identifies the two displacement laws. Copying the number is not
supplying that identity, and the `1/2` remains a hypothesis.

---

## 4. The measured part: (S2) fails, and by how much

I did not argue (S2) — I measured the quantity it is about. `scripts/`-equivalent
`sewing_multiplicity.py` enumerates every occupied configuration of `(Z/wZ) x {0..L-1}`
exactly once with weight `p^{|A|}(1-p)^{wL-|A|}`, finds components with a union-find
carrying integer **lift gains** (so "winds" means the lift has a cycle of nonzero winding),
and for each winding component counts

```
c(C) = the number of distinct rows at which C crosses a fixed reference seam,
```

where "crosses" means containing an edge of nonzero lift gain — for NN that is the
horizontal edge from column `w-1` to column `0`; for NN+NNN it also includes the wrapping
diagonals `(x,w-1)-(x+1,0)` and `(x,0)-(x+1,w-1)`.

At `p = 1/2`, exact rational arithmetic, all lengths `L = 2` to the stated maximum:

| case | lengths | `E[c]` | `P(c >= 2)` | `P(c >= 3)` | `P(c = 0)` |
|---|---|---|---|---|---|
| `w = 2`, NN | 2..9 | **1.4331** | 0.3133 | 0.0896 | 0 |
| `w = 2`, NN+NNN | 2..9 | **2.2337** | 0.6063 | 0.3244 | 0 |
| `w = 3`, NN | 2..6 | **1.5236** | 0.3982 | 0.1039 | 0 |
| `w = 3`, NN+NNN | 2..6 | **2.1667** | 0.6510 | 0.3309 | 0 |
| `w = 4`, NN | 2..5 | **1.5157** | 0.4064 | 0.0954 | 0 |
| `w = 4`, NN+NNN | 2..5 | **2.0625** | 0.6492 | 0.2996 | 0 |

What this says, and what it does not:

1. `P(c = 0) = 0` in every case, as it must be: a winding component crosses the reference
   seam at least once. That is the sanity check on the winding detection.
2. `E[c]` is **strictly greater than 1** and does **not** decay towards 1 over the widths
   where it can be measured: NN sits at `1.43, 1.52, 1.52` and NN+NNN at `2.23, 2.17, 2.06`
   for `w = 2, 3, 4`. The matching adjacency is worse throughout, which is expected — its
   wrapping diagonals give a second route across the seam.
3. Therefore a complete winding component is **not** described by one cut. Whatever marking
   a proof of (S2) chooses, it has to produce the factor `E[c]` or show that its marking is
   not this geometric one. The renewal object of (5.1) contains no such factor.
4. **Consequence for `A`, and only for `A`.** `E[c]` is `O(1)` at the measured widths. An
   `O(1)` multiplicity cannot change the power `w^{-1/2}` or the rate `kappa`; it changes the
   amplitude. So the honest statement is: even if (S1) and (S3) were supplied tomorrow and
   the renewal-loop calculation applied verbatim, the amplitude would be
   `A_renewal * E[c]`-corrected, with `E[c]` not computed here beyond `w = 4`. Whether
   `E[c]` tends to a constant, or grows, is exactly the part of `A` that remains unknown.

Caveat, stated so the number is not overread: these are finite boxes, so a component that
would be cut by the top or bottom boundary in the infinite cylinder is included. That
inflates `E[c]` for small `L`. The quantity that matters for the asymptotic is its `w -> inf`
limit at fixed subcritical `p`, which is not measured here and is not claimed.

---

## 5. Status, per the ticket's own three options

| option | outcome |
|---|---|
| (a) a cited theorem with a complete model/closure mapping | **partial.** Analyticity and strict convexity of `kappa_G(p)`, uniformly in `p < p_c`, are cited with hypotheses and mapped onto this model at `q = 1`. That closes the regularity half. The closure mapping itself is not supplied by any source. |
| (b) a derived sewing lemma with proof and the resulting `A`, `beta` | **not achieved.** (S1) has literature support for the linear case only; (S2) is false as stated and quantified above; (S3) has no source. |
| (c) a precise obstruction | **achieved.** Three named obstructions, one measured exactly. |

**What remains unknown, itemised.**

- **The power.** `1/2` stays a hypothesis. Nothing retrieved gives a component-density
  power, and the two-point OZ amplitude does not transfer without an extra identity.
- **The amplitude `A`.** Unknown, and now known to include a seam-multiplicity factor that
  the renewal-loop calculation does not contain; measured as `E[c] = 1.5` (NN) and `2.1`
  (NN+NNN) at `p = 1/2`, `w <= 4`, with the `w -> inf` limit open.
- **The remainder class.** The renewal-loop calculation gives `1 + O(1/w)`
  **conditionally** on the sewing identity. Unconditionally, nothing beyond `o(1)` is
  available for the actual site clusters.
- **`p`-regularity.** Obtained by citation: analytic on `(0, p_c)`, strictly convex in 2D
  uniformly in `p < p_c`. What is *not* obtained is the emptiness of the zero set of
  `a(d) + b(d) - 1`.

Near-critical crossover is untouched and no claim is made about it: everything above is at
fixed subcritical `p` and the uniformity statements are the cited ones, not new ones.

---

## 6. What was executed

`sewing_multiplicity.py` (new), `sewing_multiplicity.json` (its output). Exact integer and
rational arithmetic throughout; `2^wL` configurations per cell. No Monte Carlo, no GPU, no
width census, no new `p_c`.

No merge, no `docs/STATUS.md` edit, no new issue. Full repository CI has not been run.
