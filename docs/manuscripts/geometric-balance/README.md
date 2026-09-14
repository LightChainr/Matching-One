# Geometric balance: one manuscript, one new lemma

Read [manuscript.md](manuscript.md). It combines the root theorem of #735 with
#736's full-law question, and supplies the missing arbitrary-direction
staircase corridor. The resulting all-period full-law criterion is
`log N / ell -> 0`; balance-root consistency only needs `ell -> infinity`.

The new proof uses axis-aligned NN rectangles and disjoint translates in the
finite quotient group. It does not rotate the physical interaction, need an
ambient-primitive shortest vector, or import a fixed-width continuum limit.
No unrelated source/Jordan calculation is a dependency.

From the repository root:

```sh
python -m unittest discover -s tests -p 'test_oblique_winding_corridor.py' -v
python scripts/oblique_winding_corridor.py --output /tmp/oblique-corridor-new.json
```

The script refuses to overwrite an existing result. The committed result is
`results/research-control-20260913/oblique-corridor-controls.json`.
Python standard library only. Five local tests and 135,168 tiny configurations
were executed; full repository CI was not run. No Monte Carlo was performed.

The manuscript states the external RSW, site-sharpness and matching inputs,
contains the entire root and concentration arguments, and gives a bounded
closest-source comparison. It is an author-supplied proof; no independent
publication acceptance or originality certification is claimed. Existing
proofs, data, frozen designs and research branches are retained unchanged.

## 2026-09-14 integrated structural continuation

The stacked continuation on PR #771 collects consequences that reduce several
previously separate analysis directions without changing the parent #739
acceptance status.

**Start review with [`round2-claim-ledger-20260914.md`](round2-claim-ledger-20260914.md).**
It is the current claim-level index and explicitly separates `EXACT`,
`AUTHOR-PROOF`, `CONDITIONAL`, `CONJECTURE / PROGRAMME`, and external-literature
boundary statements.  It also records the shortcuts that were rejected or
corrected during self-audit.  The thematic list below is an earlier compact
index; the ledger supersedes it when statuses differ or when newer continuation
files are not listed here.

Since this compact index was first written, PR #771 has also added the
full-period first-exit torus bound, arbitrary-shape fixed-`p` homological free
energy, varying-direction centre theorem, all-direction quantitative mass
monotonicity, exact charge-neutral crossover coordinates, fixed-direction
Poisson--Gumbel extension, giant-white mean reward theorem, fixed-width charge
free energy, and a 2026 near-critical OZ/SITE-renewal boundary audit.  Those
results are organized in the claim ledger rather than duplicated exhaustively
below.

### Exact finite topology and birth structure

- [`structural-consequences-20260914.md`](structural-consequences-20260914.md): persistent 4/8 birth reflection; dual-even/odd birth coordinates; exact same-parameter `(rank,K)` count reduction; alternating black/white barrier Palm identities; marked-Poisson transport; complementary Palm score constraints; convex loop/branch frontier; the `D^{-1}=partial_yy tau` consistency relation; and deterministic directional separation in exponential elongation.
- [`finite-reflection-dominance-20260914.md`](finite-reflection-dominance-20260914.md): graph inclusion plus persistent Alexander duality gives the exact finite stochastic order `1-T2 <=st T1`, hence `P2(1-p)<=P0(p)`, `M(p)+M(1-p)<=0`, `F(p)+F(1-p)<=1`, `Q(u)+Q(1-u)>=1`, and `int_0^1 M<=0` on every honest torus.
- [`rare-charge-balance-20260914.md`](rare-charge-balance-20260914.md): exact finite factorization of the matching-root slope into rare topological-charge susceptibility and conditional endpoint-odds slope, making `balance without concentration` algebraically explicit.
- [`neutral-gas-topological-charge-20260914.md`](neutral-gas-topological-charge-20260914.md): extensive winding counts live in one neutral `(1,1)` source direction while `D=W4-W8=r-1` is a bounded subextensive charge; includes the area-under-`M` reconstruction of the two sharp birth centres.
- [`projective-homology-gas-20260914.md`](projective-homology-gas-20260914.md): exact rank-one state `(slope,K)` and the resulting projective hard-core-gas language for multi-direction competition.

### Alternating black/white geometry

- [`poisson-tessellation-consequence-20260914.md`](poisson-tessellation-consequence-20260914.md): explicit total-variation/Markov-kernel derivation of the Exp component-Palm gap, Gamma(2,1) stationary-location gap, uniform relative position, and higher Gamma spacings. The coarse two-colour limit is a Poisson interval tessellation, not two independent Poisson clouds.
- [`supercritical-white-slab-bulk-20260914.md`](supercritical-white-slab-bulk-20260914.md): exact infinite-cluster boundary-density identity `beta(q)=(1-q)theta(q)/q` and a concrete bulk LLN/CLT target for the huge complementary white component.
- [`dual-surface-excess-slope-20260914.md`](dual-surface-excess-slope-20260914.md): exact finite dual identity `E_black[q n-p b]=E_white[q B-p N]=partial_z log nu`; at regular points the giant-white `O(exp(kappa w))` bulk terms cancel, leaving the `O(w)` excess `pq v w` that directly measures the Gumbel mass slope `v=-kappa'`.

### Directional mass, first-exit bodies and large-d centres

- [`matching-enhancement-mass-gap-20260914.md`](matching-enhancement-mass-gap-20260914.md): author-level two-terminal enhancement proof, using Grimmett--Li plus the corrected square-lattice rerouting theorem of Balister--Bollobas--Riordan, giving `kappa_8(p)<kappa_4(p)` for `p<p_c(G8)` and hence the strict centre inequality `a(d)+b(d)>1`.
- [`directional-enhancement-sandwich-20260914.md`](directional-enhancement-sandwich-20260914.md): one positive `p` sprinkling is beaten by full matching enhancement uniformly over endpoint direction, yielding `tau_8,p(e)<=tau_4,p+delta(e)` for all directions on compact parameter intervals.
- [`vector-first-exit-domain-20260914.md`](vector-first-exit-domain-20260914.md): repeated first-exit skeleton + site BK proves `B_S(t)<1` implies finite exponential susceptibility; certified domains are convex, certificates from different boxes may be convex-hulled, their large-box limit exhausts compact interiors of the Wulff/exponential-moment domain, and the note separates support functions from radial intercepts.
- [`dilute-directional-mass-centres-20260914.md`](dilute-directional-mass-centres-20260914.md): elementary axial dilute bounds `kappa_4=-log p+O(p)`, `kappa_8=-log(3p)+O(p)` and sharp `d->infinity` centre asymptotics.
- [`dilute-directional-geodesic-entropy-20260914.md`](dilute-directional-geodesic-entropy-20260914.md): fixed rational directions satisfy “graph-distance cost minus geodesic entropy”; gives explicit large-`d` tilted-centre asymptotics, e.g. diagonal `p4~(1/2)e^{-d/sqrt2}` versus `p8~e^{-sqrt2 d}`.
- [`p-regularity-audit-20260914.md`](p-regularity-audit-20260914.md): bond mass p-analyticity is classical, but no direct whole-subcritical square-SITE mass p-analyticity theorem is promoted; the at-most-countable exceptional `d` set in the SITE Gumbel theorem therefore remains.

### Prefactor, sewing and finite-width locality

- [`matrix-sewing-unit-residue-20260914.md`](matrix-sewing-unit-residue-20260914.md): finite-state matrix cyclic-sewing theorem. A simple Perron band with finite Markov memory gives exactly `exp(-kappa w)/sqrt(2 pi D w)` with unit logarithmic residue; finite local memory alone cannot explain an anomalous power or amplitude.
- [`sewing-amplitude-diagnostic-20260914.md`](sewing-amplitude-diagnostic-20260914.md): separates pure log-determinant residue, multiple soft-band multiplicity and genuine insertion/mark amplitudes; gives a falsification table for future `beta,zeta,D` certificates.
- [`periodic-mass-locality-mechanism-20260914.md`](periodic-mass-locality-mechanism-20260914.md): pure transverse periodization preserves the zero Fourier Perron mode exactly; any `gamma_w-kappa` comes from wrap-sensitive pieces/decorations. Exponential decoration locality would imply `gamma_w-kappa=O(e^{-cw})` by Perron perturbation.
- [`loop-branch-linear-tail-tests-20260914.md`](loop-branch-linear-tail-tests-20260914.md): because the variational bulge saturates at `r_*<1/2`, every `A>=1/2` is in the candidate's exactly linear branch regime; in particular #762's A=1,2 test predicts `-(1/w)log[P(L>=2w)/P(L>=w)]->kappa` if the #758 rate is correct.
- [`research-frontier-20260914.md`](research-frontier-20260914.md): the homological-free-energy variational picture, a small-`p` actual-SITE transfer/local-CLT programme for the complete-component `w^{-1/2}` prefactor, and the recommended order for the remaining common-window work.

Finite controls for the persistent reflection are in
`scripts/persistent_alexander_birth_reflection.py`, with committed L=3 and L=4
outputs in `results/geometric-consistency/`. The L=3 run exhausts all 512
configurations and all 362,880 strict site orders; the L=4 run exhausts all
65,536 configurations and checks 20,000 fixed-seed site orders. A separate
zero-cost calculation of the constrained-Poisson/topology constants is in
`scripts/poisson_topology_constants.py`.

These notes deliberately separate deterministic consequences, author-level
proofs requiring independent review, deductions that use existing #739
probability inputs, and conjectural proof programmes. They do not identify a
continuum field, claim a full all-subcritical SITE sewing theorem, or turn the
near-critical crossover into an accepted square-site theorem.
