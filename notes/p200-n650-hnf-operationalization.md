# N650 HNF operationalization after the Gaussian CRT no-go

The Gaussian CRT result proves that the proposed antisymmetric path mark is a
configurationwise null.  This extension checks that conclusion in the actual
HNF label convention used by the integer-period Monte Carlo backend, for both
final designs `(23,11)` and `(17,19)`.

For every one of the 650 final labels it computes the exact quotient labels
on N130, N325, and N65.  The maps have uniform fiber sizes `5`, `2`, and `10`.
Both compositions to N65 agree entry by entry.  Within each N65 fiber the
pair `(N130 label,N325 label)` distinguishes all ten final labels, giving the
concrete CRT square rather than a merely abstract group isomorphism.

The final column-HNFs in the C++ runner convention are

```text
(23,11): [[650,593],[0,1]]
(17,19): [[650,343],[0,1]].
```

The artifact freezes SHA-256 hashes of all projection arrays in final
HNF-label order and includes their first sixteen rows as a readable fixture.

## Production decision

No N650 C++ path-flag runner is added.  The endpoint configuration has two
different, well-typed quotient labels; subtracting those integers is
meaningless, while applying the two honest connectivity joins in opposite
orders is exactly zero.

The surviving nonlinear candidate is the symmetric mixed join

\[
\Delta_{25}h=h(\Pi\vee R_2\vee R_5)-h(\Pi\vee R_2)
              -h(\Pi\vee R_5)+h(\Pi).
\]

Its algebra is canonical once the typed topology functional `h` is fixed.
The current threshold-rank histogram is not sufficient: it discards the full
connectivity partition and lifted homology state before the two relations can
be applied.  A matching-charge implementation has one further sharp issue:
ordinary binary pushdown does not preserve black/white complementarity, so
`q` needs a declared two-colour typed transport in a common ambient basis.

The next implementation should freeze `h` first, then emit the four
same-configuration rows `(h0,h2,h5,h25)` and their full `4x4` batch
covariance.  The toy partition-rank value `-4` proves this symmetric channel
is not algebraically forced to vanish; it is not evidence that physical
ambient-H1 or matching `q` couples to it.

## Sharp exploratory replacement

The most informative next object is a **typed two-colour mixed homology
defect**, not an ordinary quotient mask.  Keep black-NN and white-matching
connectivity as separate layers on the final lift.  For each of
`empty,R2,R5,R2 join R5`, add only same-colour fiber-identification edges,
with the exact deck displacement carried by the HNF map, and set

\[
h_R={1\over2}\left[\operatorname{rank}H_1(B\vee R)
                         -\operatorname{rank}H_1(W\vee R)\right].
\]

Then score `Delta25 h`.  Complement swaps the two typed layers, so this
candidate is exactly matching-odd even though a mixed fibre can appear in
both layers after quotienting.  It is deliberately **not** called the usual
matching charge of a binary quotient mask.  Before production it needs only
a tiny exact complement/winding oracle; if that gate passes, the HNF maps in
this artifact already specify every identification edge for N650.

Reproduce:

```bash
python3 scripts/p200_n650_hnf_maps.py \
  --output results/exact-cover-character-oracles/n650_hnf_maps.json
python3 -m unittest discover -s tests -p 'test_p200_n650_hnf_maps.py'
```
