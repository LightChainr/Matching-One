# N65 contact-stage held-out production

This fixed experiment uses the paired N65 Gaussian tori `(8,1)` and `(7,4)`.
The two geometries share one counter-keyed Bernoulli occupation for a common
65-displacement transversal.  `x=(0,0)` is vacant, `z=East(x)` is flipped,
and every other displacement is retained as a source `y`.

The carrier is selected before aggregation by all four conditions: alternating
four-arm collar, rank birth `0->1` or `1->2`, Bell state change, and canonical
`g16` change.  Carrier rows retain the full `stage x contact_mask(1,2,3) x
displacement x k x z-state` sufficient statistics.  Global rows retain the
source moments needed to compute `mu_a` and the pooled-root `beta_y` separately
for every displacement.  No mask-0 carrier is manufactured when it has no
support.

Build and run four shards of a 20M production:

```bash
clang++ -std=c++17 -O3 experiments/p537-contact-stage-n65-mc-20260901/producer.cpp -o /tmp/p537-contact-stage-n65
for shard in 0 1 2 3; do
  /tmp/p537-contact-stage-n65 \
    experiments/p537-landing-matrix-preflight-20260901/kernel.tsv \
    results/p537-contact-stage-n65/shard-${shard}.tsv \
    20000000 "$shard" 4 20260901537 0.5927311266364432 100 \
    frozen-N65-contact-stage &
done
wait
```

Score against the independent P45 100-batch original-U baseline:

```bash
python3 experiments/p537-contact-stage-n65-mc-20260901/score.py \
  --baseline results/server-20260828/P45-root-amplitude/n65.hist.csv \
  --tables results/p537-contact-stage-n65/shard-*.tsv \
  --output results/p537-contact-stage-n65/result.json
```

The scorer solves the common matching root, computes `M_t`, `R`, and each
displacement's pooled `beta_y`, then applies the complete Schur term before
summing the frozen `2 x 3` matrix.  P45 and new-MC delete-one factors are kept
as two independent covariance groups.  The primary collapse is
`single=mask1+mask2` versus `double=mask3`; `Delta=det(L)` and
`theta=Delta/(|L11 L22|+|L12 L21|)` are formed inside every replicate before
the two independent jackknife variances are added.  A 10k smoke is labelled
`SMOKE`; it is only an execution/schema check, not a scientific readout.

After the held-out result is fixed, the declared N25-to-N65 mechanism
fingerprint is reproduced with

```bash
python3 experiments/p537-contact-stage-n65-mc-20260901/score_scale.py \
  --n25 results/p537-one-defect-diagonal-edge/contact-stage-tensor.json \
  --n65 results/p537-contact-stage-n65/result.json \
  --output results/p537-contact-stage-n65/scale-fingerprint.json
```

This second command uses the joint N65 covariance, the exact N25 anchor and
the preregistered exposure only.  It compares the split power fingerprint
`[[3,29/8],[3,3]]` with common-power alternatives; it does not rescore or
extend the random block.
