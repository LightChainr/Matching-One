# P537 full-T transport quotient

This scorer consumes only the completed N65 P45 baseline and the four frozen
20M contact-stage shards.  The shards retained the complete `K`-resolved
global canonical-pair moments even though their primary gate selected a local
carrier.  No random counters are replayed or added.

The producer fixed one vacant source endpoint `x` and one thermal site
`z=x+e1`; its global source catalogue contains all `y` except `x,z`.  After
Bernoulli integration over the two states of `z`, the four nearest-neighbour
pair directions are C4-equivalent at every `p`.  For each of
`a,qa,Ea,aS,qaS,EaS`, the omitted `y=z=+e1` column is therefore reconstructed
without a new sample as

```text
F_full = sum_(retained 63 directions) F_d
       + (F_-e1 + F_+e2 + F_-e2)/3.
```

Here `a=g16/(16*N)`, `S=K-Np`, and the fixed-vacancy importance weight for
rest count `k` and `z=i` is

```text
(1-p)/samples * (p/p*)^k * ((1-p)/(1-p*))^(N-2-k)
              * ((1-p) if i=0 else p).
```

The complete logit jet is

```text
T_t = jY_t - R*jM_t - R_t*jM,
J_N = (N^(13/8)/2) * T_t/M_t.
```

The P45 baseline and new source block are independent 100-batch covariance
groups.  Each delete-one baseline replicate resolves the pooled root and
reweights the complete source table at that root.  The `p_ref` comparison is
paired within each group before their variances are added.  The fixed
`p_ref=0.592746050790` is the prescribed square-site reference used by the
existing spatial experiment; it is a transport stress test, not a rigorous
exact-`p_c` enclosure.

Run:

```bash
python3 experiments/p537-full-t-transport-20260901/score.py \
  --baseline results/server-20260828/P45-root-amplitude/n65.hist.csv \
  --tables results/p537-contact-stage-n65/shard-0.tsv \
           results/p537-contact-stage-n65/shard-1.tsv \
           results/p537-contact-stage-n65/shard-2.tsv \
           results/p537-contact-stage-n65/shard-3.tsv \
  --n25-result experiments/p537-landing-matrix-preflight-20260901/result.json \
  --output results/p537-full-t-transport/RESULT.json
```

The command refuses to overwrite an existing result.
