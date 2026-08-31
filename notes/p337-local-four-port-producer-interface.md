# Local four-port producer staged; execution awaits contract and GO

Execution update: both theory gates were accepted and the coordinator gave
explicit GO. The single enumeration per geometry is now complete; see
[raw outputs and receipts](../results/p337-local-four-port-insertion/README.md).
The text below records the original staging interface. No response was scored.

Base: `bea717e826df5a22518774b1725ae7bcbe2cb801`.
Root's subsequent frozen contract: `d7f15e68:analysis/p337_local_pair_insertion_contract.json`.
Producer: `scripts/p337_local_four_port_exact.cpp`.
Only source code and a syntax check are delivered at this stage. No
configuration or response has been read out. The final tensor/contract
gates and the root coordinator's GO are required before enumeration.

## Fixed origin readout

The old quotient representative construction puts the origin at vertex0.
The public rollback `root(v)` query is read-only and is used only after
checking that a port is occupied. Four ports are the physical NN directions
`N=(0,1),S=(0,-1),E=(1,0),W=(-1,0)`, reduced using the existing quotient
map. This is not a seam mark or deck-direction substitution.

All 2^25 configurations are retained. The origin-vacant condition is a
source gate, **not conditioning the sampled population**. If the origin
is vacant and the four ports are occupied, source2 is:

- `-2` for `rootN=rootS != rootE=rootW`;
- `-1` for either `rootN=rootE != rootS=rootW` or
  `rootN=rootW != rootS=rootE`;
- `0` otherwise.

The pairing inequalities enforce precisely two distinct occupied NN
components, each containing two ports. The provisional source requested
by the coordinator is `S=source2/2=-t`. The final theory gate must confirm
that definition; this staged code does not itself certify its tensor origin.

## Interface for the root's future scorer

There will be26 rows K=0,...,25 per geometry, with header

```text
k,count,sum_q,sum_e,sum_s2,sum_qs2,sum_es2
```

Divide the last three columns by2 exactly once to obtain S,qS,ES sums.
The raw moments are integers and include every configuration at that K;
no additional binomial multiplicity belongs in the scorer. q/E retain
the old white-Alexander rank formula. No Sstar/Bvac cross moments are
recomputed. The original iid root/D/U are to be imported by the root's
separate scorer.

For every translation-invariant global observable O, including functions
of q/E/K, transitivity gives `E[O S_origin]=E[O N^-1 sum_x S_x]`.
Thus this origin readout has **site-average**, not site-sum, source units.
There is no per-configuration scan over25 sites and no hidden multiplier25.

After contract/GO the producer takes the unchanged CLI `a b output.csv`.
Allowed geometries remain `(5,0)` and `(4,3)`. Execution and output hashes
will be recorded only after the authorized once-per-geometry enumeration.
