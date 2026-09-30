# Directly sample transport versus selection without a path of completion censuses

2026-09-30. Finite sampling identity for the next mechanism calculation.
The [735-step source transmission](geometric-source-window-readout-20260930.md)
is now resolved. Its missing interpretation is not more precision: how
much of passive history memory is already in safe transport, and how much
is introduced by survivor selection?

The earlier [safe-path formula](birth-selection-intervention-20260930.md)
uses W=product(1-c_i/(N-i)). Computing c_i by scanning the whole lattice at
every insertion is unnecessary. The following coupled sampler realizes
the same distinction using local hypothetical-completion queries.

## Coupled path and survival flag

Start at a specified first-birth rank-one configuration with flag I=1.
At each subsequent insertion count:

1. Visit the remaining vacancies in a fresh uniform random order until a
   safe site is encountered. Insert that site in the reference path.
2. If the **first** vacancy examined would have completed rank two, set
   I=0 permanently. Continue the reference path even after I becomes0.
3. If no safe vacancy exists, send the reference path to a cemetery and
   set I=0. Do not loop forever or silently cap attempts.

Only the first examination updates I; later rejections do not constitute
additional natural time steps. The inserted site is uniform among current
safe vacancies. Completing-site rejection is tested before altering the
geometry. Rejected sites remain vacant and may be examined again next step.

With s safe out of m vacancies, for each safe v,

    P(accepted v) = 1/s,
    P(first examination safe, accepted v) = 1/m.

Thus the first-examination flag and accepted site's identity are independent
conditional on the current geometry, with flag probability s/m. Until the
flag dies, the chosen sites therefore reproduce exactly the **killed
original uniform process**; the reference path always follows the safe
process. Conditional on its whole accepted path,

    E[I_end | safe path] = product_i s_i/m_i = W,
    E[I_end*f(A_end)] = E_safe[W*f(A_end)].                 (1)

There is no product of noisy completion-count estimates. I is a Bernoulli
realization of the full survival weight, not an approximation to W.
For a chosen endpoint cell C, condition on the reference endpoint in C:

    E[f | I_end=1,C] - E[f | C]
       = Cov(f,I_end | C)/E[I_end | C].                    (2)

Equation (2) gives selection; the unweighted safe-path means give
entrance/transport. Apply the same formula to the two fixed birth cohorts
and subtract, retaining their separate cell probabilities and normalizers.
The resulting components are correlated, not separate evidence blocks.

## Why this can be an actual bounded calculation

A vacancy's completion status is decided by its at most four occupied
neighbours and their lifted union-find potentials; the current engine
already performs this test inside its completion-count loop. Extracting
that one-site query removes the whole-lattice scan from each safe step.
The all-vacancy pair graph is only needed once, at the chosen endpoint,
if the readout is Y=2e/s. Random-order rejection terminates even at s=0;
for s>0 its expected number of probes is (m+1)/(s+1).

This is a consumer-specific algorithm, not another generic validation
framework. The saved original permutation determines each first-birth
configuration and its count. New safe continuations need a distinct RNG
stream; they are conditional branches of those same archived prefixes,
**not new independent preparations**. Direct rank0-to-rank2 births are
outside a rank-one entrance source and remain separately counted.

The next readout is the already specified endpoint completion creation
Y and its natural versus selection-neutralized early/late contrast. A
nonzero unweighted contrast would reject selection-only for that target;
weighting can add, oppose or dominate it. It would still not uniquely split
birth geometry from subsequent safe evolution. No new descriptor, cutoff,
size ladder or forced statistical significance is required.

This note gives the exact mechanism sampler. It is not a claim that the
coupled L512 continuation has already run; the completed experiment is the
different single-kick735-step response linked above.
