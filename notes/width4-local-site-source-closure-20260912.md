# Actual site-percolation sources expose information erased by a homogeneous quotient

2026-09-12. Exact finite analysis on the unmodified PR708 rank automaton,
head f782061c1a592ed2f9fd0e9dabaa45f0e54bc4e7, certificate Git blob
50b7297deefe7c50215aea2ed534ca5810461af3. Unlike P398, this is the square-site
Bernoulli model itself. New source experiments are not retroactive changes to
any frozen production or to #275's candidate contract.

## 1. A typed multivariate forward map

Let tau(s,b) be the 509-class deterministic successor for a four-bit row b.
For independent probabilities p_(y,j) at spatial row y and column j, write

    w_y(b) = product_j p_(y,j)^b_j (1-p_(y,j))^(1-b_j),
    K_y[s,t] = sum_(b:tau(s,b)=t) w_y(b).

The initializer weights row zero exactly once; the closure has no extra site
weight. For m>=2,

    M(p) = b_0^T K_1 ... K_(m-1) c,       c(s)=close_rank(s)-1.

PR708's arbitrary deterministic future equivalence guarantees this identity
for all such fields. Its rank-specific shear quotient remains valid because
history weights are not changed when equivalent histories are combined. It
still does not retain a primitive line or local-marked output automatically.

## 2. Exhaustive column-grouping classification at width four

Partition the columns into groups. All columns in a group share a probability;
different groups can vary independently in an open cube, independently at
every row. For a mask b, let k_g(b) be its population in group g. The row
weights are the products p_g^k_g(1-p_g)^(|g|-k_g). These multidegree Bernstein
basis functions are linearly independent. Therefore a strong lumping works
for the whole family exactly when, for each source block and each target
block, the counts of transitions into that target agree coefficientwise for
every k-vector.

Refining the rank-output partition by those count signatures gives the
coarsest common strong lumping of THIS finite 509-state representation. This
is a linear equalities argument, not numerical sampling over parameters.
There are 15 set partitions of four columns, with seven D4 types:

| independently varying column groups | exact refinement | final classes |
|---|---|---:|
| {0,1,2,3} | 3,35,94,94 | 94 |
| {0};{1,2,3} | 3,81,302,303,303 | 303 |
| {0,2};{1,3} | 3,62,179,179 | 179 |
| {0,1};{2,3} | 3,90,262,262 | 262 |
| {0};{1};{2,3} | 3,146,508,509,509 | 509 |
| {0};{2};{1,3} | 3,99,303,303 | 303 |
| {0};{1};{2};{3} | 3,164,509,509 | 509 |

Every class is reachable with a physical history of length >=2. For all seven
families, the resulting partition is EXACTLY the orbit partition of the D4
subgroup fixing every probability group setwise, not just the same cardinality.
The script reconstructs D4 actions by relabeling physical histories and checks
all 8*509*16=65,152 transition conjugacies, all rank readouts, permutations and
group multiplication. All final labels and actions are retained in JSON.
This is a finite verified classification, not a general theorem equating
lumping with automorphism orbits at arbitrary width.

This exposes a concrete information difference. Histories [1,1] and [2,2]
are in one 94-state homogeneous class, have equal occupation 2 and immediate
rank 1. Append one row with probabilities (q,p,p,p). Their final expected-rank
difference is EXACTLY q-p. Indeed, with the two initial singleton rows in
column i, the three-row closed rank is 1{the new mask contains i} plus
1{the new row is full}. The full-row term cancels between the two histories.
The homogeneous quotient therefore loses a local source response already
at first order for these conditional histories.

The minima above concern common strong lumpings of the specified automaton.
They are not minimal positive realizations among all encodings, nonlinear
state dimensions, continuum fields, or sample-complexity statements. Two
independently addressed ADJACENT columns suffice to recover all 509 states;
opposite addressed columns leave a reflection. Those are distinct source
families, not a contradiction.

## 3. A spatial two-source experiment with the original rank readout

On a selected row y, set

    p_(y,0)=p+epsilon_y,   p_(y,2)=p-epsilon_y,
    p_(y,1)=p_(y,3)=p.

This is a real site-probability perturbation, valid whenever its probabilities
lie in [0,1]. Reflection j -> 2-j exchanges the two perturbed columns and sends
epsilon_y to -epsilon_y while leaving the baseline law and r invariant.
A single-row source therefore has identically zero first derivative in M.
For two DISTINCT spatial rows a,b, write

    chi_ab(p) = d_(epsilon_a) d_(epsilon_b) M |_(0,0).

At p=1/2 the base integer row matrix A=16K and source matrix H=16K' have
mask weights 1 and 4(b_0-b_2), respectively. The exact mixed derivative is
obtained by two H insertions, the other rows retaining A, divided by 16^m.
In operator language the first H maps invariant functions to a nontrivial
reflection sector and the second can map back. If P is D4 averaging, PHP=0;
a source law built only by projecting this H into the 94 invariant classes
would incorrectly return zero for the separated mixed response.

For the 4x4 torus, using rows 1 and 2 or rows 1 and 3,

    chi_12(1/2) = 327/1024,
    chi_13(1/2) = 633/2048.

An independent lifted-graph census derives their entire unnormalized Bernstein
coefficient polynomials, of degree N-4. The signed configuration factor is
(r-1)(b_(a,0)-b_(a,2))(b_(b,0)-b_(b,2)); division by p^2(1-p)^2 removes the
two necessarily occupied and two necessarily empty source sites. This verifies
the source transfer calculation independently of the state transitions.

### No small-amplitude extrapolation is necessary for M

For distinct source rows the Bernoulli expectation is of degree at most two in
EACH epsilon: only two sites are changed per row. Thus its odd-odd part is
exactly epsilon_a epsilon_b chi_ab. For ANY nonzero admissible amplitudes h,k,

    [M(h,k)-M(h,-k)-M(-h,k)+M(-h,-k)]/(4 h k) = chi_ab(p).

This is an exact multiaffine polynomial identity, not an O(h^2) approximation.
The code verifies it at two unequal amplitude pairs and at a second baseline
probability. A noisy estimator still has variance/amplitude tradeoffs; exact
finite-amplitude unbiasedness does not mean zero Monte Carlo noise. This is
NOT asserted for finite differences of the roots themselves.

## 4. The finite matching root also feels the two-source response

Let p_*(epsilon_a,epsilon_b) be the local root of M in uniform background p.
Strict positivity of M_p at the unperturbed root gives the implicit function.
Single-source root derivatives vanish there. Hence

    d_a d_b p_* |0 = - chi_ab(p_*(0,0))/M_p(p_*(0,0)).

Independent 4x4 integer polynomials isolate the unperturbed root in a rational
90-bisection interval. Termwise outward rational bounds give M_p>0 and chi>0
on the WHOLE interval, not merely at a decimal evaluation. Diagnostics are

    p_*(0,0) = 0.5906721123310283...
    adjacent rows: chi=0.2875032870596704..., root mixed derivative=-0.0577415189712485...
    separated rows: chi=0.2734211524893206..., root mixed derivative=-0.0549132944707030...

These are finite 4x4 responses. No infinite-volume critical point or operator
identity is inferred. The source-specific diagonal second derivatives may
also contribute to a simultaneous finite perturbation; the numbers above
are the mixed coefficients only.

Even the sign of chi is not fixed by symmetry at other geometries: at baseline
p=1/2, adjacent rows in a 4x12 torus give the exact negative coefficient
-74769655273/4398046511104, whereas row separation two gives the positive
16319187079/2199023255552. These DP-only controls are not new independent data.

## 5. Validation and scope

69,888 physical NN configurations across 4x2,4x3,4x4 were traversed independently
in the universal-cover displacement convention. All ranks agree with the PR708
state table. Component saturation and rank-two-in-one-component properties
were also checked. Every per-occupation rank count is conserved exactly.
These finite checks validate the executable, not the all-width topology proof.

The source profile table and separated mixed responses are new executed finite
analyses; no stochastic production was used to estimate them. #275's two
candidate same-source original-U predictions are still not supplied by an
unmarked-rank response calculation. Do not change that frozen contract merely
because these additional sources are analytically accessible.
