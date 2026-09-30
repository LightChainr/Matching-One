# A distributed stochastic rank-three code, and its single-readout cost

2026-09-29. Bounded sidecar for #809/#833; consumers are the
[deterministic encoding obstruction and physical product queries](deterministic-block-encoding-obstruction-20260929.md).
This is an elementary parity-sharing construction with exact proofs, not a
novelty claim or a new full-machine positive-realization theorem.

**Answer.** Private stochastic encoding does evade the deterministic
four-state obstruction: for every `k>=2`, two independent uniform parity
shares give four disjoint physical preparation distributions, every proper
subset of blocks reveals nothing, and the **full tensor-response rank is
exactly three**. But the explicit physical product-query family has at most
`2^(1-k)` total-variation distinction for one logical bit. Even an arbitrary
single monotone Boolean readout has sharp distinction only

    b_k = binom(k-1, floor((k-1)/2)) / 2^(k-1).

Thus this is a genuine change of encoding hypothesis, not a demonstrated
reliable redundant memory. The explicit product-query family even admits
a three-label normalized positive initial-memory model, matching its linear
rank. Tensor-response rank three is **not** asserted as an upper bound on
the entire physical future-query language.

## 1. Encoding and its private-randomness contract

Logical states are `(u,v) in {0,1}^2`. Independently sample

    U uniformly on {x in {0,1}^k : xor_i x_i = u},
    V uniformly on {y in {0,1}^k : xor_i y_i = v}.

Physical block `i` receives corner `(U_i,V_i)` of the original five-column
preparation. The random samples/encoder seed are not disclosed to the query
controller. They are not additional measurements. A code distribution
contains `2^(2k-2)` equiprobable physical preparations. The four supports are
disjoint and partition all `4^k` preparations; observing the entire sampled
configuration would identify both logical bits exactly.

Every set of fewer than k physical blocks is uniformly distributed on its
four-corner configurations, independently of `(u,v)`. Indeed every proper
coordinate marginal of a uniform parity class is uniform: after fixing r<k
coordinates there are exactly `2^(k-r-1)` completions of either parity.
Independence of U and V then gives the block claim. In particular every
individual block is uninformative for `k>=2`. This is truly distributed
stochastic sharing, not a random name for the one-informative-block
deterministic exception.

The distinctions are between four **mixed preparations**. Disjoint support
does not imply that an allowed observation can identify which support was
chosen, and does not by itself impose four states on every positive model
for a restricted observation language.

## 2. Full tensor expectation and exact rank

For disjoint subsets `A,B` of `{1,...,k}`, the admitted tensor monomial is

    q_(A,B)(U,V) = product_(i in A) U_i * product_(i in B) V_i.

Disjointness expresses the local response `(1,U_i,V_i)`: at most one of the
two coordinates can appear in each block. Write `sigma_u=(-1)^u` and
`sigma_v=(-1)^v`. Uniform parity gives

    E_u[product_(i in A) U_i]
      = 2^(-|A|) [1 + 1{A=all} (-1)^k sigma_u].

For proper A this is its uniform marginal. For A=all, the all-ones vector
has parity k mod 2 and either has probability `2^(1-k)` or is excluded.
Consequently the encoded response is exactly

    E_(u,v) q_(A,B)
      = 2^(-|A|-|B|)
        [1 + 1{A=all} (-1)^k sigma_u]
        [1 + 1{B=all} (-1)^k sigma_v].                 (2.1)

Since A and B are disjoint, both cannot equal all. No `sigma_u sigma_v`
term occurs. Every tensor response lies in the three-dimensional span of
`1,sigma_u,sigma_v`. This is an upper bound for the **declared tensor
language**, including all its linear combinations.

Conversely, the constant query, the all-U product and the all-V product
give respectively

    1,
    2^(-k) [1+(-1)^k sigma_u],
    2^(-k) [1+(-1)^k sigma_v].                         (2.2)

These three columns are independent, so rank is exactly three for every
finite k. In row order `(00,10,01,11)` the only relation is

    response_00 - response_10 - response_01 + response_11 = 0.

The unseen joint logical parity requires the microscopic character
`product_i (-1)^(U_i+V_i)`, which uses both coordinates in every block and
is absent from this tensor space. This explains how disjoint underlying
supports coexist with dependent response rows.

## 3. What the actual product queries read

The prior note's section 5 constructs a **single two-new-row rank experiment**
with deterministic masks for every disjoint A,B, whose final NN occupied
rank on the `w x 5` torus, `w>=5k`, equals `q_(A,B)` pointwise. Averaging that
same experiment over this stochastic preparation yields (2.1). No new
geometric argument or graph enumeration is needed here.

For distinguishing u=0 from u=1 with v fixed, the only informative column
in this product family is `A=all, B=empty`. Its two Bernoulli output laws
have total variation

    TV_product = 2^(1-k),
    optimal equal-prior single-query error = (1-2^(1-k))/2.    (3.1)

All other product queries have identical laws for the two u values.
Choosing a query randomly independently of the logical bit, even while
reporting the query label, cannot improve this maximum. The v statement
is symmetric. Although linear inversion of (2.2) recovers `sigma_u` from
the exact response expectation, that inversion is not a single Boolean
measurement or a reliable single-preparation decoder.

### An exact three-label positive model for this single-query family

There is a stronger operational distinction than rank alone. Put
`a=1{u=k mod 2}`, `b=1{v=k mod 2}`, and `c=2^(1-k)<=1/2`. Encode the
logical mixed preparation into three latent labels with probabilities

    pi_(u,v) = (1-(a+b)/2, a/2, b/2).

For the all-U query use success probabilities `(0,2c,0)` on these labels;
for all-V use `(0,0,2c)`. Every other product query has a logical-state
independent success probability p by (2.1); use `(p,p,p)`. All encoding
and response entries are nonnegative and normalized as probabilities.
The resulting binary output law agrees with every encoded product query.
The latent labels are statistical labels, not asserted physical preparations.
Thus the minimum common normalized **initial** memory for this family of
single experiments is exactly three: the rank gives a lower bound of
three and this explicit construction attains it. Disjoint physical
supports do not preserve the old positive-versus-linear gap here.

This is not a recursive model for growth, nor a model of multiple
correlated measurements of the same sampled configuration. It also does
not treat every bounded linear combination of response expectations as
an additional allowed Boolean measurement: a signed linear combination
can take latent response probabilities outside [0,1]. The construction
serves the actual product-query columns and their randomized choices.

### Fresh independent encodings have an explicit repetition cost

One explicit operational cost: repeat the all-U product on n **fresh,
independently randomized encodings** of the same logical state. Under one
hypothesis all outcomes are zero; under the other each is Bernoulli with
success probability `c=2^(1-k)`. The optimal equal-prior error is exactly
`(1-c)^n/2`. To reach error eta<1/2 requires and suffices

    n >= ceil(log(2 eta) / log(1-c)).                    (3.2)

For fixed eta this is order `2^(k-1)`. Repeating a noiseless measurement on
the same sampled preparation is not a new independent encoding and cannot
be substituted into this calculation.

### The physical full-language boundary remains

The constructed physical product family has rank exactly three after
encoding. A larger actual future-query family has rank **at least three
and at most four**; proving that its columns all lie in the local tensor
space would be an additional upper-bound input. The earlier physical
product construction supplied inclusion, not that upper bound.

For example, the abstract monotone function
`f(U,V)=product_i U_i V_i` would have encoded expectation

    2^(-2k)[1+(-1)^k sigma_u][1+(-1)^k sigma_v],

and would restore rank four together with (2.2). It uses both coordinates
per block, and is **not** asserted to be a realizable future rank query in
the original geometry. It demonstrates why monotonicity alone cannot
prove the full-physical-language rank-three upper bound.

## 4. Sharp single monotone readout bound

Let `mu_0,mu_1` be the uniform even/odd laws on `{0,1}^k`. Let
`f:{0,1}^k -> {0,1}` be coordinatewise nondecreasing, and observe just
`Y=f(U)`. The following result holds for every k>=1:

    TV(Law_mu0(Y),Law_mu1(Y))
      = |E_mu0 f - E_mu1 f| <= b_k,
    b_k = binom(k-1,floor((k-1)/2)) / 2^(k-1).          (4.1)

The constant is attained by a Hamming-weight threshold. Hence the optimal
equal-prior error among **all abstract monotone Boolean readouts** is
`(1-b_k)/2`; the decoder may reverse its binary decision as appropriate.
Attainment in this abstract class is not a claim that the threshold can
be implemented by the allowed physical rank queries.

### Proof by layer averaging

Let `a_j` be the mean of f over vectors of weight j. Couple a uniform
j-subset to a uniform `(j+1)`-subset by adding a uniformly chosen missing
coordinate. Monotonicity implies

    0 <= a_0 <= a_1 <= ... <= a_k <= 1.

Set `d_m=a_m-a_(m-1)>=0`, so `sum_(m=1)^k d_m<=1`. The constant layer term
cancels, and the elementary alternating binomial tail identity gives

    sum_(j=0)^k (-1)^j binom(k,j) a_j
      = sum_(m=1)^k d_m sum_(j=m)^k (-1)^j binom(k,j)
      = sum_(m=1)^k d_m (-1)^m binom(k-1,m-1).         (4.2)

The identity follows directly from Pascal's rule, whose two tails cancel
except for their first term. The absolute value of (4.2) is at most the
largest coefficient `binom(k-1,floor((k-1)/2))`. Dividing by `2^(k-1)`
gives the difference of the two conditional expectations and proves (4.1).

For `f=1{|U|>=m}`, only `d_m=1` is nonzero, and its signed difference is
`(-1)^m binom(k-1,m-1)/2^(k-1)`. Choose m with m-1 a middle index to attain
the bound. No regularity, symmetry of f, or enumeration is assumed.
The same proof works for monotone functions into [0,1] interpreted as
success probabilities of a randomized binary output.

For k=2,3,4 the sharp TVs are respectively `1/2,1/2,3/8`. By the central
binomial estimate, `b_k ~ sqrt(2/(pi*k))`; even the best monotone binary
readout loses all constant advantage as k grows. This polynomial loss
should not be confused with the stronger exponential loss (3.1) for the
explicit product family.

### Conditioning on V and on independent query noise

For a fixed logical v, V is sampled independently of the **entire** U
vector under both u hypotheses. For each fixed V=z, a Boolean readout
`f(U,z)` monotone in U satisfies (4.1). Averaging proves the same upper
bound for Y alone. Even revealing z gives

    TV(Law_0(V,Y),Law_1(V,Y))
      = E_V |E_mu0 f(U,V)-E_mu1 f(U,V)| <= b_k.        (4.3)

The same argument permits independent public query randomization or
independent future-site randomness, conditional on which the output is
monotone in U. V need not itself be observed or have a product distribution
across its coordinates; its fixed-parity law is allowed.

**Important side-information condition:** merely having a marginal law
independent of u is insufficient. The side variable
`W=(U_1,...,U_(k-1))` has the same uniform law under both parities, but
observing W together with the monotone output `Y=U_k` recovers u exactly.
Conditioning here requires independence from the sampled U, as our V and
the independent query random seeds have. Revealing the encoder's private
seed would change the experiment.

### Consequence for a single final physical rank

Fix the appended row masks and V. Adding old U sites cannot decrease the
occupied NN ambient rank on that fixed final torus. Because the old first
row is empty, the rank is Boolean (0 or 1). Thus every **single final-rank,
nonadaptive** query with future-site randomness independent of the old
sample satisfies (4.1), at any fixed chosen height, even if its response
does not lie in the declared tensor span. This applies to the original
deterministic or independently noisy prescribed source experiments.

It does not bound TV of an entire vector of intermediate outputs: several
monotone coordinates together can reveal parity (the identity vector is
a simple abstract example). Nor is it a claim about unrestricted feedback
controllers: fixing their random seed need not leave a monotone final
function of U when the later commands depend on earlier observations.
Multiple preparations and their costs are separate contracts as well.

## 5. Hypothetical noisy access to all shares is not a free physical readout

Suppose, under a stronger observation model, all k U bits can be read with
independent bit errors `E_i~Bernoulli(epsilon)`, `0<=epsilon<=1/2`, and one
observes `Y=U xor E`. Parity decoding then has exact error

    Pr[xor_i Y_i != u] = [1-(1-2epsilon)^k]/2.         (5.1)

Indeed `E[(-1)^(sum E_i)]=(1-2epsilon)^k`, which is the even-minus-odd
error probability. More precisely, putting `a=(1-2epsilon)^k`, the full
noisy-share output law is

    Pr[Y=y | u] = 2^(-k)[1+(-1)^u a (-1)^(|y|)].     (5.2)

Its even/odd TV is exactly a, so parity decoding is Bayes optimal for
equal priors in this stronger model. For any fixed `epsilon>0` it tends
to a fair guess as k grows. Noiseless full-share access gives TV=1; it is
not a monotone Boolean readout.

If only a monotone Boolean function of the noisy shares may be read,
(5.2) scales every difference in section 4 by a. The sharp TV is therefore
`a*b_k`, again attained by a weight threshold in the abstract class.
This is distinct from (5.1), which uses the entire vector and a nonmonotone
parity decoder.

Independent future instruction flips in a physical query are not, without
a separate channel construction, independent BSC errors on measured old
shares. Neither (5.1) nor (5.2) is claimed as that physical channel's law.

## 6. One bounded exact control and interpretation

The [standard-library checker](../analysis/stochastic-parity-block-encoding-20260929/check.py)
performs one run at k=2,3,4. It enumerates the four finite code supports,
all `3^k` tensor monomials, and the monotone Boolean functions generated
recursively from sections `f_0<=f_1`; it does not scan a larger encoder
space or rerun the previous deterministic census. Integer counts and
`Fraction` arithmetic check the moments, rank, private marginals, sharp
TVs, threshold equality, and the separate BSC formulas at epsilon=1/10
and 1/4. All-k statements rely on the displayed proofs.
The one run checked 468 tensor expectations and 194 monotone Boolean
functions, with zero failures; the three computed tensor ranks were all 3.
The three-label positive construction above is a displayed algebraic proof,
not an extra computation or an output of this checker.

[result.json](../analysis/stochastic-parity-block-encoding-20260929/result.json)
records both the computed scope and the missing full physical-language
upper bound. Run from the repository root:

```sh
python3 analysis/stochastic-parity-block-encoding-20260929/check.py
```

The four mixed preparations have exact tensor response dimension three
at width `w>=5k` and two new query rows; no constant single-readout
reliability is gained by increasing k. The constructive lesson is narrow
but useful: **randomness evades the deterministic rank obstruction, while
a single monotone Boolean readout cannot make this particular code a
large-k reliable single-shot replacement.** This is not a general
stochastic-encoder impossibility theorem. It does not settle arbitrary
future-query rank, adaptive transcript decoding, full-language positive memory,
or a different geometry. No cloud work, navigation changes, full-machine
audit, or novelty certification accompanies this sidecar.
