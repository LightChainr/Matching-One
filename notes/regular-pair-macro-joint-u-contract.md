# Regular-pair macro joint U: frozen first-pilot contract

**Decision target.** This pilot asks one bounded question: can the fixed
canonical two-site kernel, after projection to one fixed macroscopic separation
window, transmit through the **complete original moving-root U functional** at
both `N=100` and `N=400`? The executable contract is
[`analysis/regular_pair_macro_joint_u_contract.json`](../analysis/regular_pair_macro_joint_u_contract.json).
This note is the **pre-run frozen contract**, not a current progress claim.
It was subsequently implemented and run locally; the valid pilot stopped
unresolved under its declared rule, before `D17/D21`, with no server or
production continuation. See the [completed report](../results/regular-pair-macro-joint-u/REPORT.md).

The priority is forced by the exact critical-point summability result in
[`eed2190c`](https://github.com/LightChainr/Matching-One/blob/eed2190c04b67084ab5aef5827e00377853a0bca/notes/p337-critical-spatial-summability.md):
the raw canonical pair kernel obeys an absolutely summable bound
`E_pc |g_xy| <= C d^(-2-eta)`. More raw-distance points or a fitted raw spatial
exponent therefore do not address the remaining mechanism. The unresolved gate
is whether the Bernoulli thermal derivative, pooled-root motion, and slope
normalization amplify that summable interaction into original U. Accordingly,
raw `E[g]` is retained as a diagnostic, while the complete thermal/pivotal
functional is the primary readout.

## 1. Geometry and the one immutable window

Use the two homothetic square-torus pairs

| N | k | axis | tilted |
|---:|---:|---:|---:|
| 100 | 2 | `(10,0)` | `(8,6)` |
| 400 | 4 | `(20,0)` | `(16,12)` |

For `(a,b)`, use the repository quotient basis `u=(a,b)`, `v=(-b,a)` and
the existing physical N/E/S/W port convention. Keep
`DeltaCos4=1152/625`, `A_N=N^(13/8)/2`, and canonical
`Kreg=K2+K0` with completion coefficient one.

The only spatial projection is

\[
  \mathcal W_N(x,y)=
  \mathbf 1\!\left\{\frac14\le
       \frac{d_{\mathbb T_N}(x,y)}{\sqrt N}\le\frac25\right\}.       \tag{1}
\]

Both boundaries are included. Distance is the exact Euclidean shortest
physical-torus displacement, never quotient-label or graph distance. If its
squared shortest length is `r2`, the integer membership test is

```text
16*r2 >= N  and  25*r2 <= 4*N.
```

The upper radius is strictly below the cut locus. Enumerate every nonzero
quotient class once, retain `z` and `-z` as distinct ordered displacements, and
save the sorted table and SHA256 before occupations are generated. There is no
angle weight and no distance subsampling. The choice `1/4` to `2/5` removes all
fixed-lattice contacts while leaving a broad positive-area annulus; it is frozen
for both sizes and both directions and cannot be adjusted after readout.

## 2. Fixed source, normalization, and structural split

Extend the delivered signed lookup `g16=16 g` by zero when either endpoint is
occupied. The projected source is

\[
 S_{2,\mathcal W}(A)=\frac1{N^2}
       \sum_{x\ne y}\mathcal W_N(x,y)g_{xy}(A).             \tag{2}
\]

With `A=16` sampled anchors and the full accepted displacement table `D_N`, the
configuration estimator is

\[
 H(A)=\frac1{16 A N}\sum_{x\in\mathrm{anchors}}
                    \sum_{z\in D_N}g16_{x,x+z}(A),          \tag{3}
\]

which is unbiased for (2) against every translation-invariant mark. Also save

\[
 \bar C(A)=\frac1{16 A|D_N|}\sum_{x,z}g16_{x,x+z}(A),
 \qquad H=\frac{|D_N|}{N}\bar C,                            \tag{4}
\]

so that raw `E[g]` and its K-score derivative remain inspectable without
becoming the decision target.

For every pair, count occupied exterior NN components meeting both marked port
groups. Accumulate the same readout as `s=2` and `s>=3`; the exact theorem says
`s<=1` must have `g16=0`. Require

```text
total = s=2 + s>=3
```

configurationwise, in every batch/K cell, and term by term after scoring. This
is a predeclared mechanistic decomposition of the same data, not two independent
votes and not an invitation to add further descriptors.

## 3. One Bernoulli stream and the required p jets

At each size generate 100 batches of 500 iid occupation configurations at the
fixed uint64 Bernoulli reference

```text
p_ref = 10934234699625173385 / 2^64
      = 0.592746050790... .
```

Use `mt19937_64`, the four literal seeds in the JSON contract, and one word per
site in canonical order. For each configuration select 16 distinct anchors
uniformly without replacement using the separate fixed anchor stream. Reuse the
occupation vector and anchor indices across axis and tilted geometries at the
same N; use independent streams between N100 and N400. Anchors and pairs are
within-configuration measurements, never inference units.

Retain integer sums by `(N,batch,geometry,K)` for `1,q,E` and, separately for
`total,s2,sge3`, `H,qH,EH`. Reweight this single p stream by the exact binomial
likelihood ratio. Solve the pooled baseline root

\[
 M(p)=\tfrac12(\langle q\rangle_a+\langle q\rangle_t)=0     \tag{5}
\]

inside the immutable bracket `[0.590,0.596]` by bisection to `1e-12`; do not
widen the bracket or add a second p stream. At the fitted root use

\[
 \ell_1=K/p-(N-K)/(1-p),
\]

\[
 \ell_2=\ell_1^2-K/p^2-(N-K)/(1-p)^2                       \tag{6}
\]

to retain q/E jets through order two and H/qH/EH jets through order one.
For any normalized weighted mean,

\[
 \mu_f'=\operatorname{Cov}_w(f,\ell_1),\qquad
 \mu_f''=\operatorname{Cov}_w(f,\ell_2)
             -2\langle\ell_1\rangle_w
                 \operatorname{Cov}_w(f,\ell_1).           \tag{7}
\]

Save both `E[Cbar]` and `partial_p E[Cbar]` for each geometry and stratum. The
second is the direct diagnostic of thermal/pivotal amplification of the raw
summable kernel; it does not replace the U functional below. Save the
importance-weight ESS for every full and delete-one estimate.

## 4. Primary readout: all four original-U terms

At the common root define

\[
 M=\tfrac12(m_a+m_t),\quad
 Y=(e_a-e_t)/\Delta_4,\quad D=M_p,\quad R=Y_p/D.            \tag{8}
\]

For the fixed source, separately center
`j_q=<qH>-<q><H>` and `j_E=<EH>-<E><H>` in each geometry, then form

\[
 jM=\tfrac12(j_{q,a}+j_{q,t}),\qquad
 jY=(j_{E,a}-j_{E,t})/\Delta_4.                             \tag{9}
\]

The scorer must save all four terms, not merely their sum:

\[
 \begin{aligned}
 J_{2,\mathcal W}={}&
 A_N\frac{jY_p}{D}
 -A_N\frac{Y_{pp}jM}{D^2}
 -A_NR\frac{jM_p}{D}
 +A_NR\frac{M_{pp}jM}{D^2}.                               \tag{10}
 \end{aligned}
\]

They are, in order, direct centered transmission, root motion, source-slope
transmission, and root-slope transmission. Also save the joint root tangent
`-jM/D` and

\[
 T_N=N^2J_{2,\mathcal W}.                                  \tag{11}
\]

Equation (10), not raw `E[g]`, is the pilot's decision variable. It is a fixed
bilocal projection of the original-U response, not the unfiltered homogeneous
coupling derivative, a closed-mark covariance, or a free-energy
susceptibility.

## 5. Batch covariance and validity

Delete one **paired batch index** from both same-N geometries, then recompute the
root, reweighted jets, raw diagnostics, all four terms, all strata, and `T_N`.
Use the 100 delete-one vectors for the full jackknife covariance. The two sizes
are block independent by construction. Report simultaneous familywise-95%
two-sided intervals for the two total `T_N` values with
`t_(0.9875,99)`; all other coordinates are descriptive but retain their joint
covariance.

The run is invalid if any pinned hash or exact window check fails, if any
`s<=1` pair is nonzero, if the stratum sum fails, if a full or delete-one root
is unbracketed, if `D` is not resolved positive, if a fitted root leaves
`[0.5905,0.5955]`, if any ESS falls below 90% of retained configurations, or if
a batch is dropped. An invalid interface is stopped, not repaired by altering
the window, bracket, source, seed, or sampling plan.

## 6. Pilot decision, upgrade gate, and stop rule

The thermal transmission is called resolved only when the simultaneous
intervals for **both** total `T100` and `T400` exclude zero and have the same
sign.

- If either interval contains zero or the signs disagree, report
  `thermal_tail_unresolved`, stop this field-ratio route, and do not fit an
  exponent or evaluate a field comparison.
- If both point estimates have the same nonzero sign and every validity gate
  passes, project the configurations needed at each size for a future 99%
  interval half-width of 25% of `|T_N|`:

  ```text
  n_required = n_pilot *
      [t_(0.995,99) * SE(T_N) / (0.25*abs(T_N))]^2.
  ```

  Round up to a whole 100-batch design. A production contract is eligible only
  if this projection is at most 2,000,000 paired configurations per size.
  Eligibility authorizes writing a separately frozen contract, not launching a
  job.
- If either estimate is exactly zero, signs disagree, or the projected ceiling
  is exceeded, stop the fixed-window route at the declared computational
  ceiling. Do not top up this pilot or tune the window, source, completion,
  anchor count, distances, or strata.

Only after the thermal gate passes may the same two readouts be carried into the
conditional single-field comparisons

\[
 D_{17}=T_{400}-2^{-5/4}T_{100},\qquad
 D_{21}=T_{400}-2^{-13/4}T_{100}.                           \tag{12}
\]

These fixed contrasts distinguish two declared power models under a fixed
macroscopic window, one dominant field, and nonzero size-independent loading.
They neither fit a free exponent nor identify a unique field. The pilot reports
them without accepting or rejecting a field model; a production comparison
would first need a frozen finite-size error rule and held-out dilation.

## 7. First-pilot compute envelope

The immutable pilot contains 50,000 paired occupation configurations at each
size: 100,000 unique configurations and 200,000 geometry evaluations in total.
Each configuration uses 16 anchors and every accepted displacement. The hard
ceiling is 800 million pair-kernel evaluations, one producer compilation, two
size-level paired-geometry producer runs, and one joint scorer. The execution
gate is 1,800 seconds wall time and 4 GiB peak RSS per process. This is designed
for CPU execution and makes no GPU, cloud, or automatic production request.

The fixed window is deliberately broad enough to make this a useful first
transmission measurement, while anchor sampling bounds the cost and same-N
common random numbers preserve the directional covariance. Most importantly,
the contract has a terminal outcome: it either resolves the complete
thermal/root/slope transmission cheaply enough to justify one frozen production
contract, or it stops this canonical macro-window mechanism line without adding
another distance grid.
