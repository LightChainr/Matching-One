# Full topology loading: clock center, source competition and a usable paired policy

The complete topology integral has an exact mechanism interpretation:

\[
 \boxed{\int_0^1 A_{\rm top}(p)\,dp
 =1-\frac{K_1+K_2}{N+1}=1-\frac{2C}{N+1}},\qquad
 C=(K_1+K_2)/2.
\]

The rank-one plateau integral, by contrast, is `W/(N+1)`, where
`W=K2-K1`. Thus the complete topology integral cancels the lifetime coordinate
W and reads the center of the two births. Large variation of a gated R1
completion term need not survive this cancellation. This identity also
includes direct rank-zero-to-rank-two births, for which K1=K2 and W=0.

## 1. Complete observable and checkpoint layers

For one insertion permutation define

\[
 g_p(k)=P\{\operatorname{Bin}(N,p)\ge k\},\quad
 F_1(p)=g_p(K_1),\quad F_2(p)=g_p(K_2).
\]

Here F1/F2 are permutation-level canonical contributions; their ensemble
means are the canonical rank-at-least-one/rank-two probabilities. The rank
probabilities on this permutation are `(1-F1,F1-F2,F2)`, so
`A_top=F2-(1-F1)=F1+F2-1`. The elementary beta integral gives
`integral g_p(k) dp=(N-k+1)/(N+1)`, proving the boxed identity and
`integral(F1-F2) dp=(K2-K1)/(N+1)`.

At the fixed checkpoint k0, write R=rank(k0). Its layers are

| Layer | Births measurable from the ordered prefix | Future contribution |
|---|---|---|
| R0 | Neither | Both F1 and F2 |
| R1 | K1 | F2 from K2=k0+T |
| R2 | K1 and K2 | None |

The full contribution is the exact sum

\[
 A_{\rm top}(p)=\sum_{r=0}^2 1_{R=r}
 [g_p(K_1)+g_p(K_2)-1].
\]

The subtraction of one belongs to every layer with its indicator. Replacing
the R1 F2 term alone does not complete the R1 topology term, much less global
A_top. Across orientations f,s, with layer frequencies r_i and conditional
topology means a_i, each layer contrast splits exactly as

\[
 r_f a_f-r_s a_s
 =\tfrac12(r_f-r_s)(a_f+a_s)
 +\tfrac12(r_f+r_s)(a_f-a_s).
\]

The same original shared-batch covariance must cover all three layers and
both terms. The old four-state `(1_R1,f,1_R1,s)` variance allocation does not
transfer to complete A_top; the natural global checkpoint partition is the
nine joint states `(R_f,R_s)`.

## 2. Immediate global whole-pair replacement

Let \(\mathcal F\) contain the **common ordered prefix**, both embedded
checkpoint states, and the prefix-only computation policy (including
exogenous runtime information if used). An occupied label set alone does
not make the past K1/K2 measurable. Conditioning on the complete prefix
order does: future labels remain a uniform remaining permutation, and the
existing set-based R1 clock law is still its conditional future law.

Let X be the global vector containing both orientations' F1/F2 at p_ref
and their p-integrals. The following gate uses the already available clocks:

\[
 G=1\{R_f\ge1,\ R_s\ge1,\ \text{every R1 orientation has its
 original whole-pair accepted exact clock}\}.
\]

When both orientations are R2, accept the identity replacement without any
clock solve. Operationally:

1. Any orientation R0: keep the **entire global vector X** unchanged.
2. Both R2: all components are prefix-measurable, so replacement equals X.
3. Otherwise both ranks are at least one: accept only the original
   `exact_pair` clock status; an original whole-pair fallback stays a global
   fallback, even if it saved a partially completed clock.
4. On acceptance, keep every prefix-measurable F1 and R2 F2 exactly as known;
   replace each R1 F2 by its saved full conditional clock average, for both
   canonical and integrated kernels together.

Then set

\[
 Z=G E[X\mid\mathcal F]+(1-G)X.
\]

This is a sufficient immediately executable policy, not a claim of maximum
possible replacement. It needs no new R0 solver and no joint-suffix DP.
Marginal conditional clock means are enough to obtain the whole conditional
mean vector; the shared suffix dependence is retained in the original paired
observations and their residual outer products.

In particular, an old `exact_pair` containing one R1 row and one missing row
must not automatically pass this global gate. The missing row may be R0,
whose global F1/F2 remain suffix-random. Partial substitution would leave
`Cov(replaced residual, unreplaced R0 response)` uncontrolled, so the old
R1-only covariance guarantee cannot simply be relabeled global.

## 3. Noise and cancellation identities

Write `epsilon=X-Z=G(X-E[X|F])`. Since G is prefix-measurable,
`E[epsilon|F]=0`, `E[Z]=E[X]`, and `Cov(Z,epsilon)=0`. Consequently

\[
 \boxed{\operatorname{Cov}(X)=\operatorname{Cov}(Z)
 +E[\epsilon\epsilon^T]}.
\]

Every fixed affine readout inherits this identity: full A_top, plateau,
orientation sum/difference and their original H4 normalization. This is a
population identity. Finite-sample covariance subtraction need not equal
the empirical residual second moment exactly; retain both rather than
forcing an identity in the sample.

On accepted pairs, F1 has no residual because both K1 are past. R2 F2 also
has no residual. Thus `epsilon_A,i=epsilon_F2,i` for R1 and zero for R2.
For the integrated full topology contrast,

\[
 \epsilon_{A_f-A_s}^{\rm int}
 =-\frac{G}{N+1}
 [1_{R_f=1}(T_f-E[T_f|\mathcal F])
  -1_{R_s=1}(T_s-E[T_s|\mathcal F])].
\]

Its squared paired residual gives the removed suffix-noise second moment;
do not delete the cross-orientation term or replace it by independent-clock
variances. The global gate and the denominator of the complete observable
both differ from the old R1-only result, so its previous noise-removal
percentages are not predictions for this global quantity.

For the two-birth orientation contrasts B1 and B2,
`Var(B1+B2)=Var(B1)+Var(B2)+2Cov(B1,B2)`. Both mean cancellation and covariance
cancellation are scientifically relevant. In the integral, `B1+B2` reads
the contrast of 2C and `B1-B2` reads W, up to their exact signs/scales.
Large plateau or layer terms can therefore coexist with a small center
contrast without indicating a missing physical source.

## 4. Completion sources become signed full-topology contributions

For an R1 prefix, retain the original-singleton direct set D with h=H2,
full safe counts f_j and S(j)=f_j/binomial(d,j). The already-derived laws are

\[
 p_D(j)=S(j-1)h/(d-j+1),\qquad
 p_G(j)=S(j-1)-S(j)-p_D(j).
\]

G in the subscript here denotes the collective source, not the replacement
gate. Define for either completion source s
`pi_s=sum_j p_s(j)` and `tau_s=sum_j j p_s(j)`; then
`sum_s pi_s=1` and `sum_s tau_s=E[T|F]`. These are direct coefficient
readouts with no network or reliability recomputation.

The full R1 topology contribution marked by the eventual completion source
is exactly

\[
 A_s(p|\mathcal F)
 =E[A_{\rm top}(p)1_{\text{final source}=s}|\mathcal F]
 =\pi_s[g_p(K_1)-1]+\sum_j p_s(j)g_p(k_0+j).
\]

Its integral is

\[
 \boxed{A_s^{\rm int}(\mathcal F)
 =\frac{(N+1-K_1-k_0)\pi_s-\tau_s}{N+1}}.
\]

Therefore source competition separates a winning-probability term from a
source-weighted delay term, while including the already-known first birth.
Summing gives

\[
 E[A_{\rm top}^{\rm int}|\mathcal F,R1]
 =1-\frac{K_1+k_0+E[T|\mathcal F]}{N+1}
 =1-\frac{2k_0-\text{age}+E[T|\mathcal F]}{N+1}.
\]

This is the direct age-versus-completion-wait competition read by the full
clock center. The earlier positive F2 source loading omits
`pi_s*(g_p(K1)-1)` and cannot determine the sign of this full marked topology
term. The correct source remains defined by the original checkpoint D;
collective includes new gates made during subsequent safe insertions.

Apply these marked conditional means only on globally accepted R1 pairs.
Keep R2's known-past component separately and put all global fallback vectors
in their full unclassified layer. This preserves additivity without asserting
a microscopic source for an unsolved R0 path. Full-topology source/fallback
terms can be signed, so the earlier nonnegative-F2 unknown-allocation
envelopes must not be reused for global A_top.

## 5. First scientific products, without another solver campaign

1. Full F1/F2/A_top mean contrasts and one shared covariance: does the
   completion source survive, reinforce or cancel the first-birth term?
   Read integrated C and W simultaneously to expose center/lifetime mixing.
2. All nine checkpoint-state contributions with their cross-covariances:
   identify which R0/R1/R2 layers create the global mean and noise, without
   importing the previous gated four-state percentage.
3. On accepted R1 pairs, `(pi_D,pi_G,tau_D,tau_G)` and the boxed full marked
   topology contributions: is a source contrast driven by winning probability,
   waiting time, or cancellation against the first birth?

These readouts belong in the same original 20-batch covariance per size,
with gate counts and unclassified mass visible. They require the already
planned replay of original counters for full K1/K2 and the stored exact
clock coefficients; no new MC or DP follows from this note.

Scientific card: this advances from an R1 F2 contribution to the complete
topological observable and gives an executable paired conditional policy.
The structural claims are exact finite-permutation identities. No numerical
global finding, field identity or asymptotic result is asserted yet. Source
lineage is the original e81dd59f paired populations, complete clock archive
0d1e586d, direct/collective decomposition 32270fa2 and shared C/L covariance
6133b39d. Their resulting readouts are dependent, not new independent streams.
