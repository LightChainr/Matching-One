# P398: weak-source amplification, not a larger bare feedback

Parent: `39e06607ec3a353b1130acebf770da591acaf340`.
This decomposes its unchanged width-eight memory kernels. The 6.81225 ratio
remains correct **per unit source variance**, but its geometric interpretation
is now sharper: the plus source is near-dark, not more strongly coupled in
absolute covariance units.

## The normalization explains the enhancement

| Quantity | psi-minus | psi-plus |
|---|---:|---:|
| Source variance C | 3.295348526 | .2641691831 |
| Bare initial feedback C k(0) | 1.487710352 | .8124372905 |
| Bare integrated feedback C integral(k) | .3994196956 | .2076301378 |
| Bare first moment C integral(t k) | .1086992449 | .0562790772 |

The plus bare initial feedback is only **.5460991041** of the minus value,
while the inverse source variance is larger by **12.47438663**. Their product
is the previously measured **6.812251365** normalized-curvature ratio.

In A/L coordinates the relative off-diagonal coherence is .85157024 for
the source covariance, but only .29357814 for initial feedback, .31593709
for integrated feedback and .31773973 for the first moment. The fixed plus
combination cancels about 85.16% of the diagonal source variance but only
29.36% of the diagonal initial feedback. This unequal cancellation, not
extra raw interaction strength, amplifies the response of the weak source.

## RR/RT/TR/TT: specify the projection and the left force

Let P be the stationary-L2 projection onto **both** A and L, Q=I-P, and
D=Q M Q on its range, where M=-G. The right hidden source columns are
`(QMA,QML)=-(QR,QT2)`. On the left use the adjoint source forces
`(QM^dagger A,QM^dagger L)`, not the same right forces: the process is
nonreversible. We call these R-left/T-left and R-right/T-right. Thus

\[
W_n=F^\dagger\Pi M Q D^{-n}Q M F,\qquad n=0,1,2,
\]

where D^0 is the identity and n=0,1,2 give k(0), integral(k), integral(t k),
respectively. For `v_s=(1,s exp(-i*pi/4))/sqrt(2)`, each reported entry is
`conjugate(v_s[a]) W_n[a,b] v_s[b]/Var(psi_s)`. This is the before-ray
decomposition; all four signed entries sum to the scalar memory moment.

| Ray | Moment | RR = TT | RT = TR | Total |
|---|---|---:|---:|---:|
| minus | k(0) | .1744995730 | .05122925969 | .4514576654 |
| plus | k(0) | 2.176775141 | -.6390535923 | 3.075443097 |
| minus | integral(k) | .04605353793 | .01455002076 | .1212071174 |
| plus | integral(k) | .5744896380 | -.1815025845 | .7859741072 |
| minus | integral(t k) | .01251599951 | .003976830318 | .03298565966 |
| plus | integral(t k) | .1561294170 | -.04960851896 | .2130417961 |

There is no basis for saying that R alone beats T2. Kreweras complement
gives `R(Kx)=-i T2(x)` and `T2(Kx)=R(x)`. If each individual term is first
projected into one protected ray, `P_s R` and `P_s[s exp(-i*pi/4)T2]` are
identical. In that convention all four force contributions become exactly
one quarter of the same kernel. The before-ray table retains the cancelling
opposite-ray components and must not be mistaken for four independent
mechanisms or evidence blocks.

## A physical generator split: contact joins and pair detachments

Write M=M_J+M_D, with J the adjacent joins and D the one-site detachments.
The right hidden forces obey

```text
Q M_D A = 0;   Q M_J L = 0;
Q M_J A = -Q R;   Q M_D L = -Q T2.
```

This makes the elementary contact/size-two origin literal. Splitting both
left and right generator insertions gives signed JJ/JD/DJ/DD terms. At t=0:

```text
minus: JJ=DD=-.03802135622; JD=DJ=.2637501889
plus:  JJ=DD= .1503233989;  JD=DJ=1.387398150.
```

The mixed join-detach terms account for 90.22% of plus k(0), 82.33% of its
integrated memory, and 75.69% of its first moment. For minus they exceed the
net total because the same-generator terms subtract. These are algebraic
response contributions, not literal event probabilities or a demonstration
of path-order memory.

## Which directed geometry distinguishes the rays?

For each actual state transition x->y, assign the signed contribution

\[
\frac{\pi_x\overline{\psi_s(x)}[v_n(x)-v_n(y)]}{C_s},
\quad v_n=QD^{-n}QM\psi_s,\quad n=0,1,2.
\]

Summing this exactly specified budget reproduces the three memory moments.
The classes below are disjoint. Join contact multiplicity counts adjacent
edges between the two merging frontier blocks; detach size is the original
block size, not a newly fitted covariate.

| Directed class | minus k(0) | plus k(0) | minus integral(k) | plus integral(k) |
|---|---:|---:|---:|---:|
| Detach from size 2 | .137321417 | .451383738 | .033330334 | .166017779 |
| Detach from size 3 | -.199179882 | .089227995 | -.029699922 | .039935857 |
| Detach from size >=4 | .287587297 | .997109815 | .056973147 | .187033417 |
| Join at a single contact | .088407416 | 1.086337810 | .027273225 | .226969274 |
| Join with multiple contacts | .137321417 | .451383738 | .033330334 | .166017779 |

The equality of size-two detachment and multiple-contact join budgets is
not independent corroboration: the exact Kreweras event map pairs them.
The finite enumeration maps all 3,432 two-contact join site-events to
size-two detachments; one-contact joins map to sizes 3 through 8.

**The nontrivial discriminator is the embedding of triplets.** Their net
contribution changes sign. For minus, detaching a member of a size-three
block with exactly one same-block nearest neighbor contributes -.21634034
to k(0). The corresponding plus contribution is only +.00705328. Plus
instead has a +.25496021 contribution from detaching a member with two
same-block nearest neighbors, partly cancelled by -.17278549 from zero
such neighbors. The integrated-memory entries retain this difference:
one-neighbor triplet detachment is -.04089546 versus +.01707276.
Hence triplet boundary incidence is an identifiable next geometry mark;
R versus T2 inside a protected ray is not.

These signed budgets are not causal effects of deleting a transition class:
such a rate change would also change the stationary law and projection.

## The next level is an explicit coagulation/singleton-chipping hierarchy

Let `S_C=sum_{j in C} i^j` and
`T_m=sum_{blocks C of size m} S_C`. In particular T1=L. The full-circle
character gives **sum_m T_m=0**, with no additional m factor.
For m>=2, detaching each of the m members removes a size-m block, while
detaching one of m+1 members leaves total charged weight m S_C in size m.
Therefore `G_D T_m=m(T_{m+1}-T_m)`.

For every adjacent edge joining distinct blocks A,B, define

\[
\mathcal Q_m=\sum_{edges}(S_A+S_B)\mathbf1\{|A|+|B|=m\},\quad
B_m=\sum_{edges}\left[S_A\mathbf1\{|A|=m\}+S_B\mathbf1\{|B|=m\}\right].
\]

All actual adjacent edges are counted once; repeated contacts between the
same two blocks retain their physical join-rate multiplicity. Thus

\[
\boxed{GT_m=m(T_{m+1}-T_m)+\mathcal Q_m-B_m\quad(m\ge2).}
\]

For m=1, `B1=2T1`, there is no join gain, and singleton chipping gives
`G_D T1=T2-T1`; hence `GT1=T2-3T1` as before. In particular,

\[
\boxed{GT_2=2T_3-2T_2+S_{11}-B_2,}
\]

where S11 is the charge of adjacent singleton-singleton joins, and B2 is
the charge of each size-two block times its cut-boundary edge count. The
script constructs these motifs and the m=2 identity on the existing finite
state set; the general m statement above is a direct algebraic derivation,
not a new testing campaign.

This hierarchy is not closed by size counts alone: pair adjacency and
boundary multiplicity remain. It explains both the success of the first
R/T2 force correction and why it cannot give an exact two-variable closure.
"Chipping" here is the frontier detach operation; it is not deletion of a
bulk cluster articulation site with arbitrary fragmentation of the remainder.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  /Users/lc/python-envs/research-py311/bin/python scripts/p398_width8_memory_motifs.py
```

One deterministic run produces all budgets and checks their additivity;
no new MC, width extension, fitting, server task or broad verification suite.
The interpretation stays within this finite positive frontier process.
