# One shared spatial perturbation, invisible to both immediate Euler/rank states

The preceding two-own-policy mixture is not a paired-H4 perturbation. Here
one next-label law acts on both geometries simultaneously. This definition
is recorded before its new response is read from the saved trajectories.

At a paired ordered prefix Z, let A_a be the vacant labels that preserve
both current ambient ranks and have joint occupied contact degrees
a=(e_first,e_second). Write pi_a=|A_a|/d. Define

```
L_f(u) = 1[old_rank_f=0] (e_f(u)-c_f(u)),
L_s(u) = 1[old_rank_s=0] (e_s(u)-c_s(u)),
g_plus=(L_f+L_s)/2,   g_minus=(L_f-L_s)/2.
```

Each input g defines its own scalar policy (and the two can be combined):

```
q_t(u|Z) = pi_a exp(t*pi_a*g(u))/sum_{v in A_a} exp(t*pi_a*g(v)),
```

with all labels outside these classes unchanged at 1/d. The remaining
suffix is uniform, and the two geometries share its complete permutation.
Every class retains its exact probability. Therefore, for every finite t,
the full joint distribution of both immediate ranks and both Euler
increments 1-e_f,1-e_s is unchanged. Prefixes with neither rank zero have
zero marks and are unchanged.

For either future observable Y the exact derivative is

```
H(g,Y|Z) = sum_a pi_a^2 Cov(g, E[Y|Z,u] | Z,A_a).
```

The existing iid U,V labels give an unbiased masked half-difference
estimator: retain only pairs in the same A_a, and average
`(g(U)-g(V))*(Ybar_U-Ybar_V)/2`. Each Ybar averages the two saved suffixes.
No estimated small class probability, extra census or new suffix is needed.

Responses include both birth CDFs and their full A=F1+F2-1, E=1-F1+F2
combinations. Orientation output coordinates are

```
S=(Y_f+Y_s)/2,
D=(Y_f-Y_s)/(cos4theta_f-cos4theta_s).
```

The two marks by two outputs form a response matrix. The plus-to-D entry
is the principal new question; it is not inferred from the old positive
own-orientation result. D is a finite two-geometry contrast with H4
normalization, not by itself a proof that the entire response is pure H4.
The raw difference changes sign on relabeling orientations while its
normalizing denominator also changes sign; the normalized D is invariant.

We retain all original prefix batches and separately preserve the five
rank cells containing R0:00,01,02,10,20. Other cells contribute exact zero.
These are the original40k prefixes and existing auxiliary suffixes, not
new independent data. Outputs retain all20 aligned batch vectors per N.

Source: e32a85939279b8574278024d647b56d2d1485247 conditional forks and
959a7fa26677c416b874d272f1ba66523fb38f73 contact marks. No new Monte Carlo,
DP, remote session or package installation is needed.
