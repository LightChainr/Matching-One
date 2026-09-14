# Perron score formula for the fixed-width charge root

2026-09-14.  Exact finite-state refinement of `fixed-width-charge-free-energy-20260914.md`.  Once the safe no-horizontal-homology transfer matrix is built, both the infinite-length charge coexistence point and its slope can be read from Perron data without numerical differentiation of exponentially small torus probabilities.

## 1. Safe transfer matrix with row-resolved weights

Fix circumference `w` and graph `G` (NN or matching).  Let

\[
T^0_{G,w}(p)                                                  \tag{1.1}
\]

be the finite nonnegative transfer matrix on safe frontier states, rejecting a transition as soon as horizontal homology is created.

Resolve each matrix entry by the newly added row mask `M`:

\[
T^0_{ij}(p)
=\sum_M A_{ij}(M)
 p^{k(M)}(1-p)^{w-k(M)},                                     \tag{1.2}
\]

where `A_ij(M)` is `0/1` (or a finite multiplicity if the state representation intentionally aggregates equivalent microscopic transitions) and

\[
k(M)=\text{number of occupied sites in the added row}.       \tag{1.3}
\]

Let

\[
\lambda^0_{G,w}(p)                                            \tag{1.4}
\]

be the Perron root of the primitive reachable/co-reachable safe block, with positive left/right eigenvectors `l,r` normalized by

\[
l^Tr=1.                                                       \tag{1.5}
\]

The no-homology free energy is

\[
I^0_{G,w}(p)=-\log\lambda^0_{G,w}(p).                         \tag{1.6}
\]

## 2. Perron edge/mask measure

Define a probability measure on safe one-row transitions by

\[
\boxed{
\mathbb Q^0_{G,w,p}(i,M,j)
=\frac{l_i A_{ij}(M)
 p^{k(M)}(1-p)^{w-k(M)}r_j}{\lambda^0}.}                     \tag{2.1}
\]

Normalization follows from

\[
l^TTr=\lambda^0 l^Tr=\lambda^0.                              \tag{2.2}
\]

This is the stationary one-step law of the Perron/Doob process conditioned to survive indefinitely in the safe homology sector, with the microscopic added row retained as a mark.

Write

\[
\bar K^0_{G,w}(p)
=E_{\mathbb Q^0}[k(M)].                                      \tag{2.3}
\]

This is **not** the ordinary Bernoulli mean `wp`; it is the row occupation under the quasi-stationary no-winding phase.

## 3. Exact derivative of the void free energy

Differentiate the Perron root:

\[
(\lambda^0)'=l^T(T^0)'r.                                     \tag{3.1}
\]

For one row mask,

\[
\partial_p\log[p^k(1-p)^{w-k}]
=\frac{k}{p}-\frac{w-k}{1-p}
=\frac{k-wp}{p(1-p)}.                                        \tag{3.2}
\]

Therefore

\[
\boxed{
\partial_p\log\lambda^0_{G,w}(p)
=\frac{\bar K^0_{G,w}(p)-wp}{p(1-p)}.}                       \tag{3.3}
\]

Since `I^0=-log lambda^0`,

\[
\boxed{
(I^0_{G,w})'(p)
=\frac{wp-\bar K^0_{G,w}(p)}{p(1-p)}.}                       \tag{3.4}
\]

Thus strict increase of the no-homology free energy is equivalent to the intuitive quasi-stationary depletion inequality

\[
\boxed{\bar K^0_{G,w}(p)<wp.}                                \tag{3.5}
\]

Avoiding horizontal winding biases the surviving transfer phase toward fewer occupied sites.

## 4. Exact slope of the charge free energy

For black NN density `p`, put

\[
q=1-p.                                                        \tag{4.1}
\]

Recall

\[
\Theta_w(p)
=I^0_{4,w}(p)-I^0_{8,w}(q).                                  \tag{4.2}
\]

Differentiate with respect to `p`:

\[
\Theta_w'(p)
=(I^0_{4,w})'(p)+(I^0_{8,w})'(q).                            \tag{4.3}
\]

Apply (3.4) to the two safe processes.  Since `p(1-p)=pq`,

\[
\boxed{
\Theta_w'(p)
=\frac{
wp-\bar K^0_{4,w}(p)
+wq-\bar K^0_{8,w}(q)
}{pq}.}                                                       \tag{4.4}
\]

Equivalently,

\[
\boxed{
\Theta_w'(p)
=\frac{w-ar K^0_{4,w}(p)-\bar K^0_{8,w}(1-p)}{p(1-p)}.}    \tag{4.5}
\]

At the fixed-width charge-coexistence point `p_w^ch`, this is the exact thermodynamic charge-field slope.

The strict monotonicity of `Theta_w` implies the nontrivial paired quasi-stationary inequality

\[
\boxed{
\bar K^0_{4,w}(p)+\bar K^0_{8,w}(1-p)<w.}                    \tag{4.6}
\]

This compares two different safe transfer phases; it is not the trivial identity `wp+w(1-p)=w` of the unconditioned product measures.

## 5. Second derivative / susceptibility formula

Because the safe matrix is finite and analytic in the interior parameter interval, Perron perturbation gives higher derivatives as ordinary additive-functional cumulants of the Doob chain.

Let

\[
S(M)=\frac{k(M)-wp}{pq}.                                      \tag{5.1}
\]

Then

\[
\partial_p\log\lambda^0=E_{\mathbb Q^0}S.                   \tag{5.2}
\]

Differentiating again produces

\[
\partial_p^2\log\lambda^0
=E_{\mathbb Q^0}[\partial_p S]
+\sum_{n\in\mathbb Z}\operatorname{Cov}_{\mathbb Q^0}(S_0,S_n),\tag{5.3}
\]

with the absolutely convergent Green--Kubo covariance series for the finite primitive Markov chain.

Hence the curvature of the fixed-width charge free energy has a direct quasi-stationary susceptibility interpretation.  This is potentially useful for an `m^{-1/2}` or `m^{-1}` crossover analysis if the first derivative becomes small with width.

No such large-`w` curvature scaling is asserted here.

## 6. Perron amplitudes give the finite-m root correction

For a primitive safe transfer block, standard Perron asymptotics give

\[
P^4_{0;w,m}(p)
=A_{4,w}(p)[\lambda^0_{4,w}(p)]^m
[1+O(\gamma_{4,w}(p)^m)],                                    \tag{6.1}
\]

\[
P^8_{0;w,m}(q)
=A_{8,w}(q)[\lambda^0_{8,w}(q)]^m
[1+O(\gamma_{8,w}(q)^m)],                                    \tag{6.2}
\]

on compact parameter intervals where the subleading spectral ratio is uniformly below one.  The amplitudes depend on the torus cut/closure vectors.

Digital Alexander gives

\[
P^4_2(p)=P^8_0(q).                                            \tag{6.3}
\]

Thus the finite charge fugacity has expansion

\[
\boxed{
\theta_{w,m}(p)
=m\Theta_w(p)+B_w(p)+O(\gamma_w^m),}                          \tag{6.4}
\]

where

\[
B_w(p)=\log A_{8,w}(1-p)-\log A_{4,w}(p).                    \tag{6.5}
\]

Let `p_w^ch` solve `Theta_w=0` and suppose `Theta_w'(p_w^ch)>0`.  The finite torus root `p^*_{w,m}` therefore satisfies

\[
\boxed{
p^*_{w,m}
=p_w^{ch}
-\frac{B_w(p_w^{ch})}{m\Theta_w'(p_w^{ch})}
+O(m^{-2})+O(\gamma_w^m/m),}                                 \tag{6.6}
\]

provided the amplitudes and Perron data are twice differentiable in a neighbourhood of the coexistence point.

This gives a controlled alternative to fitting root drift directly in extremely small endpoint probabilities.

## 7. Why the finite root CDF can be almost flat

At a finite matching root,

\[
P_0=P_2=\chi/2,                                               \tag{7.1}
\]

and the exact charge-coordinate identity gives

\[
F'(p^*)=\frac14\chi(p^*)\theta'(p^*).                         \tag{7.2}
\]

For fixed width and large `m`, (6.1)--(6.4) imply

\[
\chi(p^*)
=e^{-mI_w^*+O(1)},                                            \tag{7.3}
\]

where

\[
I_w^*=I^0_{4,w}(p_w^{ch})=I^0_{8,w}(1-p_w^{ch})>0,            \tag{7.4}
\]

while

\[
\theta'(p^*)
=m\Theta_w'(p_w^{ch})+O(1).                                  \tag{7.5}
\]

Therefore

\[
\boxed{
F'(p^*_{w,m})
=\frac{m\Theta_w'(p_w^{ch})}{4}
 e^{-mI_w^*+O(1)}.}                                           \tag{7.6}
\]

The birth CDF is exponentially flat at the balance root as `m->infinity`, even though the root location itself converges at order `1/m` to the unique free-energy crossing.

This is the quantitative fixed-width version of **balance without concentration**.

## 8. Practical computation

A fixed-width implementation no longer needs to estimate tiny `P0,P2` directly.

For each graph:

1. build the safe no-horizontal-homology transfer block;
2. compute its Perron root `lambda^0(p)` and left/right eigenvectors;
3. solve
   \[
   \lambda^0_{4,w}(p)=\lambda^0_{8,w}(1-p)                   \tag{8.1}
   \]
   for `p_w^ch`;
4. compute `Theta_w'` from the row-occupation formula (4.5), not finite differences;
5. if desired, compute the cut/closure Perron amplitudes for the `1/m` correction (6.6).

This is a potentially high-precision new `p_c` approach through the iterated limit `m->infinity` then `w->infinity`, but no production estimate is claimed until the safe transfer implementation is independently cross-checked.

## 9. Boundary

Equations (2.1)--(4.5) are exact finite-state Perron identities once the safe primitive block is correctly constructed.  The `1/m` correction uses a standard spectral-gap expansion and should be stated with the usual compact-interval simplicity/gap assumptions.  No large-width convergence rate for `p_w^ch->p_c` is asserted.