# Exact single-A-vacancy endpoint coefficients

`first.csv` is parent (5,5), corresponding to child (5,0).
`second.csv` is parent (1,7), corresponding to child (4,3).
In each parent, origin A is vacant, the other 24 A vertices occupied, and
all 2^25 free-B configurations are enumerated once on the actual N50 graph.

Rows k=0..25 store count and sums of q, E, Sstar, q*Sstar and E*Sstar.
The free-coordinate degree is **25**, Ktotal=24+k; integer sums already
include multiplicity. Weight them by p^k*(1-p)^(25-k), not an extra binomial.

Sstar=CB_NN+CW_matching+F+Bvac, Bvac=100−4*Ktotal+Bocc.
All fixed A neighbors are active independently of vertex IDs. Faces close
at their last free B activation. No conventional child graph is substituted
for the defective parent. Only the origin defect was computed; saturated
coefficients must come from the preexisting N25 complement dictionary.

The declared finite endpoint question is H_gain: U(s,t)=g(s)U_child(t),
with zero determinant R=U*U_st−U_s*U_t at s=1,t=0. This producer does not
score R or any response; root-owned rational scoring remains separate.

Contract, code/source hashes, compiler and exact enumeration receipts are
in run.json. No Monte Carlo, old-script replay, tests or cloud jobs ran.
