# Integrated source competition is an exact collective-acceleration readout

The source partition at `32270fa2` has a closed integrated form. Conditional
on a rank-one prefix, write d=N-k0, h=H2, and mu=E[T], where T is the number
of subsequent insertions until rank two. Then

\[
\boxed{I_D=\frac{h\mu}{N+1},\qquad
I_G=\frac{d+1-(h+1)\mu}{N+1}.}
\]

D means final insertion in an original-checkpoint singleton gate; G is all
other collective completion. These are the p-integrals of the corresponding
canonical F2 contributions. In particular, this source-integral partition
does not require higher moments of the waiting law once h and mu are known.

## Derivation from the complete physical stopping law

Let S(j)=P(T>j) conditional on the prefix, S(0)=1 and S(d)=0. By monotonicity,
the h original gates remain uninserted whenever the event has not occurred.
The already derived source identity is

\[
P(T=j,V\in D)=S(j-1)\frac{h}{d-j+1}.
\]

For K2=k0+j the integrated binomial-tail kernel is
(d-j+1)/(N+1). It cancels the denominator exactly:

\[
I_D=\sum_{j=1}^d S(j-1)\frac{h}{d-j+1}
       \frac{d-j+1}{N+1}
    =\frac{h}{N+1}\sum_{j=1}^d S(j-1)
    =\frac{h\mu}{N+1}.
\]

The total integral is (d+1-mu)/(N+1); subtracting I_D gives I_G.
No continuum approximation, fit or independent-source assumption is used.

## The collective part measures a coupled stopping-time advance

On the same remaining permutation, let T_D be the first position of any
original gate, ignoring every other way of forming rank two. For h=0 set
T_D=d+1, the after-deck no-direct-event sentinel. For all h>=0,

\[
E[T_D]=\frac{d+1}{h+1},\qquad T\le T_D.
\]

Consequently,

\[
\boxed{I_G=\frac{h+1}{N+1}E[T_D-T].}
\]

The collective integrated loading is an exactly normalized advance of the
physical birth clock relative to the original-gates-only stopping rule.
When a direct gate wins, T=T_D. Collective completion is precisely where
the physical stopping time can precede that reference. The comparison
changes the stopping rule on the same permutation; it is not a claim about
physically deleting sites or a unique causal intervention on the lattice.

For h=0 the first formula gives I_D=0 and the sentinel form gives the full
integral as collective. With no collective advance, mu=(d+1)/(h+1) and
I_G=0. These are algebraic boundary cases, not additional fitted regimes.

## Why a stronger collective mechanism can reduce direct loading

Across two prefixes with the same N,k0,h, a mean-clock change delta_mu obeys

\[
\delta I_D=\frac{h\,\delta\mu}{N+1},\quad
\delta I_G=-\frac{(h+1)\delta\mu}{N+1},\quad
\delta I_{\rm total}=-\frac{\delta\mu}{N+1}.
\]

Thus faster collective completion (delta_mu<0) raises the total integral
while lowering its direct-source part. It pre-empts future original gates.
Opposing source changes are structurally expected even when every source
loading is nonnegative. At fixed h, the source changes cancel down to only
1/(h+1) of the collective change in the total readout. This relation is not
an explanation of every cross-geometry contrast, because h and prefix
membership can also change there.

## Minimal population coordinates and dependency

With R indicating the declared rank-one stratum, a fully solved population
requires only three aggregate coordinates for its integral and source split:

\[
r=E[R],\qquad a=E[R\mu],\qquad b=E[Rh\mu].
\]

Its total/direct/collective integrals are respectively
`[(d+1)r-a]/(N+1)`, `b/(N+1)`, and `[(d+1)r-a-b]/(N+1)`.
The new source information beyond prevalence and the first clock moment is
the H2-weighted first moment b, not a new independent higher-order mode.
For the actual hybrid archive use R times the exact-pair selector for these
classified formulas and retain every original fallback in U. The three
channels then remain the same jointly dependent decomposition already
delivered at `32270fa2`.

At a fixed canonical p the kernel no longer cancels d-j+1. That source
readout still probes the shape of the entire survival law. This distinguishes
an integrated collective-acceleration coordinate from finite-p shape loading,
and supplies a concrete way to avoid counting their algebraic reexpressions
as independent scientific evidence.

Scientific card: exact finite-prefix identity; original-gate/collective
sectors of the R1 contribution to F2; all N and all h under the stated
uniform-permutation law. No new random source, DP, fitting family or numerical
verification run. The same N325/N425 archive and original batches apply when
the identity is used for population loading. Global A_top and other checkpoint
strata are not completed by this identity.
