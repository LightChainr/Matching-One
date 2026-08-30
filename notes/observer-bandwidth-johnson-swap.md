# Exact Johnson/slice clock spectrum

On the uniform `k`-subset slice, one swap step chooses an occupied site and an
empty site uniformly and exchanges them. For the generator `L=P_swap-I`, the
degree-`j` slice-harmonic eigenvalue is

`-j(N-j+1)/(k(N-k))`.

The exact oracle builds the full rational 70-by-70 transition matrix for
`N=8,k=4`. For each `j=0,1,2,3,4`, it computes the nullity of the shifted
matrix and recovers multiplicity

`C(N,j)-C(N,j-1)`.

The five multiplicities sum to all 70 slice states. Independently, evaluation
spaces of monomials of degree at most `j` have exact ranks
`C(8,j) = 1,8,28,56,70`; adjoining their images under the swap matrix does not
increase rank. Thus a degree-four observer has no slice-clock modes beyond
`j=4`. The endpoint slices `k=0,N` are handled as one-state absorbing chains,
without dividing by zero.

These rates are finite sampling-algebra eigenvalues. They must not be fitted or
renamed as scaling dimensions, and they are not Bernoulli-noise powers.
