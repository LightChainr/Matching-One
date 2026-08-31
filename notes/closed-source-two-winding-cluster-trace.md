# An actual two-winding-cluster torus trace survives the regular endpoint zero

**Result.** On a fixed square torus, two separated occupied essential rows
give an honest closed colour-transfer contribution with a nonzero `[2]`
trace. In the ordered two-cluster space `V tensor V`, `[2]` occurs **once**.
Its exact contribution to the declared projected partition is

\[
\boxed{w_{[2]}(a,Q)=a^{2L}Q^{2L-1/2}\frac{Q(Q-3)}2.}\tag{1}
\]

This includes the existing occupied-edge weight, activity and ambient-rank
projection. The configuration has `q=0,E=0`, so this trace enters the
denominator only. Its finite normalization and thermal derivative give
an explicit route into normalized q/E and their original pooled-root
ratio. This route is not annihilated by the regular endpoint identity.
Noncancellation in the whole angular population is a separate question.

This note is based on root `2690f665`, the projector formula in
`d006f9c1:notes/p262-confluent-potts-projector-tomography.md`, and the
regular-endpoint distinction in
`93206494:notes/global-matching-spin4-selection-rule.md`.
It addresses only the fixed two-winding-cluster witness requested by
`notes/weak-q-paths-and-regular-selection.md` in the overview. There is no
new enumeration, coupling point, source fit or simulation.

## 1. Fix one realizable occupation configuration

Take `Z_L x Z_L`, L>=4, N=L², with horizontal coordinate x. Occupy exactly
the rows `y=0` and `y=2`; leave every other site vacant. These two rows
are not NN adjacent, even across the periodic seam. Each is one essential
NN cycle in the x direction. Consequently

\[
K=2L,\quad B=2L,\quad C_B=2,\quad r=1,\quad F_4=0.
\tag{2}
\]

The white matching graph has two connected vacant bands: the row between
the two occupied rows and the band on the other side. A NN or diagonal
step cannot cross a fully occupied row. Thus `C_W=2`, and independently
the original graph observer is

\[
q=C_B-C_W-(K-B+F_4)=0,\qquad E=q^2=0.
\]

The vacant NN edge count and the fixed source are

\[
B_{\rm vac}=2L^2-6L,\qquad
S_*=C_B+C_W+F_4+B_{\rm vac}=2L^2-6L+4.
\tag{3}
\]

This is one specified subset. No translation/orientation multiplicity is
silently attached to it.

## 2. Its true closure is the identity seam on two ordered colour lines

For integer Q>=4 let V be the Q-dimensional permutation representation
of S_Q. Cut the torus across the x seam. The two occupied rows occupy
different fixed spatial positions, so their colour carrier is the
**ordered** space `V tensor V`. The two disconnected clusters may have
the same colour: equality of colours is not equality of graph components.

The local gas edge weight is Q times equality on occupied endpoints and1
when an endpoint is vacant. One column of this fixed occupation pattern
therefore transports the two colours by

\[
T_{\rm col}=a^2Q^2\,I_{V\otimes V}.
\]

Closing L columns identifies each outgoing row colour with its own
incoming colour, with no row exchange. Before the topological projection,

\[
\operatorname{Tr}T_{\rm col}^{L}
 =a^{2L}Q^{2L}\operatorname{Tr}I_{V\otimes V}
 =a^{2L}Q^{2L+2}.
\tag{4}
\]

Every colouring of this subset has r=1. The stipulated `m^-r=Q^-1/2`
therefore multiplies (4) as a scalar, giving its complete reduced weight
`a^(2L) Q^(2L+3/2)`. This is exactly the occupation-family weight
`a^K Q^(B+C_B-r/2)`. All Q^B factors are retained.

The carrier is closed for this **fixed occupation sequence**. Summing
other sequences can change cluster number; no claim that this is an
invariant block of the entire fluctuating transfer matrix is required.

## 3. S_Q decomposition and the genuine `[2]` coefficient

Write `[]=[Q]`, `[1]=[Q-1,1]`, `[2]=[Q-2,2]` and
`[1,1]=[Q-2,1,1]`. Splitting ordered colour pairs into equal-colour,
symmetric unequal-colour and antisymmetric unequal-colour spaces gives

\[
\begin{aligned}
V_{\rm diagonal}&=[]\oplus[1],\\
V_{\rm distinct,sym}&=[]\oplus[1]\oplus[2],\\
V_{\rm distinct,antisym}&=[1]\oplus[1,1].
\end{aligned}
\]

Therefore

\[
V\otimes V=2[]\oplus3[1]\oplus[2]\oplus[1,1],
\tag{5}
\]
\[
Q^2=2+3(Q-1)+\frac{Q(Q-3)}2+\frac{(Q-1)(Q-2)}2.
\]

The repository's unordered distinct-pair projector P_[2] embeds into the
symmetric unequal-colour summand using normalized vectors
`(|a,b>+|b,a>)/sqrt2`. Its multiplicity is one. Since the physical seam
and the fixed-pattern transfer are identities on this summand,

\[
Q^{-1/2}\operatorname{Tr}(P_{[2]}T_{\rm col}^L)
 =a^{2L}Q^{2L-1/2}\operatorname{Tr}P_{[2]},
\quad \operatorname{Tr}P_{[2]}=Q(Q-3)/2,
\]

which proves (1). No factor2 or quotient by exchange of the spatial rows
is appropriate. Imposing distinct colours and discarding the other
summands would change the physical partition; here they all remain.
For m=2, Q=4, the full closure has16 colour states and the `[2]` trace
has dimension2. It supplies exactly1/8 of this fixed subset's positive
weight at that Q.

In contrast, the regular invariant endpoint ell in the unordered carrier
satisfies `ell P_[2]=0`. A closed trace contracts outgoing and incoming
indices with the identity seam, not with this rank-one invariant endpoint.
Even though the identity operator has the same scalar contraction in some
unprojected examples, replacing its seam by endpoint contractions changes
its P_[2] contraction. The witness gives a concrete instance:
endpoint zero, trace `Q(Q-3)/2`, with actual local lattice weights.

## 4. The closed-source Q path and its finite contracted continuation

Along the already fixed path `a=y Q^(-5/2)`, equation (1) is

\[
w_{[2]}(y,Q)=y^{2L}Q^{-3L-1/2}d_2(Q),
\quad d_2(Q)=Q(Q-3)/2.
\tag{6}
\]

The omitted common factor of the full odds partition is `Q^(N+1/2)`.
Restoring it makes the **total** subset weight
`y^(2L)Q^(N-3L+2)`, as required by (3). Normalized observables do not
depend on this common factor; below we consistently use the reduced
partition `Z=sum_A a^K Q^(B+C_B-r/2)`.

The prescribed diagram/occupation continuation gives

\[
w_{[2]}(y,1)=-y^{2L},\qquad
\partial_Qw_{[2]}(y,Q)|_{Q=1}=3L\,y^{2L}.
\tag{7}
\]

Both are finite contracted expressions. The first is negative, so it is
not a probability or a surviving positive colour-sector dimension at
Q=1. The actual integer-colour decomposition was established at Q>=4;
its trace coefficients, and the specified positive occupation family,
provide the continuation. In particular the singular projector is not
treated as a standalone physical operator at Q=1. Other summands in (5)
must be retained; they restore the total colour weight Q² and its correct
derivative. No equal-dimension continuum collision follows from this.

## 5. An explicit denominator-to-thermal-response route

Use an additive derivative bookkeeping channel with `delta Z=h(y,Q)`
equal to this trace term's Q derivative, while `delta N_q=delta N_E=0`.
The latter zeros are the actual rank weights of the chosen subset, not
an endpoint-selection assertion. For either O=q or E, put `f=h/Z`.
Differentiating the normalized observer gives the exact contributions

\[
\delta\langle O\rangle=-f\langle O\rangle,\qquad
\delta\partial_p\langle O\rangle
 =-f\partial_p\langle O\rangle-f_p\langle O\rangle.
\tag{8}
\]

At the iid Q=1 endpoint, (7) and `Z=(1+y)^N` make this kernel explicit:

\[
f(p,1)=3L\,p^{2L}(1-p)^{N-2L},\qquad
f_p=f\frac{2L-Np}{p(1-p)}.
\tag{9}
\]

It is nonzero for `0<p<1`; its thermal derivative vanishes only at
`p=2/L`. No generic `(Q-1)trP_[2]` factor has been appended. Equation (9)
comes from this configuration's true closure coefficient and the tied
Q path, in the declared reduced normalization. It is one additive term
in the full Q response, not a separately measurable positive sector at
Q=1 or a newly proposed source intervention.
For the original source parameter `t`, `Q=exp(2t)` makes this Q-derivative
contribution twice as large at Q=1; the distinction is retained.

For completeness its route into **original pooled-root U** can be written
without an unproved noncancellation. Let the witness affect geometry a;
write `M=(q_a+q_b)/2`, `Y=(E_a-E_b)/Delta`, `D=M_p`, and take `M(p0)=0`.
The contribution (8) gives

\[
M_s=-f q_a/2,\quad Y_s=-f E_a/\Delta,\qquad
p_s=fq_a/(2D),
\]
\[
M_{ps}=-(f_pq_a+f q_{a,p})/2,\qquad
Y_{ps}=-(f_pE_a+f E_{a,p})/\Delta,
\]
\[
\boxed{\frac{U_s}{A_N}
 =\frac{Y_{ps}+Y_{pp}p_s}{D}
  -\frac{Y_p}{D^2}(M_{ps}+M_{pp}p_s).}\tag{10}
\]

Here s labels the additive derivative channel, not a sampled coupling.
For a contribution in geometry b use its corresponding angular sign.
All derivatives are at the original pooled root. A pooled zero does not
set `q_a=q_b=0`; the root-motion term therefore remains. Even in a
single-geometry root quotient, the common `-f` terms cancel but the
`-f_p E/D` term remains. Normalization-only support consequently does
not supply a general thermal-ratio selection rule.

The attribution of individual derivative terms uses the fixed partition
normalization above. The complete sum, including all irreps and occupation
configurations, is the invariant physical Q response. Formula (10) gives
its exact available channel; it does not prove that the channel survives
all other terms or the actual angular difference.

## Interface closed, interface still open

The proposed regular-endpoint mechanism stays excluded. The stronger
claim that the same exclusion automatically removes every ordinary torus
trace is disproved by this realizable two-winding-cluster closure.
The surviving interface has an explicit local carrier, multiplicity-one
P_[2] coefficient, source/rank weights and normalization-to-U map.

What remains open is a nonzero **whole-population** projection onto a
named continuum contribution. This conditioned transfer is a scalar
identity in colour space; it does not produce distinct spectral branches,
identify H4, establish a logarithmic field, or prove an activated size
law. No such conclusion is inferred from the nonzero trace alone.
