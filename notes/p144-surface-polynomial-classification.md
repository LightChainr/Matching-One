# Matching defect as a site-Krushkal rank-image quotient

Status: exact classification for Issue 144 after PR 229 and the unrestricted
Issue 269 theorem.

## Decision

The finite matching polynomial is **not** a direct specialization or source
derivative of the multivariate Tutte, Bollobas-Riordan, Las Vergnas, Krushkal,
topological transition, or 2025 embedded-graph vertex polynomial under the
natural contract relevant here:

1. one fixed embedded object for a lattice quotient;
2. one independent ground element per site;
3. product local weights for the two site states;
4. only variable substitution or source differentiation after the state sum.

The obstruction is not merely that site occupation differs verbally from an
edge subset. One site must realize a typed four-terminal junction `J4`.
An ordinary ribbon edge has two ends, and a four-valent transition state is a
pairing. Neither local state alphabet contains `J4`. A derivative can reweight
states already present but cannot add the missing terminal partition.

There is nevertheless an exact positive identification. The matching
polynomial is a relative-rank derivative of the **rank-image quotient of the
typed incidence-spine site state sum**. Its surface readout is exactly
Krushkal-like: on a torus, the difference between the genera owned by the two
complementary carriers is the matching charge. The state family is new; the
topological exponent is familiar.

## 1. What the named polynomials sum over

The multivariate Tutte/random-cluster polynomial sums over subsets of the
edges of a fixed graph. Bollobas-Riordan and Las Vergnas retain an edge-subset
ground set while enriching the exponent by ribbon or matroid-perspective
data. The Krushkal polynomial has the particularly relevant state sum

```text
K_(G subset Sigma)(x,y,a,b)
 =sum_(A subset E) x^(r(G)-r(A)) y^kappa(A)
                    a^(s(A)/2) b^(s_perp(A)/2),
```

where `s(A)` and `s_perp(A)` are the Euler genera of a regular neighborhood of
the spanning subgraph and of its complement. This is the closest named
readout to the matching problem. See the primary state-sum summary in
[Huggett--Moffatt, arXiv:2212.14233](https://arxiv.org/abs/2212.14233) and the
original orientable-surface construction in
[Krushkal, arXiv:0903.5312](https://arxiv.org/abs/0903.5312).

The 2025 embedded-graph vertex polynomial also remains edge-indexed. Its
two-state form is a sum over partial duals indexed by `A subset E`; its more
general forms partition the ribbon edges among twist/dual states. It is
equivalently described using boundary-component enumeration and the
topological transition polynomial. See
[Yan--Deng--Metsidik, arXiv:2506.07522](https://arxiv.org/abs/2506.07522).
The Las Vergnas extension likewise obtains its recursion from deletion and
contraction of an edge in a matroid perspective; see
[Ellis-Monaghan--Moffatt, arXiv:1311.3762](https://arxiv.org/abs/1311.3762).

These families differ in which ranks, boundary counts, genera or dual ranks
they attach to an edge subset. They share the local problem below.

## 2. Minimal local no-go

The doubled-lattice incidence spine has two site states:

```text
B = J_edge: one block joining N,E,S,W edge-midpoint ports,
W = J_face: one block joining the four face-center ports.
```

For four terminals an absent/present ordinary edge produces either four
singletons or one two-terminal block plus two singletons. The three standard
four-valent transition states are pairings:

```text
(NE)(SW), (NW)(ES), (NS)(EW).
```

None equals the required one-block partition `(NESW)`. Exact enumeration also
shows that at least three independent ordinary edges are needed to connect
four terminals. Such a three-edge gadget has eight subsets. Keeping only its
all-off and all-on macro-states requires deleting six mixed subsets by a block
projector or coefficient extraction. That projector is additional structure;
it is not an ordinary product-local edge specialization.

This proves the no-go for one ground element per site. It also identifies the
minimal repair: promote each site to one **typed partition element** whose
binary values are `J_edge` and `J_face`. Equivalently, a multi-edge gadget may
be used only if it comes with an explicit block-correlation projector.

The claim intentionally excludes vacuous encodings in which a different
unrelated graph is reverse-engineered from each already-known scalar
polynomial. It also does not say that every possible gadget with cancellations
is impossible. It says those gadgets are no longer specializations of the
declared independent edge-subset state sum.

## 3. Exact one-sided rank quotient

For site weights `a,b`, define

```text
R_site(a,b;z)
 =sum_(B subset V) a^|B| b^(N-|B|) z^r_B(B).
```

The unrestricted Issue 269 theorem gives configurationwise, for every finite
integer-period quotient including loops and self-identified faces,

```text
r_B+r_W=2,
q=(r_B-r_W)/2=r_B-1.
```

Therefore

```text
M(a,b)=(z d_z-1) R_site(a,b;z)|_(z=1).              (1)
```

This is the smallest rank-image quotient: complementary rank need not be
stored. Equivalently use the relative Laurent source

```text
Z_rel(a,b;Q)=sum_B a^|B| b^(N-|B|) Q^q,
M(a,b)=Q d_Q Z_rel|_(Q=1).                           (2)
```

For fixed `N`, homogeneity removes one occupation scale, so the independent
terminal variables are only `t=a/b` and one topology source `Q`. In a fully
multivariate local form, replace `a^|B|b^(N-|B|)` by
`product_(i in B) a_i product_(i notin B) b_i`; the topology source remains
one-dimensional.

This is a quotient of the detailed typed partition transfer state. It does
not imply that the connectivity frontier itself has three states.

## 4. Exact Krushkal-type genus derivative

On an honest quotient, let `U_B,V_B` be the complementary cellwise carriers
from digital Alexander duality and let `g_B,g_W` be their genera. For an
arbitrary degenerate quotient, use the exact rank-determined extension
`g_B=1[r_B=2]`, `g_W=1[r_W=2]`. The only possibilities are

```text
(r_B,r_W)  (g_B,g_W)  q
(0,2)       (0,1)     -1
(1,1)       (0,0)      0
(2,0)       (1,0)     +1.
```

Thus

```text
q=g_B-g_W.
```

The site-Krushkal quotient

```text
K_site(a,b;X,Y)
 =sum_B a^|B| b^(N-|B|) X^g_B Y^g_W
```

obeys

```text
M(a,b)=(X d_X-Y d_Y)K_site|_(X=Y=1).                (3)
```

Equation (3) is why Krushkal is the closest named polynomial. Its two genus
variables are exactly the right surface statistic on honest carriers, and the
unrestricted rank theorem extends that readout to every quotient. But
`K_site` is indexed by typed site junction states, not by spanning edge
subsets of a fixed ribbon graph, so it is a Krushkal-type rank-image quotient
rather than a specialization of `K_G`. In particular, this does not claim that
a self-identified base CW presentation itself is a Krushkal ribbon subgraph.

## 5. Exact oracles and minimality

The axis `L=3` oracle exhausts all 512 masks. Equations (1)--(3) reproduce the
same Bernstein coefficients

```text
[-1,-9,-36,-78,-90,-36,36,36,9,1]
```

from black rank, relative rank and carrier-genus derivatives separately.

An independent tiny search establishes that a topology source is necessary.
For every quotient of fixed order the occupation-only site sum is just
`(a+b)^N`. Already at order two, two HNF quotients have the same occupation
counts `[1,2,1]` but different defect Bernstein vectors:

```text
P=((1,0),(0,2)): [-1,0,1],
P=((2,1),(0,1)): [-1,-2,1].
```

Therefore occupation alone cannot recover matching. One extra topology
source is sufficient by (1) or (2), so the terminal extension is minimal.

## Scientific boundary

- Exact: the local `J4` no-go, the rank/genus derivatives, and the order-two
  minimal-source witness.
- Classification: no direct named edge-subset specialization under the stated
  naturality contract.
- Mechanism: the matching defect is a Krushkal-type side-genus derivative of
  a typed **site** state sum.
- Not claimed: a width-independent deletion-contraction algorithm, an
  impossibility theorem for arbitrary cancellation gadgets, or equality with
  the 2025 vertex polynomial after an undeclared block projector.
