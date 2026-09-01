# Why global component labels cannot be the P537 four-arm landing state

## Exact obstruction

Let `z` be vacant on a square-lattice torus and suppose its four nearest
neighbours alternate occupied/vacant.  Define

- `b_z=1` when the two occupied neighbours belong to distinct occupied
  components before `z` is added, and `b_z=0` otherwise;
- `w_z=1` when the two vacant neighbours belong to distinct components of
  the matching graph after `z` is removed from that graph, and `w_z=0`
  otherwise.

For the digital matching rank `q=r-1`, direct component/Euler bookkeeping
gives the one-site identity

\[
 q(A\cup\{z\})-q(A)=1-b_z-w_z.                              \tag{1}
\]

Indeed, adding `z` joins the two occupied arms exactly when `b_z=1`, so the
occupied component count changes by `-b_z`.  Removing the formerly vacant
`z` splits its matching component exactly when `w_z=1`, so the matching
component count changes by `+w_z`.  The alternating mask adds one vertex and
two occupied edges but no occupied square, hence its occupied Euler
characteristic changes by `-1`.  Substitution in

\[
 q=C_{occ}-C_{match}-\chi_{occ}
\]

proves (1).

The homology rank `r` is monotone under addition of occupied sites.  Therefore
the left side of (1) is nonnegative and

\[
                         \boxed{b_z+w_z\le1}.                \tag{2}
\]

In particular, the state with two globally distinct occupied branches *and*
two globally distinct vacant matching separators is empty on every finite
torus.  This is not a small-N accident.

## N25 manifestation

The exact fixed-`x,z` traversals on both `(5,0)` and `(4,3)` N25 tori contain
`2^20=1,048,576` alternating backgrounds apiece.  Their global landing types
partition as follows:

| geometry | occupied distinct only | neither distinct | vacant distinct only | both distinct |
|---|---:|---:|---:|---:|
| `(5,0)` | 522,883 | 419,177 | 106,516 | 0 |
| `(4,3)` | 534,472 | 404,816 | 109,288 | 0 |

The counts are a finite illustration of (2), not the proof.

## Consequence for the landing matrix

The ordinary four-arm label in Issue #537 must live at a **finite landing
cut** (a collar or annulus around `z`) before the arms reconnect in the outer
torus.  Replacing those cut labels by global occupied/matching component IDs
makes the intended row identically empty.  Conversely, retaining only local
occupation alternation while ignoring cut identities defines a nonempty
`near_block`, but it is a relaxation and cannot be called the canonical
ordinary-four-arm block.

This changes the computation target sharply: the next exact producer should
propagate four labelled arms only to the first collar boundary, record their
partition there, and then attach the outer connectivity as a separate state.
The outer attachment is precisely where the rank transition occurs; merging
it into the landing label erases the state that the minor is meant to test.
