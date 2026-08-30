# Boundary-regular two-carrier systems

The torus data suggested two aggregate inequalities. This note separates their
formal consequence from the topology-specific problem of proving them.

## Basic system

Start with a finite Boolean ground set `E`, a monotone rank map
`r:2^E->{0,1,2}`, and a contiguous nonempty rank-one sector `F`. Define the
complementary carrier by

`r*(W)=2-r(E\W)`.

For `S in F_k`, let `b,d,x,u` be its birth, internal-down, exit, and
internal-up degrees. Thus `b+d=k` and `x+u=N-k`.

The basic axioms alone imply the current identities, complement duality and
edge-pivotal nesting. They do **not** imply ULC.

## Two moment axioms

For an adjacent pair `F_k,F_(k+1)`, write `A_j=|F_j|`, `I_k` for the number of
internal edges, and let all unlabelled sums over `T` run through `F_(k+1)`.

Boundary association is

`BA_k: A_(k+1) sum b(T)x(T) >= (sum b(T))(sum x(T))`.

Transport moment domination is

`TM_k: (N-k)A_k sum d(T)x(T)`

`      >= (N-k-1)I_k sum_(S in F_k)x(S)`.

BA is exactly nonnegative aggregate covariance of upper birth and exit
hazards. TM is exactly

`E_(upper edge)[h_x] >= E_(lower uniform)[h_x]`.

## Exact theorem

If BA and TM hold on every adjacent layer, then

`E_(upper uniform)[h_x] >= E_(upper edge)[h_x]`

`                            >= E_(lower uniform)[h_x]`.

So exit hazard is nondecreasing. If BA and TM hold for both complementary
carriers, complement duality makes the primal birth hazard nonincreasing. The
exact identity

`q_(k+1)/q_k=(1-xi_k)/(1-beta_(k+1))`

then proves ULC. This is a genuine theorem for the new abstract class.

## Why both new axioms are necessary for this modular proof

All nested-upset systems were enumerated exactly through `N=4`: 5 admissible
systems at `N=2`, 111 at `N=3`, and 7,076 at `N=4`. No system through `N=3`
violates BA, TM, exit-hazard monotonicity or ULC.

At `N=4`, each boundary becomes sharp:

- Sector `[1,6,9,10,12,14]` satisfies TM `6>=4`, violates BA `32<35`, and its
  exit hazard falls by `1/24`.
- Sector `[1,2,6,10,14]` satisfies BA with equality `4=4`, violates TM
  `12<16`, and its exit hazard falls by `1/6`.
- Sector `[1,5,9,13,14]` directly violates ULC:
  `q_1=1/4`, `q_2=1/3`, `q_3=1/2`, with margin `-1/72`.

Thus monotonic rank, complement duality, order-convexity and continuity are
not enough even together. BA and TM are logically independent aggregate
inputs. They are sufficient rather than individually necessary for every
possible proof, but neither can simply be deleted from this two-step theorem.

## Topological proof targets

Every one of the 984 existing torus carrier-layer pairs satisfies BA and TM.
What remains is no longer an amorphous covariance problem:

1. prove `A sum bx >= (sum b)(sum x)` by an aggregate pairing of discordant
   and concordant boundary pairs;
2. prove the integer TM inequality by double-counting length-two internal/exit
   paths.

The symbolic counterexamples show that either proof must use genuine homology
carrier geometry, not only the Boolean order axioms.
