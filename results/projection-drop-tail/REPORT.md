# Removing the topological projection reverses the strong-coupling U tail

The fixed N25 histogram gives **U_drop/A = 625/384 lambda^(42/5)
+ higher terms**. Its leading coefficient is positive. The already proved
projected law has **U_star/A = -625/1152 lambda^11 + O(lambda^13)**.
These two fixed, bulk-pressure-equivalent laws have opposite eventual signs
for the original root/slope-normalized angular observer.

No coupling point, series fit or configuration was added. The minimum of
eta=g-r+2K/N was extracted exactly from each saved integer histogram. The
axis leading support has K=5,g=9,count=10; after partition normalization
its E_d coefficient is75/2. Dividing by Q_d=25/2 and Delta=1152/625 gives
1875/1152=625/384. Root motion begins later and cannot cancel this leading term.

The comparison law is the previously defined Sdrop=Sstar+r, obtained solely
by omitting m^(-r). The thermal variable d=p/((1-p)m)*m^(2/N) is common to
both geometries, so its Jacobian cancels in U. The result is not a new
rank-fugacity fit, an arbitrary observer or a fixed-root substitute.

## Cross-size consequences are theory, not new measurements

For an axis L×L torus with L>=5 and same-area companion ell1>=L+2,

`U_drop/A_N = (L-2)/Delta * lambda^(2L-2+2/L) + O(lambda^(2L-2+4/L))`.

Compared with the original negative law,
`U_drop/U_star ~ -(L-2)/(L²-6L+6) * exp((3-2/L)t)` at each fixed L.
N100 and N225 coefficients/exponents are recorded in latest.json as unmeasured
combinatorial predictions. See the [proof](../../notes/topological-projection-reverses-global-u-tail.md).

The first stripe widths now have unequal source costs. Their reciprocal
occupation symmetry no longer cancels the lowest normalized thermal slope.
This identifies how a topological projection can change U even though its
pressure-density discrepancy is bounded by2t/N. That pressure statement takes
fixed t with growing N; it is not an interchange of the two limits.

The finite-t sign-change location, quantitative useful-coupling window and
thermodynamic behavior remain open. These are deterministic views of the
same N25 complete populations, not independent evidence or a continuum field
assignment. P154/P334/F4 decisions stay unchanged.

Reproduce: python scripts/analyze_projection_drop_tail.py --output-dir NEW_DIRECTORY.
