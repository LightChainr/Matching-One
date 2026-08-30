# P250 N101 resolution-only refreeze

The first independent N101 40k stream completed, but its weakest real pair row
had `z=3.919702`, below the frozen `z>=5` requirement.  The scorer therefore
did not construct any cylinder, exponential, held-out-d5, or deck-phase
residual.  The locked score hash is
`f670fe0d0879d4a098ffa455182c17bb8c34f84bdd597562de02d10d58af9134`.

The predeclared square-root rule requires about 65,087 total replicas to reach
five standard errors.  The next simple grid point is an independent fresh 80k
stream, whose projected weakest resolution is 5.543.  It uses a new seed and
does not reuse or append the locked pilot counters.

Nothing in the scientific score changes.  In particular, the primary
parameter-free target remains the thermal `2x=5/4` sine law frozen before any
N101 collection; the old-N65 common sine exponent and exponential mass also
remain byte-for-byte unchanged.  The three models may be evaluated once only
if the fresh 80k resolution gate passes.
