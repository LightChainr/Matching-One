# Typed serial reversal certificate

For a four-terminal partition in the fixed order `(L0,L1,R0,R1)`, left-right
port reversal is the terminal permutation `(2,3,0,1)`.  The certificate
exhaustively checks all 15 partition states and all 225 ordered serial
products.  It proves

`rev(rev(x)) = x` and `rev(a ∘ b) = rev(b) ∘ rev(a)`.

The committed JSON also freezes the induced state-index permutation and its
fixed-point/two-cycle decomposition.  This is a typed serial anti-involution;
it is not a planar dual, complement dual, periodic gluing construction, or
reliability claim.
