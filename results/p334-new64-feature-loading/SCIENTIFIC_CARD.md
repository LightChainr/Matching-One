# New64 separates source energy from contact-geometry loading on original00

- **Changes:** the fixed four-feature contact model captures78.2–80.4% of
  own-source signed clock loading across four size/receiver rows; physical
  source energy alone captures46.3–55.4%. The remaining roughly20% is retained.
- **After baseline clocks:** own safe degree and joint-safe mass keep positive
  partial covariance with own-source center response in every row. Raw safe
  loop and lifetime-response patterns are less uniform.
- **Boundary:** no variance-fraction, complete-closure, causal-allocation,
  independent-population, or field-count claim. Same original prefixes.
- **Source/geometry:** original00 at N325/N425; separate physical sources and
  receivers. Old8 clocks and New64 tangents use conditional cross-stream
  products; clock Gram uses old8 distinct-quartet estimates.
- **Dependency:** `375cd3a1`, `1cfa4ae8`, `8ad30617`; fixed model `011f50e3`;
  original20batch factor and full20000-prefix denominator.
- **Artifact:** `score.json`, raw sufficient rows, partial contact moments,
  LOO/factor and paired differences. No new sampling, raw-fork replay,
  determinant or finite-policy computation; no automatic integration into
  root's frozen analysis.
