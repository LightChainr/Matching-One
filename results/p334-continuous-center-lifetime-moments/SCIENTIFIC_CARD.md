# The joint center–lifetime moment is now in the same common-label response block

- **Changes:** one read of the existing 1.28M saved tails adds continuous
  `C,W,CW,C²,W²`, two cubic endpoint moments, and the rank-one plateau second
  thermal moment. It retains the exact common order-statistic sampling variance.
- **New joint coordinate:** `E(tau1 tau2)=E(C²)-E(W²)/4`; `E(CW)` alone is an
  endpoint-marginal contrast. Raw second-moment movement still includes center
  translation and is not automatically dispersion response.
- **Observer/sector/source/geometry:** unperturbed baseline and the existing
  common-label Euler-invisible `g±` tangents, each orientation separately;
  N325/N425 original paired quotient geometry and fixed CRN labels.
- **Dependency:** forks `e32a8593`, contacts `959a7fa2`, policy `4db356e1`;
  exactly the same twenty original prefix batches per N. No new MC, DP, raw
  curves, or independent evidence block.
- **Product:** immutable raw-moment lock `f4682eb379b5709a2840faf92beef44ff27f6f23`,
  `batch_moments.json`, with exact integer numerators and source SHA256 hashes.
- **Next readout:** p158 performs pooled, orientation-specific connected
  covariance derivatives with delete-one recentering; root extracts the
  lifetime-weighted plateau centroid/width from this same batch source.
