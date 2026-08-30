# P250 augmented map result

This directory contains the single completed score under freeze `2bca045` and
scorer `be65bd5`.  `influences.npz` preserves centered delete-one influences
separately for the old 80k and fresh 1.2M sources.  For each candidate and
source, `influence.T @ influence` reconstructs that source's covariance
contribution.  `score.json` contains their sum and the full cross-candidate
covariance.  No new simulation was run.
