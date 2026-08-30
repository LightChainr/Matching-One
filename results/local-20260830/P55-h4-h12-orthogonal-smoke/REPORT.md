# P55 H4/H12 target-blind variance smoke

The 20,000-replica-per-design smoke completed locally from source commit
`c1a353a0718d86894ebf49f7b7200152e402ad09`. It used distinct N305/N325 RNG
domains and common random fields only inside each signed orientation pair.
The threshold-rank engine exact self-test passed before acquisition.

The primary output is `variance_power.json`. It deliberately contains no
observed target mean or target score. Centered variance extrapolation freezes
600M paired replicas per design because 300M gives conditional H4-only versus
equal-amplitude H12 Mahalanobis distance 2.768, while 600M gives 3.915 and is
the first preregistered grid point at or above 3.

No Huawei machine was contacted or occupied.
