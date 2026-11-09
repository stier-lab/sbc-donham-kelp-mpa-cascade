# Review disposition

The completed 999-replicate main output passes all 19 independent tests without skips. The checked model source is recorded in ../independent_checks.json. The earlier implementation review correctly distinguished code verification from inferential validation; its three pending-output skips have now been resolved. The same count covariance is used across calendar-correlation scenarios.

All 17,982 main bootstrap fits converged, with two checked optimizer retries. Four simulated Abalone purple-urchin fits at rho=0.6 reached the computational variance bound and remain flagged; none of the observed-data fits did. No replicate was silently discarded or replaced.

The conditional simulation table and its Monte Carlo uncertainty determine the coverage warning in the report. No working-model or nonlinear coverage certificate is claimed. Fixed allometry, sampling variances, count covariance and expansion weights remain limits. The six histories are fitted separately: cross-species pooled effects require additional covariance work.
