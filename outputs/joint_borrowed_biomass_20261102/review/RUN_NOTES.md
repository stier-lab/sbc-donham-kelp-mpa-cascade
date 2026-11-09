# Execution and review notes

The initial 199-replicate calculation was expanded to 999 replicates for the main correlation scenarios and second-seed check, and 499 for pooling/timing sensitivities. A larger-run simulated fit exhausted its first two line-search attempts and stopped the model run. The optimizer was then given a longer line search and a real constrained SLSQP retry. An independent forced-failure test checks that retries are real optimizations and that exhaustion still raises an error. Final outputs must have a completion marker and source hashes before use.

Count-only sampling covariance is held identical across correlation scenarios, so those comparisons isolate the process-correlation assumption. The initial coarse simulation screen used 20 datasets and 39 inner replicates; final calibration uses 20 datasets and 199 inner replicates with independent per-dataset seeds. This remains a small calibration exercise with reported Monte Carlo bounds, not a certificate of nominal coverage.

No primary model, manuscript, supplement, canonical raw source, researcher exclusion or protected-site borrowing policy was changed. The working model changes the representation of already supported size borrowing: it uses a within-role hierarchical size process rather than treating fixed donor means as known.
