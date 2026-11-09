# Independent diagnostics review

The ecological-stats diagnostics reviewer checked the implementation twice. The first review identified annual multiplication before visit matching, zero-count certainty, out-of-cohort calibration contributions, and a shared-transect weight that canceled because the source contains one transect per site. The final implementation replaces these with visit products, a conditional count-rate model, retained-cohort calibration, and explicit survey-replication limits. Five comparison-years lack full visit-level size support under the seven-day rule.

The second review verified the 52 comparison-years, source hashes, annual products, and reported interval/boundary conclusions. Follow-up corrections distinguish unique fish/size rows/size visits from reused visit applications, use a separate random stream for count–size dependence scenarios, and describe the count prior as a Gamma shape increment (rather than a literal half fish).

Neither review validates the roving survey as representative of the count transect. Working-model intervals remain conditional; variance boundaries and interval-method sensitivity preclude treating these as settled primary protection effects.
