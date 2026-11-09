# Compare conditional conclusions and quantiles across independent random seeds.
suppressPackageStartupMessages(library(data.table))
second <- 'local_archive/primary_size_uncertainty_2026-10-30/second_seed'
inputs <- fread(file.path(second,'table_s_size_uncertainty_inputs.csv'))
stopifnot(all(unname(tools::md5sum(inputs$path))==inputs$md5))
a <- fread('tables/table_s_size_uncertainty_pooled_summary.csv')
b <- fread(file.path(second,'table_s_size_uncertainty_pooled_summary.csv'))
z <- merge(a,b,by=c('scheme','Taxa','Resp'))
stopifnot(nrow(z)==18,all(abs(z$fraction_q_lt05.x-z$fraction_q_lt05.y) <=
            4 * 0.5 * sqrt(1/z$n_draws.x + 1/z$n_draws.y)),
          all(z$fraction_positive.x==z$fraction_positive.y))
fwrite(z,'tables/table_s_size_uncertainty_seed_comparison.csv')
cat('PASS: second seed preserves pooled directions; detection fractions agree within four worst-case binomial Monte Carlo SEs; differences exported.\n')
