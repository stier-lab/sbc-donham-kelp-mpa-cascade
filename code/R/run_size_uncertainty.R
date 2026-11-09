# Optional raw-source sensitivity: run from the analysis repository root.
# It writes supplemental diagnostics, never harmonized inputs or primary fits.
source('code/R/55_primary_size_inventory.R')
status <- system2(file.path(R.home('bin'),'Rscript'), 'code/R/prepare_size_uncertainty_domain.R')
if(status!=0L)stop('Full annual-domain reconstruction failed; no size refits run.')
source('code/R/57_primary_size_uncertainty.R')
source('code/tests/check_primary_size_uncertainty.R')
status <- system2(file.path(R.home('bin'),'Rscript'),'code/R/57_primary_size_uncertainty.R',
  env=c('DONHAM_SIZE_DRAWS=250','DONHAM_SIZE_SEED=20261031',
    'DONHAM_SIZE_OUTPUT_DIR=local_archive/primary_size_uncertainty_2026-10-30/second_seed'))
if(status!=0L)stop('Second-seed size sensitivity failed.')
source('code/tests/check_primary_size_uncertainty_seeds.R')
