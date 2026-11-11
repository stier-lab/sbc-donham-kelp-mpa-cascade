# Joint hierarchical density model (stage 1)

Design: `docs/JOINT_MODEL_DESIGN.md`. Validation report: `outputs/joint_model_20261109/REPORT.md`. Review-only; not used by the paper.

**Status, 2026-11-10.** The stage-1 model (the default in `model.py`) converges and reproduces, but it fails its placebos:
a regional trend shared by protected and reference stations leaks into the protection estimate (`REPORT.md`, gate 6).
Experimental options that fix or probe this are behind flags (`--drift`, `--regions`, `--inc-scale`, `--zi`); results and
a recommendation are in `outputs/joint_model_20261109/REPORT_structure_experiments.md`. Nothing is adopted: the
structural choices are Adrian's.

- `data.py` builds model inputs from the integrated review's survey-unit table (private, not in git). Options: time placebo
  (`placebo_year`), space placebo (`space_placebo_seed`), leaving out MPAs or stations, and a cap on PISCO counts.
- `model.py` NumPyro model: negative-binomial unit counts; station baseline, region-year effect, station AR(1);
  protection = siting difference + response curve by years since designation + MPA rate deviation + optional MPA level shift
  + shared protected-year shock.
- `run.py simulate` simulates data on the real design at known values (independent NumPy simulator) and refits.
- `run.py fit` fits real data. Sampler settings that met the convergence gate:
  - urchins: `--marg-b --step0 0.01 --accept 0.99 --warmup 3000 --samples 6000` (station baselines integrated out);
  - lobster, sheephead: `--noncentered --level-dev --step0 0.01 --accept 0.99 --warmup 3000 --samples 6000`.
  `--ppc 200` adds posterior predictive checks. `--cpus 4` runs the four chains in parallel.
- `report_x.py` writes the structure-experiment report from `validation/x?/`.
- `report.py` writes the validation report from `outputs/joint_model_20261109/validation/<gate>/*.json`
  (unpack `validation/simulation_replicates.tar.gz` there first).
- `hpc/` runs everything on the UCSB GRIT cluster: `setup_env.slurm` (virtual environment on scratch), `stage.sh` (code and
  the one data file), `make_jobs.py` (task lists), `array.slurm` (task-list array, skip-existing with validation),
  `fetch.sh` (results back over sftp), `pull.sh` (results back over rsync), `put.sh` (single files up, safe on a full
  file system), `g.sh` (ssh with a hard timeout). `hpc/GRIT_NOTES.md` lists what went wrong.

## Stage 2: all five taxa in one model (experimental, started 2026-11-11)

- `data_multi.py` builds the single-taxon inputs for lobster, sheephead, purple urchin, red urchin and kelp, plus index
  maps from each target station-year to the same station in a source taxon one year earlier.
- `model_multi.py` keeps each taxon's single-taxon model (recommended structure) and adds six lagged links between
  latent states at the same station: predators -> urchins, urchins -> kelp. Reports direct, link-carried and total
  protection effects at ten years.
- `run_multi.py fit|simulate` fits or simulates-and-refits. Each chain runs as its own process: NumPyro's parallel
  chains in one process leak memory with five taxa (20 GB in 25 minutes).
- `test_multi.py` checks that, with links off, each taxon's log density equals `model.py`'s, and checks the link maps.
- `report_m.py` writes `outputs/joint_model_20261109/REPORT_multi_taxon.md` from `validation/{m_sim,m_fit}`.
- On GRIT: `--export=ALL,SCRIPT=run_multi.py,...` with `array.slurm`; 4 cores and 24G per fit; task lists `jobs/m_*.tsv`.

Requires Python with `numpyro`, `jax`, `pandas`, `pyarrow`, `scipy`. A converged fit takes 25 to 40 minutes on four cores.
Set `XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1"` to keep each chain on one core.
