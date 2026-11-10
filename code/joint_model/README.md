# Joint hierarchical density model (stage 1)

Design: `docs/JOINT_MODEL_DESIGN.md`. Validation report: `outputs/joint_model_20261109/REPORT.md`. Review-only; not used by the paper.

**Status, 2026-11-10.** The model converges and reproduces, but it fails its placebos: a regional trend shared by protected
and reference stations leaks into the protection estimate (report, gate 6). Do not present its E10 values as protection
effects until that is fixed and the space placebo centres on zero.

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
- `report.py` writes the validation report from `outputs/joint_model_20261109/validation/<gate>/*.json`
  (unpack `validation/simulation_replicates.tar.gz` there first).
- `hpc/` runs everything on the UCSB GRIT cluster: `setup_env.slurm` (virtual environment on scratch), `stage.sh` (code and
  the one data file), `make_jobs.py` (task lists), `array.slurm` (task-list array, skip-existing with validation),
  `fetch.sh` (results back over sftp), `g.sh` (ssh with a hard timeout). `hpc/GRIT_NOTES.md` lists what went wrong.

Requires Python with `numpyro`, `jax`, `pandas`, `pyarrow`, `scipy`. A converged fit takes 25 to 40 minutes on four cores.
Set `XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1"` to keep each chain on one core.
