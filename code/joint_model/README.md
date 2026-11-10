# Joint hierarchical density model (stage 1)

Design: `docs/JOINT_MODEL_DESIGN.md`. Review-only; not used by the paper.

- `data.py` builds model inputs from the integrated review's survey-unit table (private, not in git).
- `model.py` NumPyro model: negative-binomial unit counts; station baseline, region-year effect, station AR(1);
  protection = siting difference + response curve by years since designation + MPA rate deviation + shared protected-year shock.
- `run.py simulate` simulates data on the real design at known values (independent NumPy simulator) and refits.
- `run.py fit` fits real data (`--curve linear` for comparison with the two-stage slope; `--placebo-year` for the KFM placebo;
  `--noncentered` for taxa whose station process SD is near zero).
- `report.py` writes `outputs/joint_model_20261109/REPORT.md`.

Requires Python with `numpyro`, `jax`, `pandas`, `pyarrow`. A fit takes 5-10 minutes on one core for 4 chains x 800-1000 iterations.
Set `XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1"` to run one fit per core.
