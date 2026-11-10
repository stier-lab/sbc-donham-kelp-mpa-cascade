#!/bin/bash
#SBATCH --job-name=mpa-setup
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/home/hpc-scratch/astier/mpa-joint/logs/setup-%j.log
# Build the Python environment on scratch (home is quota-limited). Run once.
set -uo pipefail
B=/home/hpc-scratch/astier/mpa-joint
export PIP_CACHE_DIR=$B/pip-cache XDG_CACHE_HOME=$B/cache TMPDIR=$B/tmp
mkdir -p "$PIP_CACHE_DIR" "$XDG_CACHE_HOME" "$TMPDIR"
hostname; python3 --version
[ -x $B/venv/bin/python ] || python3 -m venv $B/venv || { echo "FATAL: venv creation failed"; exit 1; }
$B/venv/bin/python -m pip install --upgrade pip >/dev/null
$B/venv/bin/python -m pip install numpyro jax pandas pyarrow scipy || { echo "FATAL: pip install failed"; exit 1; }
$B/venv/bin/python - <<'PY' || { echo "FATAL: import check failed"; exit 1; }
import numpyro, jax, pandas, pyarrow, scipy, numpy
print("VERSIONS", "numpyro", numpyro.__version__, "jax", jax.__version__, "pandas", pandas.__version__,
      "pyarrow", pyarrow.__version__, "scipy", scipy.__version__, "numpy", numpy.__version__)
import jax.numpy as jnp
print("JAX_OK", float(jnp.arange(4.0).sum()))
PY
$B/venv/bin/python -m pip freeze > $B/requirements.lock
echo SETUP_DONE
