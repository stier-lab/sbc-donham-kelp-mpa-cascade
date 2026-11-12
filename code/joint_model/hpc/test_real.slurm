#!/bin/bash
# Run code/joint_model/test_multi.py on the real inputs (needs JAX, so it runs on GRIT, not on the Mac).
#   sbatch --exclude=hpc-13.grit.ucsb.edu code/joint_model/hpc/test_real.slurm
#SBATCH --partition=grit_nodes
#SBATCH --job-name=v-test
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/home/hpc-scratch/astier/mpa-joint/logs/%x-%j.log
B=/home/hpc-scratch/astier/mpa-joint
source "$B/venv/bin/activate"
export XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1" OMP_NUM_THREADS=1
cd "$B/code" && python3 joint_model/test_multi.py "$B/data/supported_measurement_uses.parquet" 2>&1 | grep -v Warning > "$B/out/test_multi_real.log"
