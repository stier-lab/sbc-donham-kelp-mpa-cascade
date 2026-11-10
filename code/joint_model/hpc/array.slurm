#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=12:00:00
#SBATCH --output=/home/hpc-scratch/astier/mpa-joint/logs/%x-%A_%a.log
# Generic task-list array. Submit from $B with:
#   sbatch --job-name=NAME --array=0-(N-1)%LIMIT --export=ALL,JOBLIST=$B/jobs/NAME.tsv,OUT=$B/out/NAME code/joint_model/hpc/array.slurm
# JOBLIST lines: <expected result file name><TAB><arguments to run.py>
# A task whose result file already exists and validates is skipped, so reruns only fill gaps.
set -uo pipefail
B=/home/hpc-scratch/astier/mpa-joint
source $B/venv/bin/activate
export XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$JOBLIST")
[ -n "$LINE" ] || { echo "FATAL: no line $SLURM_ARRAY_TASK_ID in $JOBLIST"; exit 1; }
RESULT="$OUT/$(printf '%s' "$LINE" | cut -f1)"
ARGS="$(printf '%s' "$LINE" | cut -f2-)"
valid() { [ -s "$1" ] && python3 -c 'import json,sys; r=json.load(open(sys.argv[1])); assert "E10" in r or "ok" in r' "$1" 2>/dev/null; }
mkdir -p "$OUT"
if valid "$RESULT"; then echo "SKIP existing valid $RESULT"; exit 0; fi
echo "RUN $(hostname) $(date -Is) :: $ARGS"
cd $B/code
eval "python3 joint_model/run.py $ARGS --uses $B/data/supported_measurement_uses.parquet --out $OUT"
valid "$RESULT" || { echo "FATAL: missing or invalid result $RESULT"; exit 1; }
echo "DONE $(date -Is) $RESULT"
