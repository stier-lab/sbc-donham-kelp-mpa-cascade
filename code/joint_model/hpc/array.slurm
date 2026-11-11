#!/bin/bash
# Generic task-list array for the joint model. One line of $JOBLIST per fit:
#   <result file name><TAB><run.py arguments>
# Submit from /home/hpc-scratch/astier/mpa-joint, for example:
#   sbatch --exclude=hpc-13.grit.ucsb.edu --job-name=g9 --cpus-per-task=4 --mem=12G --array=0-55%40 \
#       --export=ALL,JOBLIST=$PWD/jobs/g9.tsv,OUT=$PWD/out/g9 code/joint_model/hpc/array.slurm
# PACK=n runs n consecutive lines per array task side by side (array task i takes lines i*n .. i*n+n-1); ask for
# n * (cores per fit) with --cpus-per-task. GRIT caps RUNNING JOBS per user (150, and the ssh gateway counts as one),
# so pack single-core fits instead of running one job each, and keep the total of all array throttles under about 130.
# A fit whose result exists and validates is skipped; the task exits non-zero if any of its results is missing or invalid.
# hpc-13 has an older CPU on which NumPy 2.x wheels refuse to load: always exclude it.
#SBATCH --partition=grit_nodes
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=12:00:00
#SBATCH --output=/home/hpc-scratch/astier/mpa-joint/logs/%x-%A_%a.log
set -uo pipefail
B=/home/hpc-scratch/astier/mpa-joint
source "$B/venv/bin/activate"
export XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PACK=${PACK:-1}
mkdir -p "$OUT"
# Fallback for result JSON when scratch refuses writes (scratch was 100% full on 2026-11-10). Adrian approved up to
# 30 MB of JSON under his home directory; run.py enforces the cap. Results found there count as done.
export MPA_JOINT_FALLBACK=/home/astier/mpa-joint-out
valid() { [ -s "$1" ] && python3 -c "import json,sys; r=json.load(open(sys.argv[1])); sys.exit(0 if ('E10' in r or r.get('ok')) else 1)" "$1" 2>/dev/null; }
one() {
  local LINE NAME ARGS
  LINE=$(sed -n "$(($1 + 1))p" "$JOBLIST")
  NAME=$(printf '%s' "$LINE" | cut -f1); ARGS=$(printf '%s' "$LINE" | cut -f2-)
  [ -n "$NAME" ] || return 0
  if valid "$OUT/$NAME" || valid "$MPA_JOINT_FALLBACK/$NAME"; then echo "SKIP $NAME (valid result exists)"; return 0; fi
  echo "RUN $(hostname) $(date -Is) :: $NAME :: $ARGS"
  ( cd "$B/code" && eval "python3 joint_model/run.py $ARGS --name ${NAME%.json} --uses $B/data/supported_measurement_uses.parquet --out $OUT" )
  if valid "$OUT/$NAME" || valid "$MPA_JOINT_FALLBACK/$NAME"; then echo "DONE $NAME $(date -Is)"; return 0; fi
  echo "FAILED $NAME: result missing or invalid"; return 1
}
FAIL=0; PIDS=()
for k in $(seq 0 $((PACK - 1))); do one $((SLURM_ARRAY_TASK_ID * PACK + k)) & PIDS+=($!); done
for p in "${PIDS[@]}"; do wait "$p" || FAIL=1; done
exit $FAIL
