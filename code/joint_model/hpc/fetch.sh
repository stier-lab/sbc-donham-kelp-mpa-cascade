#!/bin/bash
# Pull result JSON files (not posterior draws) for one or more gates from GRIT scratch into
# outputs/joint_model_20261109/validation/<gate>/. Uses sftp, which is not wrapped in srun on the GRIT login
# node, so it still works when the per-user running-job limit blocks ssh commands.
# Usage: bash code/joint_model/hpc/fetch.sh g3 g4a [logs]
set -uo pipefail
B=/home/hpc-scratch/astier/mpa-joint
V=outputs/joint_model_20261109/validation
CMDS=$(mktemp)
for g in "$@"; do
  if [ "$g" = logs ]; then mkdir -p "$V/logs"; echo "get $B/logs/*.log $V/logs/" >> "$CMDS"
  else mkdir -p "$V/$g"; echo "get $B/out/$g/*.json $V/$g/" >> "$CMDS"; fi
done
perl -e 'alarm shift; exec @ARGV' "${GRIT_TIMEOUT:-70}" sftp -q -i ~/.ssh/grit -o BatchMode=yes -o ConnectTimeout=12 -b "$CMDS" astier@hpc.grit.ucsb.edu > /dev/null 2>&1
RC=$?; rm -f "$CMDS"
for g in "$@"; do [ "$g" = logs ] || echo "$g: $(ls "$V/$g" 2>/dev/null | wc -l | tr -d ' ') files"; done
exit $RC
