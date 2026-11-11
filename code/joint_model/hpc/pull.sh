#!/bin/bash
# Pull result JSON for the named gates with rsync (only new files travel; use when ssh works). Also collects any
# results that went to the home-directory fallback. Usage: bash code/joint_model/hpc/pull.sh xA xC xE
V=outputs/joint_model_20261109/validation
SSH="ssh -i $HOME/.ssh/grit -o BatchMode=yes -o ConnectTimeout=12"
for g in "$@"; do
  mkdir -p "$V/$g"
  perl -e 'alarm 150; exec @ARGV' rsync -a --include='*.json' --exclude='*' -e "$SSH" astier@hpc.grit.ucsb.edu:/home/hpc-scratch/astier/mpa-joint/out/$g/ "$V/$g/" 2>&1 | grep -v -E "dyld|^srun"
  perl -e 'alarm 60; exec @ARGV' rsync -a --include="${g}_*.json" --exclude='*' -e "$SSH" astier@hpc.grit.ucsb.edu:/home/astier/mpa-joint-out/ "$V/$g/" 2>/dev/null
  find "$V/$g" -name '*.json' -size 0 -delete
  echo "$g: $(ls "$V/$g" | wc -l | tr -d ' ') files"
done
