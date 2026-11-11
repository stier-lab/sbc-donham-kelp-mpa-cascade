#!/bin/bash
# Stage code (and, with "data", the one restricted input file) to GRIT scratch. Run from the repository root on the Mac.
set -euo pipefail
B=/home/hpc-scratch/astier/mpa-joint
RSH="ssh -i $HOME/.ssh/grit -o BatchMode=yes -o ConnectTimeout=12"
TAR=$(mktemp -t mpa_joint_code).tar.gz
COPYFILE_DISABLE=1 tar --no-xattrs --exclude='._*' --exclude='.DS_Store' --exclude='__pycache__' -czf "$TAR" -C code joint_model
FILES=("$TAR")
[ "${1:-}" = data ] && FILES+=("outputs/integrated_review_20261103/private/supported_measurement_uses.parquet")
rsync -a -e "$RSH" "$TAR" astier@hpc.grit.ucsb.edu:$B/code.tar.gz
[ "${1:-}" = data ] && rsync -a -e "$RSH" outputs/integrated_review_20261103/private/supported_measurement_uses.parquet astier@hpc.grit.ucsb.edu:$B/data/
[ -d code/joint_model/jobs ] && rsync -a -e "$RSH" code/joint_model/jobs/ astier@hpc.grit.ucsb.edu:$B/jobs/
# Extract beside the live copy, check it, then swap: scratch was 100% full on 2026-11-10 and a failed extract must not
# leave the cluster without code.
SUM=$(cat code/joint_model/run.py code/joint_model/model.py code/joint_model/data.py | md5 -q)
"$(dirname "$0")/g.sh" "cd $B && rm -rf code.new && mkdir -p code.new code && tar xzf code.tar.gz -C code.new && find code.new -name '._*' -delete && [ -s code.new/joint_model/hpc/array.slurm ] && [ \"\$(cat code.new/joint_model/run.py code.new/joint_model/model.py code.new/joint_model/data.py | md5sum | cut -c1-32)\" = \"$SUM\" ] && { rm -rf code.old; [ -d code/joint_model ] && mv code/joint_model code.old; mv code.new/joint_model code/joint_model && rmdir code.new && echo STAGED_OK; } || echo STAGE_FAILED_live_copy_untouched; df -h /home/hpc-scratch | tail -1 | awk '{print \"scratch free:\", \$4}'"
rm -f "$TAR"
