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
"$(dirname "$0")/g.sh" "cd $B && rm -rf code/joint_model && tar xzf code.tar.gz -C code && find code -name '._*' -delete; ls code/joint_model | tr '\n' ' '; ls -la data | tail -n +2 | awk '{print \$5, \$9}'"
rm -f "$TAR"
