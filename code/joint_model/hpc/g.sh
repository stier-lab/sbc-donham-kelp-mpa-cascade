#!/bin/bash
# Run one short command on the GRIT login host with a hard wall-clock cap (see the grit-cluster skill).
# Usage: hpc/g.sh 'squeue -u astier'      (run from this Mac; needs campus network or the UCSB VPN)
exec perl -e 'alarm shift; exec @ARGV' "${GRIT_TIMEOUT:-45}" ssh -i ~/.ssh/grit -o BatchMode=yes -o ConnectTimeout=12 astier@hpc.grit.ucsb.edu "$1" 2> >(grep -v '^srun:' >&2)
