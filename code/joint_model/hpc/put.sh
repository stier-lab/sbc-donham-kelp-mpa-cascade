#!/bin/bash
# Copy single files to the staged code on GRIT over sftp, safely on a full file system: upload under a new temporary
# name each try (BeeGFS pins a file to a storage target at creation), check the size, then rename over the live file.
# Usage (repository root): bash code/joint_model/hpc/put.sh code/joint_model/run.py code/joint_model/hpc/array.slurm
B=/home/hpc-scratch/astier/mpa-joint
RC=0
for f in "$@"; do
  L=$(wc -c < "$f" | tr -d ' '); R="$B/code/${f#code/}"; OK=0
  for i in 1 2 3 4 5 6 7 8; do
    T="$R.$$.$i.part"; C=$(mktemp)
    printf 'put %s %s\nls -l %s\n' "$f" "$T" "$T" > "$C"
    S=$(perl -e 'alarm 25; exec @ARGV' sftp -q -i ~/.ssh/grit -o BatchMode=yes -b "$C" astier@hpc.grit.ucsb.edu 2>/dev/null | awk -v t="$T" 'NF>=9 && $NF==t {print $5}' | tail -1)
    if [ "$S" = "$L" ]; then
      printf 'rename %s %s\n' "$T" "$R" > "$C"   # OpenSSH sftp rename replaces the target (posix-rename extension)
      perl -e 'alarm 25; exec @ARGV' sftp -q -i ~/.ssh/grit -o BatchMode=yes -b "$C" astier@hpc.grit.ucsb.edu >/dev/null 2>&1 && OK=1
    else
      printf 'rm %s\n' "$T" > "$C"; perl -e 'alarm 25; exec @ARGV' sftp -q -i ~/.ssh/grit -o BatchMode=yes -b "$C" astier@hpc.grit.ucsb.edu >/dev/null 2>&1
    fi
    rm -f "$C"; [ $OK = 1 ] && break
  done
  [ $OK = 1 ] && echo "put ok: $f ($L bytes, try $i)" || { echo "PUT FAILED: $f"; RC=1; }
done
exit $RC
