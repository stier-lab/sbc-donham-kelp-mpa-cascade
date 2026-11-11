# Running Python (NumPyro + JAX) on GRIT: what I learned, 2026-11-10

For Adrian to fold into the `grit-cluster` skill. Verified on this project's runs unless marked otherwise.

## Limits that bite
- **Per-user QOS `grit_normal`: 300 cores, 1,200 GB memory, 150 running jobs, 5,000 submitted.**
  Check with `sacctmgr -n -P show qos format=Name,MaxTRESPU,MaxJobsPU,MaxSubmitPU`.
- **The ssh gateway counts against both limits.** Every `ssh host 'cmd'` is an `srun` job named `ssh-cmd` (1 core, 1 job).
  If your own arrays fill 150 jobs or 300 cores, your next ssh command pends behind them (`QOSMaxJobsPerUserLimit` or
  `QOSMaxCpuPerUserLimit`) and you are locked out of `sbatch`, `scancel` and `squeue` until the arrays drain.
  Keep the sum of all array throttles (`%N`) at or below about 130 jobs and 290 cores.
- **sftp is not wrapped in srun.** `sftp -b cmds.txt astier@hpc.grit.ucsb.edu` works during a lock-out: it can list, read and
  write files on scratch (logs, results, task lists). It cannot run Slurm commands. It is also the dependable way to pull
  results (`get dir/*.json`), but it re-downloads everything each time; use `rsync` for thousands of files when ssh works.
- **Getting out of a lock-out without waiting.** My arrays read their task list at start-up, so I edited the task list over
  sftp to make a few upcoming tasks run `scancel --state=PENDING <my array id>`. That cancelled only my own pending tasks
  and freed the gateway. Simpler prevention: never submit an unthrottled array.
- Pack single-core fits: one array task that runs 8 fits side by side on `--cpus-per-task=8` uses one job slot instead of
  eight (`PACK=8` in `hpc/array.slurm`). The core limit, not the job limit, then sets throughput.

## Node problems
- **`hpc-13` cannot load NumPy 2.x wheels** (`RuntimeError: NumPy was built with baseline optimizations (X86_V2) but your
  machine doesn't support (X86_V2)`). Tasks sent there fail in two seconds, so the node drains an array fast: it ate about
  370 of my tasks. Use `--exclude=hpc-13.grit.ucsb.edu`. Node names must be fully qualified; `--exclude=hpc-13` is
  rejected as an invalid node name. The other grit_nodes hosts I landed on (06, 07, 08, 09, 14 to 18) were fine.
- Nodes disagree on time zone in `date` output (some UTC, some Pacific). Do not compare log timestamps across nodes.

## Python environment
- `/usr/bin/python3` is 3.12.3 on the nodes and `python3 -m venv` works. Building a venv with numpyro, jax, pandas, pyarrow
  and scipy on `/home/hpc-scratch/astier` took 4 minutes in an sbatch job; compute nodes reach PyPI.
- Put the pip cache and `TMPDIR` on scratch too (`PIP_CACHE_DIR`, `TMPDIR`); `/home/astier` is full.
- Versions that worked: numpyro 0.22.0, jax 0.12.0.1, numpy 2.5.4, pandas 3.1.1, pyarrow 26.0.1, scipy 1.18.1.
- One chain per core: `XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1"` plus
  `OMP_NUM_THREADS=1`. For parallel chains call `numpyro.set_host_device_count(n)` before importing jax and ask Slurm for
  `--cpus-per-task=n`; four parallel chains ran roughly 3.5 times faster than four sequential chains on comparable fits (not a controlled timing).
- The purple urchin fit reproduced to the fourth decimal between my cloud sandbox and GRIT with the same package versions.

## Workflow notes
- `sbatch` copies the batch script at submission, so restaging code while an array runs is safe for running tasks, but
  tasks that start later pick up the new code. Record which code version each batch used.
- A task-list array (`<result name><TAB><arguments>` per line) with skip-if-valid made every rerun cheap: resubmit the
  whole array and only missing results run.
- Make every task exit non-zero when its result file is missing or fails to parse; `sacct` state alone hid nothing then.
- Slurm `MaxArraySize` is 5001.

## When scratch fills up (it did, 2026-11-10 about 13:50 to at least 17:00)
- `/home/hpc-scratch` is one shared 53 TB BeeGFS volume. When it reaches 100%, writes fail at random with
  `No space left on device` even though `df` still shows a few GB free. Symptoms: empty Slurm logs, zero-byte or
  truncated result files, jobs failing in one second with exit code `0:53`, `sbatch: Batch script is empty`.
- **A file is pinned to a storage target when it is created.** Rewriting the same file retries the same full target and
  fails forever. Retry by writing a NEW file name each attempt, check its size, then rename it into place
  (`write_retry` in `run.py`; `hpc/put.sh` does the same for uploads over sftp).
- Save the small result (JSON) before any large optional output. My first experiment batch lost 27 finished fits because
  the posterior-draw `.npz` was written first and its failure raised before the JSON was saved.
- Never extract a tarball over the live code directory: a partial extract left an empty `array.slurm` and one bad
  staging run could have removed the code under running arrays. `stage.sh` now extracts beside the live copy, checks
  checksums and sizes, then swaps.
- Fallback that Adrian approved for this project: result JSON only, capped at 30 MB, to `/home/astier/mpa-joint-out`
  after five failed scratch attempts (`MPA_JOINT_FALLBACK` in `array.slurm`). In the end it was not needed once the
  retry used fresh file names.
- sftp to the login node can itself hang for a minute or more while the file systems are in this state.
- Other sessions under the same account share the 300-core and 150-job limits (jobs named `dafl-dd` appeared during
  these runs). Leave headroom for them and for the ssh gateway.

## Not GRIT, but it cost time
- Copying files from the Claude cloud sandbox to the Mac right after writing them delivered stale versions twice. Check the
  checksum on the Mac after every copy.
- The Mac bridge drops any command that runs longer than about 80 seconds: run long ssh or rsync calls in the background
  and read a log file.
