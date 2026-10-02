#!/bin/bash
# collect_benchmarks.sh PROJECT -- write the Slurm accounting of PROJECT's jobs to
# $RUNS/PROJECT/slurm_accounting.csv ('|'-separated, sacct -P), which 09 merges into its benchmarks.
#
# The job ids come from $RUNS/PROJECT/jobs.txt (written by submit_all.sh). MaxRSS and TotalCPU are
# measured by Slurm for the whole job (cgroup), so -- unlike the Julia process's own Sys.maxrss -- they
# include C~scape's and CoralPop's worker processes and ReefMod's MATLAB subprocesses.
set -euo pipefail
source "${CMR_BASE:-$HOME/CoralModelReview}/jobs/common.sh"

PROJECT="${1:-$CMR_PROJECT_NAME}"
JOBS_FILE="$RUNS/$PROJECT/jobs.txt"
[ -f "$JOBS_FILE" ] || { echo "No $JOBS_FILE (was the run submitted with submit_all.sh?)" >&2; exit 1; }

# jobs.txt lines: <name> <jobid>; skip the analysis job itself (still running when this is called)
ids=$(awk '$1 != "analyse" {print $2}' "$JOBS_FILE" | paste -sd, -)
[ -n "$ids" ] || { echo "No job ids in $JOBS_FILE" >&2; exit 1; }

out="$RUNS/$PROJECT/slurm_accounting.csv"
sacct -j "$ids" -P --format=JobID,JobName,State,ExitCode,Elapsed,TotalCPU,CPUTime,MaxRSS,AllocCPUS,ReqMem,NodeList,Start,End > "$out"
echo "Wrote $out ($(($(wc -l < "$out") - 1)) job steps)"
