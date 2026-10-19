#!/bin/bash
# submit_all.sh PROJECT [options] -- submit every engine job plus the analysis job that waits for them.
#
#   cd ~/CoralModelReview && bash jobs/submit_all.sh TestDomain_N8
#
# Options:
#   --engines "kora corset"   only these engines (default: all seven)
#   --reefmod-fallback        run ReefMod with host MATLAB outside the container (3-stage pipeline:
#                             reefmod_A/B/C) instead of inside it -- see RUNNING_ON_HPC.md section 6e
#   --no-analyse              don't submit the analysis job
#
# Every CMR_* variable set in your shell is passed to the jobs (sbatch exports the environment), e.g.
#   CMR_N_SAMPLES=16 CMR_DHW_IDX=19 bash jobs/submit_all.sh TestDomain_N16_dhw19
#
# Writes $RUNS/PROJECT/jobs.txt (name, job id) and prints how to follow progress. Slurm log files
# (<jobname>_<jobid>.out) go to ~/CoralModelReview/logs.
set -euo pipefail
source "${CMR_BASE:-$HOME/CoralModelReview}/jobs/common.sh"

PROJECT="${1:?usage: submit_all.sh PROJECT [--engines \"...\"] [--reefmod-fallback] [--no-analyse]}"
shift
ENGINES="coralblox cscape coralpop kora corset reefmod coconet"
FALLBACK=0
ANALYSE=1
while [ $# -gt 0 ]; do
    case "$1" in
        --engines) ENGINES="$2"; shift 2 ;;
        --reefmod-fallback) FALLBACK=1; shift ;;
        --no-analyse) ANALYSE=0; shift ;;
        *) echo "unknown option $1" >&2; exit 1 ;;
    esac
done

export CMR_PROJECT_NAME="$PROJECT"
N_SAMPLES="${CMR_N_SAMPLES:-8}"
mkdir -p "$RUNS/$PROJECT" "$CMR_BASE/logs"
cd "$CMR_BASE/logs"
JOBS_FILE="$RUNS/$PROJECT/jobs.txt"
# Keep the ids of engines not being (re-)submitted now, so a later analysis still collects their Slurm
# accounting; drop the ones being replaced (reefmod_A/B/C count as reefmod) and the old analysis.
touch "$JOBS_FILE"
awk -v e=" $ENGINES " '{ k = $1; sub(/_[ABC]$/, "", k) } !(index(e, " " k " ") || $1 == "analyse")' \
    "$JOBS_FILE" > "$JOBS_FILE.tmp" && mv "$JOBS_FILE.tmp" "$JOBS_FILE"

# Per-engine resources for the serial test run (1 CPU each). Memory/time are generous first guesses
# for the 10-site Test domain; the run's own benchmarks (K_benchmarks) tell you what to request next.
#   cscape:  one Distributed worker + this process, each holding a scenario
#   kora:    individual-based, but one 250 m² section per site: ~1 s and ~3 GB per scenario
#   reefmod: one MATLAB at a time, ~1-2 min MATLAB start-up per scenario
resources() {
    case "$1" in
        coralblox) echo "--mem=8G  --time=01:00:00" ;;
        cscape)    echo "--mem=24G --time=06:00:00" ;;
        coralpop)  echo "--mem=8G  --time=02:00:00" ;;
        kora)      echo "--mem=8G  --time=02:00:00" ;;
        corset)    echo "--mem=8G  --time=02:00:00" ;;
        reefmod)   echo "--mem=16G --time=12:00:00" ;;
        coconet)   echo "--mem=8G  --time=01:00:00" ;;
    esac
}

ids=()
for e in $ENGINES; do
    if [ "$e" = "reefmod" ] && [ "$FALLBACK" = "1" ]; then
        a=$(sbatch --parsable --job-name=cmr_reefmod_A --export=ALL "$JOBS_DIR/reefmod_A.slurm")
        b=$(sbatch --parsable --job-name=cmr_reefmod_B --dependency=afterok:$a --array=1-"$N_SAMPLES"%"${REEFMOD_ARRAY_CONCURRENCY:-4}" \
            --export=ALL "$JOBS_DIR/reefmod_B.slurm")
        c=$(sbatch --parsable --job-name=cmr_reefmod_C --dependency=afterok:$b --export=ALL,REEFMOD_STAGE_B_JOB=$b "$JOBS_DIR/reefmod_C.slurm")
        printf "reefmod_A %s\nreefmod_B %s\nreefmod_C %s\n" "$a" "$b" "$c" >> "$JOBS_FILE"
        ids+=("$c")
        echo "reefmod (fallback): A=$a -> B=$b (array 1-$N_SAMPLES) -> C=$c"
        continue
    fi
    # shellcheck disable=SC2046
    j=$(sbatch --parsable --job-name="cmr_$e" $(resources "$e") --export=ALL,ENGINE="$e" "$JOBS_DIR/run_engine.slurm")
    echo "$e $j" >> "$JOBS_FILE"
    ids+=("$j")
    echo "$e: job $j"
done

if [ "$ANALYSE" = "1" ] && [ ${#ids[@]} -gt 0 ]; then
    dep=$(IFS=:; echo "${ids[*]}")
    a=$(sbatch --parsable --dependency=afterany:"$dep" --export=ALL "$JOBS_DIR/analyse.slurm")
    echo "analyse $a" >> "$JOBS_FILE"
    echo "analyse: job $a (starts when all of the above have finished, successful or not)"
fi

cat <<EOF

Submitted project $PROJECT ($N_SAMPLES scenarios per engine). Job ids: $JOBS_FILE
  squeue -u \$USER                                   # queued / running
  sacct -j $(awk '{print $2}' "$JOBS_FILE" | paste -sd, -) --format=JobName,State,Elapsed,MaxRSS
  tail -f $CMR_BASE/logs/cmr_<engine>_<jobid>.out    # one engine's log
Results: $RUNS/$PROJECT/   (FIGURES.md + images/ once cmr_analyse has finished)
EOF
