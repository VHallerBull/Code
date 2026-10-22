# sa_common.sh -- sourced by every Sobol'-run job script in this folder (not run on its own).
#
# Layout on the HPC (RUNNING_SA_ON_HPC.md):
#   ~/CoralModelReview/sa/
#       images/<model>_sa.sif         one image per model
#       jobs/                         these scripts
#       domains/<domain>/             the shared synthetic domains (with sa_inputs/)
#       <ModelFolder>/model_run/<domain>/   designs in, chunks/results/figures out
#       logs/                         Slurm .out files
#
# Slurm runs a COPY of each job script from its spool directory, so they source this file from
# $SA_BASE/jobs (sa_submit.sh exports SA_BASE).

SA_BASE="${SA_BASE:-$HOME/CoralModelReview/sa}"
MODEL="${MODEL:?set MODEL (e.g. coralpop)}"
SIF="${SIF:-$SA_BASE/images/${MODEL}_sa.sif}"
SINGULARITY_MODULE="${SINGULARITY_MODULE:-singularity/singularity.module}"
SA_CODE=/opt/sa/Code/sa                      # the scripts inside the image
export SA_BASE MODEL SIF
export CMR_SA_DOMAINS_DIR="$SA_BASE/domains"  # read by sa_common.jl inside the container
export CMR_SA_RUNS_DIR="$SA_BASE"

# model name -> its folder (same as MODEL_FOLDERS in sa_common.jl)
model_folder() {
    case "$1" in
        coralpop) echo CoralPOP ;; coralblox) echo CoralBlox ;; cscape) echo CScape ;; kora) echo KORA ;;
        corset) echo CORSET ;; reefmod) echo ReefMod ;; coconet) echo CoCoNet ;;
        *) echo "unknown model $1" >&2; return 1 ;;
    esac
}
run_dir() { echo "$SA_BASE/$(model_folder "$MODEL")/model_run/$1"; }

load_singularity() {
    command -v singularity >/dev/null 2>&1 || module load "$SINGULARITY_MODULE"
    [ -f "$SIF" ] || { echo "ERROR: image not found: $SIF" >&2; exit 1; }
}

# container <cmd...> -- run a command inside the model's image with $SA_BASE bound (the CMR_SA_* variables
# pass through).
container() { singularity exec --bind "$SA_BASE" "$SIF" "$@"; }

# timed <logfile> <cmd...> -- /usr/bin/time -v (wall time, peak RSS) around the singularity command.
timed() {
    local log="$1"; shift
    if [ "$1" = "container" ]; then shift; set -- singularity exec --bind "$SA_BASE" "$SIF" "$@"; fi
    if [ -x /usr/bin/time ]; then /usr/bin/time -v -o "$log" "$@"; else "$@"; fi
}
