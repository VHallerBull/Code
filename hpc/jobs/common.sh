# common.sh -- sourced by every job script in this folder (not run on its own).
#
# Slurm runs a COPY of each job script from its spool directory, so the scripts cannot find this file
# relative to themselves; they source it from $CMR_BASE/jobs (submit_all.sh exports CMR_BASE; the
# default below matches RUNNING_ON_HPC.md's layout).

CMR_BASE="${CMR_BASE:-$HOME/CoralModelReview}"
SIF="${SIF:-$CMR_BASE/coral_engines.sif}"
RUNS="${RUNS:-$CMR_BASE/runs}"
JOBS_DIR="${JOBS_DIR:-$CMR_BASE/jobs}"
SINGULARITY_MODULE="${SINGULARITY_MODULE:-singularity/singularity.module}"
MATLAB_MODULE="${MATLAB_MODULE:-matlab}"
export CMR_BASE SIF RUNS JOBS_DIR
export CMR_MODEL_RUNS_DIR="$RUNS"             # read by paths.jl inside the container
export CMR_PROJECT_NAME="${CMR_PROJECT_NAME:-TestDomain_N8}"

EXTRA_BINDS=""

load_singularity() {
    command -v singularity >/dev/null 2>&1 || module load "$SINGULARITY_MODULE"
    [ -f "$SIF" ] || { echo "ERROR: image not found: $SIF" >&2; exit 1; }
}

# container <cmd...> -- run a command inside the image with the runs folder (and MATLAB, if set up)
# bound. The host environment (CMR_*, licence variables, REEFMOD_MATLAB_EXE) passes through.
container() {
    local binds="$RUNS"
    [ -n "$EXTRA_BINDS" ] && binds="$binds,$EXTRA_BINDS"
    singularity exec --bind "$binds" "$SIF" "$@"
}

# setup_matlab -- load the host MATLAB and make it usable INSIDE the container: bind its install
# folder (and any licence file outside it) and point ReefModInterface at its executable. The licence
# variables the module sets (MLM_LICENSE_FILE / LM_LICENSE_FILE) pass through unchanged.
setup_matlab() {
    module load "$MATLAB_MODULE"
    local exe
    exe="$(readlink -f "$(command -v matlab)")" || { echo "ERROR: no matlab after 'module load $MATLAB_MODULE'" >&2; exit 1; }
    MATLAB_ROOT="${MATLAB_ROOT:-$(dirname "$(dirname "$exe")")}"
    export MATLAB_ROOT
    export REEFMOD_MATLAB_EXE="$MATLAB_ROOT/bin/matlab"
    EXTRA_BINDS="$MATLAB_ROOT"
    local lic entry
    for lic in "${MLM_LICENSE_FILE:-}" "${LM_LICENSE_FILE:-}"; do
        IFS=':' read -ra entries <<< "$lic"
        for entry in "${entries[@]}"; do
            # port@host entries need no bind; files/folders outside MATLAB_ROOT do
            if [ -e "$entry" ] && [[ "$entry" != "$MATLAB_ROOT"* ]]; then
                EXTRA_BINDS="$EXTRA_BINDS,$(dirname "$entry")"
            fi
        done
    done
    echo "MATLAB: $REEFMOD_MATLAB_EXE (binds: $EXTRA_BINDS; MLM_LICENSE_FILE=${MLM_LICENSE_FILE:-} LM_LICENSE_FILE=${LM_LICENSE_FILE:-})"
}

# timed <logfile> <cmd...> -- run under /usr/bin/time -v when available (wall time, peak RSS).
timed() {
    local log="$1"; shift
    if [ -x /usr/bin/time ]; then /usr/bin/time -v -o "$log" "$@"; else "$@"; fi
}
