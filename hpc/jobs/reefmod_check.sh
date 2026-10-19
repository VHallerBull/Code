#!/bin/bash
# reefmod_check.sh -- the first run of ReefMod against a real MATLAB, interactively, BEFORE the batch run.
#
#   srun -p defq -c 1 --mem=16G -t 01:00:00 --pty bash ~/CoralModelReview/jobs/reefmod_check.sh
#
# 1. Can the host MATLAB start inside the container (libraries, licence)?
# 2. Does a 2-scenario ReefMod run on the Test domain complete through ADRIA's :reefmod engine?
# 3. Inspect the raw output against what ReefModInterface assumes (hpc/reefmod_inspect.jl): the items
#    ReefMod-GBR/LOCAL_CHANGES.md lists under "Still to verify once a real run exists".
# Output goes to $RUNS/ReefModCheck (separate from the real project). If step 1 fails, use the
# fallback route (submit_all.sh --reefmod-fallback); step 2's MATLAB error message says why.
set -euo pipefail
source "${CMR_BASE:-$HOME/CoralModelReview}/jobs/common.sh"
load_singularity
setup_matlab

echo "== 1. host MATLAB on the host"
matlab -nodisplay -nosplash -batch "disp(version)"
echo "== 1. host MATLAB inside the container"
if container "$REEFMOD_MATLAB_EXE" -nodisplay -nosplash -batch "disp(version); disp(license('test','Statistics_Toolbox'))"; then
    echo "   -> MATLAB starts inside the container: the primary route works."
else
    echo "   -> MATLAB does NOT start inside the container. Use: bash jobs/submit_all.sh <project> --reefmod-fallback"
    exit 1
fi

echo "== 2. two ReefMod scenarios through ADRIA's :reefmod engine"
export CMR_PROJECT_NAME=ReefModCheck CMR_N_SAMPLES=2 CMR_WARMUP=false
mkdir -p "$RUNS/ReefModCheck"
timed "$RUNS/ReefModCheck/time_reefmod.log" container julia --startup-file=no /opt/coral/Code/08_RUN_ENGINE.jl reefmod

echo "== 3. inspect the raw output"
container julia --startup-file=no /opt/coral/Code/hpc/reefmod_inspect.jl "$RUNS/ReefModCheck/reefmod"
grep -E "Elapsed|Maximum resident" "$RUNS/ReefModCheck/time_reefmod.log" 2>/dev/null || true
