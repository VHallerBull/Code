#!/bin/bash
# run_engine.slurm -- run ONE engine (08_RUN_ENGINE.jl) inside the container.
#
# Normally submitted by submit_all.sh, which sets the job name (cmr_<engine>) and resources per engine.
# By hand, from ~/CoralModelReview/logs:
#   sbatch --job-name=cmr_kora --export=ALL,ENGINE=kora,CMR_PROJECT_NAME=TestDomain_N8 ../jobs/run_engine.slurm
#
# ENGINE=reefmod additionally loads the host MATLAB module and binds it into the container
# (setup_matlab in common.sh).
#
# The resource lines below are the test-run defaults (1 CPU = engines run serially, the fair basis for
# the per-scenario benchmarks); submit_all.sh overrides them per engine.
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.out
# --output is flat on purpose: Slurm opens it before this script runs, so it must be a folder that
# exists -- the one you submit from (submit_all.sh submits from ~/CoralModelReview/logs).

set -euo pipefail
source "${CMR_BASE:-$HOME/CoralModelReview}/jobs/common.sh"

ENGINE="${ENGINE:?set ENGINE=coralblox|cscape|coralpop|kora|corset|reefmod}"
echo "[$(date '+%F %T')] $ENGINE on $(hostname), job ${SLURM_JOB_ID:-none}, project $CMR_PROJECT_NAME, ${SLURM_CPUS_PER_TASK:-1} CPU"

load_singularity
[ "$ENGINE" = "reefmod" ] && setup_matlab
mkdir -p "$RUNS/$CMR_PROJECT_NAME"

timed "$RUNS/$CMR_PROJECT_NAME/time_$ENGINE.log" \
    container julia --startup-file=no -t "${SLURM_CPUS_PER_TASK:-1}" /opt/coral/Code/08_RUN_ENGINE.jl "$ENGINE"

echo "[$(date '+%F %T')] $ENGINE finished"
