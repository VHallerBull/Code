#!/usr/bin/env bash
# Local pick-and-freeze pilot: designs, then N_PROCS parallel Julia processes for the walks and then for
# the zero-DHW reference (CoralPop runs through
# Python, so one process per share, no threads), then collect. Run from the CoralModelReview folder:
#
#   bash Code/sa_pf/run_pilot.sh [n_base=256] [n_procs=12]
#
# Bands 0.1, 0.2, 0.5 on both domains. A rerun resumes: finished chunk files are skipped.
set -euo pipefail
N_BASE=${1:-256}
N_PROCS=${2:-12}
ENVP=Code/sa/shapley
BANDS="0.1 0.2 0.5"
LOGS=CoralPOP/model_run_pf_pilot/logs
mkdir -p "$LOGS"

julia +release --project=$ENVP Code/sa_pf/02_design.jl coralpop "$N_BASE" $BANDS

JOBS=""
for dom in SingleReef_SobolDHW NineReefs_SobolDHW; do
  for b in $BANDS; do JOBS="$JOBS $b:$dom"; done
done

pids=()
for p in $(seq 1 "$N_PROCS"); do
  julia +release --project=$ENVP Code/sa_pf/03_run_chunk.jl coralpop "$p" "$N_PROCS" $JOBS \
    > "$LOGS/proc_$(printf %02d "$p").log" 2>&1 &
  pids+=($!)
done
fail=0
for pid in "${pids[@]}"; do wait "$pid" || fail=1; done
[ "$fail" -eq 0 ] || { echo "a run process failed: see $LOGS"; exit 1; }

# zero-DHW reference: each walk's A and B parameter sets without heat stress (2 runs per base sample)
pids=()
for p in $(seq 1 "$N_PROCS"); do
  julia +release --project=$ENVP Code/sa_pf/03b_run_zero.jl coralpop "$p" "$N_PROCS" $JOBS     > "$LOGS/zero_$(printf %02d "$p").log" 2>&1 &
  pids+=($!)
done
for pid in "${pids[@]}"; do wait "$pid" || fail=1; done
[ "$fail" -eq 0 ] || { echo "a zero-DHW process failed: see $LOGS"; exit 1; }

julia +release --project=$ENVP Code/sa_pf/04_collect.jl coralpop $BANDS
echo "pilot done"
