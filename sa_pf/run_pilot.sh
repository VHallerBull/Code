#!/usr/bin/env bash
# Local pick-and-freeze pilot: the designs, then N_PROCS parallel Julia processes (one process per share:
# CoralPop runs through Python, no threads), then collect. Run from the CoralModelReview folder:
#
#   bash Code/sa_pf/run_pilot.sh [n_base=256] [n_procs=12] [sets="main zero"]
#
# MODEL (default coralpop), BANDS (default "0.1 0.2 0.5") and DOMAINS (default both) are environment
# variables, e.g. `MODEL=cscape BANDS=0.5 DOMAINS=SingleReef_SobolDHW bash Code/sa_pf/run_pilot.sh 16 4`.
# Sets: `main` (with the DHW inputs, d + 1 runs per base sample: 177 CoralPop, 84 C~scape) and `zero` (no
# heat stress, no DHW inputs: 169 / 76). The zero set's walk end points are the zero-DHW reference of the
# output figures; without the zero set, 03b_run_zero.jl runs that reference (2 runs per base sample)
# instead. A rerun resumes: finished chunk files are skipped.
set -euo pipefail
N_BASE=${1:-256}
N_PROCS=${2:-12}
SETS=${3:-"main zero"}
MODEL=${MODEL:-coralpop}
BANDS=${BANDS:-"0.1 0.2 0.5"}
DOMAINS=${DOMAINS:-"SingleReef_SobolDHW NineReefs_SobolDHW"}
case "$MODEL" in
  coralpop) FOLDER=CoralPOP; ENVP=Code/sa/shapley ;;
  cscape) FOLDER=CScape; ENVP=Code/sa/models/cscape ;;
  *) echo "run_pilot.sh: unknown MODEL $MODEL"; exit 1 ;;
esac
LOGS=$FOLDER/model_run_pf_pilot/logs
mkdir -p "$LOGS"

julia +release --project=$ENVP Code/sa_pf/02_design.jl "$MODEL" "$N_BASE" $BANDS --sets="${SETS// /,}"

JOBS=""
MAIN_JOBS=""
for set in $SETS; do
  for dom in $DOMAINS; do
    for b in $BANDS; do
      if [ "$set" = main ]; then JOBS="$JOBS $b:$dom"; MAIN_JOBS="$MAIN_JOBS $b:$dom"; else JOBS="$JOBS $b:$dom:$set"; fi
    done
  done
done

run_parallel() {  # run_parallel <script> <log prefix> <jobs...>
  local script=$1 prefix=$2; shift 2
  local pids=() fail=0
  for p in $(seq 1 "$N_PROCS"); do
    julia +release --project=$ENVP "Code/sa_pf/$script" "$MODEL" "$p" "$N_PROCS" "$@" \
      > "$LOGS/${prefix}_$(printf %02d "$p").log" 2>&1 &
    pids+=($!)
  done
  for pid in "${pids[@]}"; do wait "$pid" || fail=1; done
  [ "$fail" -eq 0 ] || { echo "a $script process failed: see $LOGS/${prefix}_*.log"; exit 1; }
}

run_parallel 03_run_chunk.jl proc $JOBS
if [[ " $SETS " != *" zero "* && -n "$MAIN_JOBS" ]]; then
  run_parallel 03b_run_zero.jl zero $MAIN_JOBS     # zero-DHW reference only (no zero design)
fi

julia +release --project=$ENVP Code/sa_pf/04_collect.jl "$MODEL" $BANDS
echo "pilot done: next 05_analyse.jl $MODEL $BANDS, then 06_figures.jl $MODEL"
