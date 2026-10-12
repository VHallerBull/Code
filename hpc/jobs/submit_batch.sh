#!/bin/bash
# submit_batch.sh -- submit_all.sh for every domain x DHW member: one project per pair, written to
# $RUNS/<domain>/DHW<member>/, so a domain's climate members sit side by side.
#
#   cd ~/CoralModelReview
#   bash jobs/submit_batch.sh --domains "SingleReef ThreeReefsLine" --dhw "1 2 17 46" --dry-run
#   CMR_N_SAMPLES=8 CMR_WAVE_IDX=1 CMR_CYCLONE_IDX=1 bash jobs/submit_batch.sh --domains "SingleReef ThreeReefsLine" --dhw "1 2 17 46"
#
# Options:
#   --domains "A B"   domain folder names under $DOMAINS_DIR (default ~/CoralModelReview/domains)
#   --dhw "1 2 17"    DHW members; each becomes its own project <domain>/DHW<member>
#   --dry-run         print the submit_all.sh commands, submit nothing
# Anything else (--engines "...", --reefmod-fallback, --no-analyse) is passed on to submit_all.sh.
#
# Every CMR_* variable in your shell applies to every project (CMR_N_SAMPLES, CMR_WAVE_IDX,
# CMR_CYCLONE_IDX, ...); CMR_DOMAIN_PATH and CMR_DHW_IDX are set per project here. Re-submitting a
# project replaces each engine's previous output in it (as with submit_all.sh).
set -euo pipefail
source "${CMR_BASE:-$HOME/CoralModelReview}/jobs/common.sh"
DOMAINS_DIR="${DOMAINS_DIR:-$CMR_BASE/domains}"

usage='usage: submit_batch.sh --domains "A B" --dhw "1 2 3" [--dry-run] [submit_all.sh options]'
DOMAINS=""
DHWS=""
DRY=0
PASS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --domains) DOMAINS="$2"; shift 2 ;;
        --dhw) DHWS="$2"; shift 2 ;;
        --dry-run) DRY=1; shift ;;
        --engines) PASS+=("$1" "$2"); shift 2 ;;
        *) PASS+=("$1"); shift ;;
    esac
done
[ -n "$DOMAINS" ] && [ -n "$DHWS" ] || { echo "$usage" >&2; exit 1; }
for d in $DOMAINS; do
    [ -d "$DOMAINS_DIR/$d" ] || { echo "no domain folder $DOMAINS_DIR/$d" >&2; exit 1; }
done
for m in $DHWS; do
    [[ "$m" =~ ^[0-9]+$ ]] || { echo "DHW member '$m' is not a number" >&2; exit 1; }
done

n=0
for d in $DOMAINS; do
    for m in $DHWS; do
        project="$d/DHW$m"
        n=$((n + 1))
        if [ "$DRY" = "1" ]; then
            echo "CMR_DOMAIN_PATH=$DOMAINS_DIR/$d CMR_DHW_IDX=$m bash $JOBS_DIR/submit_all.sh $project ${PASS[*]+${PASS[*]}}"
            continue
        fi
        echo "=== $project"
        # ${PASS[@]+...}: an empty array is "unbound" under set -u in bash < 4.4
        CMR_DOMAIN_PATH="$DOMAINS_DIR/$d" CMR_DHW_IDX="$m" \
            bash "$JOBS_DIR/submit_all.sh" "$project" ${PASS[@]+"${PASS[@]}"}
    done
done
[ "$DRY" = "1" ] && echo "($n projects; nothing submitted)" || echo "Submitted $n projects. Results: $RUNS/<domain>/DHW<member>/"
