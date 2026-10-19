#!/bin/bash
# batch_status.sh -- how every project under runs/ is doing: each job's Slurm state (from jobs.txt and
# sacct), then, once the analysis has written FIGURES.md, the scenarios each engine completed.
#
#   cd ~/CoralModelReview
#   bash jobs/batch_status.sh                          # every project
#   bash jobs/batch_status.sh SingleReef NineReefsLine # only these domains' (or projects') folders
#   bash jobs/batch_status.sh --since 3h               # only projects (re)submitted in the last 3 hours
#   bash jobs/batch_status.sh --since 2d SingleReef    # ... combined with a domain list
#
# --since N[m|h|d] (minutes, hours or days; a bare number means hours) keeps the projects whose jobs.txt was
# written within that time, i.e. the ones submit_all.sh / submit_batch.sh last touched, including a partial
# re-submit of one engine.
#
# States: ok = COMPLETED, run = RUNNING, wait = PENDING, FAIL / OOM (out of memory) / TIMEOUT /
# CANCEL = needs a re-run (6d), ? = not in sacct (too old, or not submitted). An engine job can be "ok"
# with fewer scenarios than requested (ReefMod skips failed scenarios), hence the second line.
set -euo pipefail
source "${CMR_BASE:-$HOME/CoralModelReview}/jobs/common.sh"
cd "$RUNS"
since=""; roots=()
while [ $# -gt 0 ]; do
    case "$1" in
        --since) [ $# -ge 2 ] || { echo "--since needs a value, e.g. 3h" >&2; exit 2; }; since=$2; shift 2 ;;
        --since=*) since=${1#--since=}; shift ;;
        -h|--help) sed -n '2,17p' "$0"; exit 0 ;;
        *) roots+=("$1"); shift ;;
    esac
done
[ ${#roots[@]} -gt 0 ] || roots=(.)
mins=""
if [ -n "$since" ]; then
    [[ $since =~ ^[0-9]+$ ]] && since=${since}h
    [[ $since =~ ^([0-9]+)([mhd]?)$ ]] || { echo "--since: use a number with m, h or d, e.g. 90m, 3h, 2d" >&2; exit 2; }
    case "${BASH_REMATCH[2]:-h}" in
        m) mins=${BASH_REMATCH[1]} ;;
        h) mins=$((BASH_REMATCH[1] * 60)) ;;
        d) mins=$((BASH_REMATCH[1] * 1440)) ;;
    esac
fi

short() {
    case "$1" in
        COMPLETED) echo ok ;; RUNNING) echo run ;; PENDING) echo wait ;; OUT_OF_MEMORY) echo OOM ;;
        CANCELLED) echo CANCEL ;; "") echo "?" ;; *) echo "$1" ;;
    esac
}

n=0; n_done=0; n_bad=0
while read -r f; do
    p=$(dirname "${f#./}")
    ids=$(awk '{print $2}' "$f" | paste -sd, -)
    [ -n "$ids" ] || continue
    n=$((n + 1))
    unset st; declare -A st=()
    # one state per job; array jobs (reefmod_B) report per task: keep any task that is not COMPLETED
    while IFS='|' read -r id state; do
        k=${id%%_*}; state=${state%% *}
        [ -z "${st[$k]:-}" ] || [ "${st[$k]}" = COMPLETED ] && st[$k]=$state
    done < <(sacct -X -n -P -j "$ids" --format=JobID,State 2>/dev/null || true)
    line=""; bad=0; busy=0
    while read -r name id; do
        s=$(short "${st[$id]:-}")
        line+="$name $s  "
        case "$s" in ok|"?") ;; run|wait) busy=1 ;; *) bad=1 ;; esac
    done < "$f"
    if [ -f "$p/FIGURES.md" ] && [ "$busy" = 0 ]; then
        status=DONE; n_done=$((n_done + 1))
    elif [ "$busy" = 1 ]; then status=busy
    else status="NO FIGURES"; fi
    [ "$bad" = 1 ] && { status="$status, CHECK"; n_bad=$((n_bad + 1)); }
    printf '%-26s %-20s %s\n' "$p" "[$status]" "$line"
    if [ -f "$p/FIGURES.md" ]; then
        sc=$(grep -m1 '^- Scenarios:' "$p/FIGURES.md" | sed 's/^- Scenarios: //') || true
        miss=$(grep -m1 'Not included' "$p/FIGURES.md" | sed 's/^- //; s/\*\*//g') || true
        nfail=$(awk '/^## Figures that failed/ { f = 1; next } f && /^- / { n++ } END { print n + 0 }' "$p/FIGURES.md")
        printf '%-47s %s%s%s\n' "" "scenarios: ${sc:-?}" "${miss:+; $miss}" \
            "$( [ "${nfail:-0}" -gt 0 ] && echo "; $nfail failed figure(s)" )"
    fi
done < <(find "${roots[@]}" -name jobs.txt ${mins:+-mmin -$mins} | sort)

echo "$n projects${since:+ submitted in the last $since}: $n_done done, $n_bad with a failed/cancelled job (re-run with: bash jobs/submit_all.sh <project> --engines <engine>, same CMR_* settings)"
