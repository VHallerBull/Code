#!/usr/bin/env bash
# stage_sources.sh -- assemble the container build folder from the local clones (run in WSL).
#
#   bash "/mnt/c/Users/vhaller/OneDrive - Australian Institute of Marine Science/DecisionSupport/CoralModelReview/Code/hpc/stage_sources.sh"
#
# Writes (default BUILD_DIR=~/coral-build, on WSL's own ext4 disk -- NOT under /mnt/c or OneDrive):
#
#   ~/coral-build/
#       coral_engines.def                 the image definition (copied from this folder)
#       coral/                            becomes /opt/coral in the image
#           repos/<14 engine repos>       code only: no .git, sandbox/, Outputs/, docs/
#           Code/                         this repo, Project/Manifest [sources] rewritten to ../repos/<name>
#           data/Test_domain/             the ADRIA Test domain (folder name must stay Test_domain)
#           data/cscape_demog/            C~scape's six converted demog_julia_* files
#           data/calibrated_params.nc     the CoralBlox calibration (CALIB_SRC; required)
#           STAGED_SOURCES.md             branch / commit / uncommitted-changes per repo
#           patches/<repo>.patch          the uncommitted changes themselves, for repos that have any
#       upload/                           what goes to the HPC besides the .sif (see RUNNING_ON_HPC.md)
#           jobs/                         the Slurm scripts, LF line endings
#           reefmod_host/                 ReefMod-GBR + driver, ONLY for the fallback ReefMod route
#
# Override any of the *_SRC / BUILD_DIR variables below from the environment.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE_SRC="${CODE_SRC:-$(cd "$HERE/.." && pwd)}"
GITHUB_SRC="${GITHUB_SRC:-/mnt/c/Users/vhaller/Documents/GitHub}"
DOMAIN_SRC="${DOMAIN_SRC:-$GITHUB_SRC/ADRIA.jl/ADRIA/test/data/Test_domain}"
DEMOG_SRC="${DEMOG_SRC:-$GITHUB_SRC/CscapeInterface.jl/sandbox/ADRIAIntegration/projects/AdjustedRecruitment10_Disturbance/data}"
# The 2026-08-11 CoralBlox calibration. Required: without it ADRIA silently falls back to its built-in
# parameters (every run before 2026-10-12 did, by mistake). engines.jl picks the baked-in copy up by itself.
CALIB_SRC="${CALIB_SRC:-/mnt/c/Users/vhaller/Australian Institute of Marine Science/AIMS Decision Support - General/01_PROJECTS/RRAP M&DS/ADRIA/ADRIA Domain data/CoralBlox-ADRIA GBR Calibrated Parameters/2026-08-11/calibrated_params.nc}"
BUILD_DIR="${BUILD_DIR:-$HOME/coral-build}"

REPOS=(ADRIA.jl ADRIAIndicators.jl CScape.jl CscapeInterface.jl CoralPopInterface.jl CoralPop-Engine
       Kora.jl KoraInterface.jl CORSET.jl CORSETInterface.jl ReefModInterface.jl ReefMod-GBR
       CoCoNet-python CoCoNetInterface.jl)

# Never copied into the image. ReefModInterface keeps its (tiny) sandbox: the fallback ReefMod route
# runs its 02_prepare_hpc_job.jl / 03_load_hpc_results.jl.
EXCLUDES=(--exclude=.git/ --exclude=.claude/ --exclude=.vscode/ --exclude=Outputs/ --exclude=docs/ --exclude=data-raw/
          --exclude='*.pptx' --exclude='*.ipynb' --exclude='*.docx' --exclude=__pycache__/)

# Windows git sees these clones exactly as you do (line endings, file modes); WSL's git may not.
GIT=git; command -v git.exe >/dev/null 2>&1 && GIT=git.exe
gitc() { local repo="$1"; shift; if [ "$GIT" = git.exe ]; then git.exe -C "$(wslpath -w "$repo")" "$@"; else git -C "$repo" "$@"; fi; }

strip_cr() { find "$@" -type f \( -name '*.sh' -o -name '*.slurm' -o -name '*.def' \) -exec sed -i 's/\r$//' {} +; }

for d in "$CODE_SRC" "$GITHUB_SRC" "$DOMAIN_SRC" "$DEMOG_SRC"; do
    [ -d "$d" ] || { echo "ERROR: not found: $d" >&2; exit 1; }
done
[ -f "$CALIB_SRC" ] || { echo "ERROR: CoralBlox calibration not found: $CALIB_SRC (set CALIB_SRC)" >&2; exit 1; }
[ "$(basename "$DOMAIN_SRC")" = "Test_domain" ] || echo "NOTE: domain folder is $(basename "$DOMAIN_SRC"); CMR_DOMAIN_PATH must point at it."

echo "Staging into $BUILD_DIR"
mkdir -p "$BUILD_DIR/coral/repos" "$BUILD_DIR/coral/data" "$BUILD_DIR/coral/patches" "$BUILD_DIR/upload"
rm -rf "$BUILD_DIR/coral/patches/"*

# ---- engine repos, with their provenance ----
SRCMD="$BUILD_DIR/coral/STAGED_SOURCES.md"
{
    echo "# Sources staged into this image"
    echo
    echo "Staged $(date '+%Y-%m-%d %H:%M') from \`$GITHUB_SRC\` and \`$CODE_SRC\`."
    echo "Repos with uncommitted changes have them saved in \`/opt/coral/patches/<repo>.patch\` (tracked files)."
    echo
    echo "| Repo | Branch | Commit | Uncommitted changes |"
    echo "|---|---|---|---|"
} > "$SRCMD"

record() {   # record <name> <path>
    local name="$1" path="$2" branch commit dirty
    branch="$(gitc "$path" rev-parse --abbrev-ref HEAD 2>/dev/null || echo '?')"
    commit="$(gitc "$path" log -1 --format='%h %s' 2>/dev/null || echo '?')"
    if [ -n "$(gitc "$path" status --porcelain 2>/dev/null)" ]; then
        dirty="**yes**"
        gitc "$path" diff HEAD > "$BUILD_DIR/coral/patches/$name.patch" || true
        gitc "$path" status --porcelain >> "$BUILD_DIR/coral/patches/$name.status" || true
    else
        dirty="no"
    fi
    echo "| $name | $branch | ${commit//|/-} | $dirty |" >> "$SRCMD"
}

for r in "${REPOS[@]}"; do
    src="$GITHUB_SRC/$r"
    [ -d "$src" ] || { echo "ERROR: missing clone $src" >&2; exit 1; }
    echo "  repo $r"
    extra=(--exclude=sandbox/)
    [ "$r" = "ReefModInterface.jl" ] && extra=()
    # CoCoNet: of its 166 MB legacy/ folder (NetLogo model, GBR reef table) only the coastline is read
    # (CMR_COCONET_KM_OFFSHORE=coastline); its JS viewers aren't needed either.
    [ "$r" = "CoCoNet-python" ] && extra+=(--include=legacy/coastline.csv --exclude='legacy/*' --exclude=viz/ --exclude=gen-charts/ --exclude=documentation/)
    rsync -a --delete "${EXCLUDES[@]}" "${extra[@]}" "$src/" "$BUILD_DIR/coral/repos/$r/"
    record "$r" "$src"
done

# ---- Code/ (this repo), with [sources] pointing at ../repos ----
echo "  Code"
rsync -a --delete "${EXCLUDES[@]}" --exclude=.CondaPkg/ "$CODE_SRC/" "$BUILD_DIR/coral/Code/"
record "CoralModelReview-Code" "$CODE_SRC"
python3 - "$BUILD_DIR/coral/Code" <<'PY'
import re, sys, pathlib
code = pathlib.Path(sys.argv[1])
# path = "../../../../Documents/GitHub/<repo>/..."  ->  "../repos/<repo>/..."  (in both files; a
# separator written as an escaped Windows backslash, "\\", is accepted too)
SEP = r'(?:/|\\\\)'
pat = re.compile(r'(path\s*=\s*")(?:\.\.' + SEP + r')+Documents' + SEP + r'GitHub' + SEP + r'([^"]*)(")')
for name in ("Project.toml", "Manifest.toml"):
    f = code / name
    s = f.read_text(encoding="utf-8")
    n = 0
    def fix(m):
        global n; n += 1
        return m.group(1) + "../repos/" + re.sub(SEP, "/", m.group(2)) + m.group(3)
    s2 = pat.sub(fix, s)
    left = re.findall(r'path\s*=\s*"[^"]*Documents[^"]*"', s2)
    if left:
        sys.exit(f"ERROR: {name}: unrewritten paths remain: {left}")
    f.write_text(s2, encoding="utf-8", newline="\n")
    print(f"    {name}: rewrote {n} source paths")
PY

# ---- data ----
echo "  data"
rsync -a --delete "$DOMAIN_SRC/" "$BUILD_DIR/coral/data/$(basename "$DOMAIN_SRC")/"
mkdir -p "$BUILD_DIR/coral/data/cscape_demog"
cp -f "$DEMOG_SRC"/demog_julia_*.rds "$BUILD_DIR/coral/data/cscape_demog/"
n_demog=$(ls "$BUILD_DIR/coral/data/cscape_demog"/demog_julia_*.rds | wc -l)
[ "$n_demog" -ge 6 ] || { echo "ERROR: expected 6 demog_julia_*.rds in $DEMOG_SRC, found $n_demog" >&2; exit 1; }
cp -f "$CALIB_SRC" "$BUILD_DIR/coral/data/calibrated_params.nc"
echo "    CoralBlox calibration staged ($(sha256sum "$BUILD_DIR/coral/data/calibrated_params.nc" | cut -c1-12)...; engines.jl uses it by default)"
{
    echo
    echo "Data: domain \`$(basename "$DOMAIN_SRC")\` from \`$DOMAIN_SRC\`; $n_demog C~scape demography files from \`$DEMOG_SRC\`."
    echo "CoralBlox calibration: \`$CALIB_SRC\` (sha256 $(sha256sum "$BUILD_DIR/coral/data/calibrated_params.nc" | cut -c1-16)...)."
} >> "$SRCMD"

# ---- the definition file and the HPC upload folder ----
cp -f "$HERE/coral_engines.def" "$BUILD_DIR/coral_engines.def"
rsync -a --delete "$HERE/jobs/" "$BUILD_DIR/upload/jobs/"
mkdir -p "$BUILD_DIR/upload/reefmod_host"
rsync -a --delete "${EXCLUDES[@]}" "$GITHUB_SRC/ReefMod-GBR/" "$BUILD_DIR/upload/reefmod_host/ReefMod-GBR/"
rsync -a --delete "$GITHUB_SRC/ReefModInterface.jl/driver/" "$BUILD_DIR/upload/reefmod_host/driver/"
strip_cr "$BUILD_DIR/coral_engines.def" "$BUILD_DIR/coral/Code/hpc" "$BUILD_DIR/upload/jobs"
chmod +x "$BUILD_DIR/upload/jobs/"*.sh "$BUILD_DIR/upload/jobs/"*.slurm 2>/dev/null || true

echo
cat "$SRCMD"
echo
du -sh "$BUILD_DIR/coral" "$BUILD_DIR/upload"
echo "Next: bash \"$HERE/build_image.sh\""
