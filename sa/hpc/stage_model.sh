#!/usr/bin/env bash
# stage_model.sh MODEL -- assemble the build folder of one model's Sobol'-run image (run in WSL).
#
#   bash /mnt/c/Users/vhaller/DecisionSupportProjects/CoralModelReview/Code/sa/hpc/stage_model.sh coralpop
#
# Writes (default BUILD_DIR=~/sa-build/MODEL, on WSL's own ext4 disk -- NOT under /mnt/c):
#
#   ~/sa-build/MODEL/
#       MODEL_sa.def                the image definition (Code/sa/models/MODEL/hpc/MODEL_sa.def)
#       sa/                         becomes /opt/sa in the image
#           repos/<repo>            the clones listed in Code/sa/models/MODEL/hpc/repos.txt (code only)
#           Code/sa/                the Code/sa scripts; models/MODEL's [sources] rewritten to ../repos
#           STAGED_SOURCES.md       branch / commit / uncommitted-changes per repo
#           patches/<repo>.patch    the uncommitted changes themselves, for repos that have any
#       upload/jobs/                the Slurm scripts (LF line endings), for the HPC
#
# Override GITHUB_SRC / SA_SRC / BUILD_DIR from the environment.
set -euo pipefail

MODEL="${1:?usage: stage_model.sh MODEL (e.g. coralpop)}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SA_SRC="${SA_SRC:-$(cd "$HERE/.." && pwd)}"
CODE_SRC="$(cd "$SA_SRC/.." && pwd)"
GITHUB_SRC="${GITHUB_SRC:-/mnt/c/Users/vhaller/Documents/GitHub}"
BUILD_DIR="${BUILD_DIR:-$HOME/sa-build/$MODEL}"
MODEL_DIR="$SA_SRC/models/$MODEL"
DEF="$MODEL_DIR/hpc/${MODEL}_sa.def"

for f in "$MODEL_DIR/Project.toml" "$MODEL_DIR/Manifest.toml" "$MODEL_DIR/hpc/repos.txt" "$DEF" "$SA_SRC/models/$MODEL.jl"; do
    [ -f "$f" ] || { echo "ERROR: not found: $f" >&2; exit 1; }
done
mapfile -t REPOS < <(sed -e 's/#.*//' -e 's/[[:space:]]*$//' "$MODEL_DIR/hpc/repos.txt" | grep -v '^$')

EXCLUDES=(--exclude=.git/ --exclude=.claude/ --exclude=.vscode/ --exclude=Outputs/ --exclude=docs/ --exclude=data-raw/
          --exclude=sandbox/ --exclude=.CondaPkg/ --exclude='*.pptx' --exclude='*.ipynb' --exclude='*.docx'
          --exclude=__pycache__/)

# Windows git sees these clones exactly as you do (line endings, file modes); WSL's git may not.
GIT=git; command -v git.exe >/dev/null 2>&1 && GIT=git.exe
gitc() { local repo="$1"; shift; if [ "$GIT" = git.exe ]; then git.exe -C "$(wslpath -w "$repo")" "$@"; else git -C "$repo" "$@"; fi; }

echo "Staging $MODEL into $BUILD_DIR"
mkdir -p "$BUILD_DIR/sa/repos" "$BUILD_DIR/sa/patches" "$BUILD_DIR/upload"
rm -rf "$BUILD_DIR/sa/patches/"*

SRCMD="$BUILD_DIR/sa/STAGED_SOURCES.md"
{
    echo "# Sources staged into ${MODEL}_sa.sif"
    echo
    echo "Staged $(date '+%Y-%m-%d %H:%M') from \`$GITHUB_SRC\` and \`$SA_SRC\`."
    echo "Repos with uncommitted changes have them saved in \`/opt/sa/patches/<repo>.patch\` (tracked files)."
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
        gitc "$path" diff HEAD > "$BUILD_DIR/sa/patches/$name.patch" || true
        gitc "$path" status --porcelain > "$BUILD_DIR/sa/patches/$name.status" || true
    else
        dirty="no"
    fi
    echo "| $name | $branch | ${commit//|/-} | $dirty |" >> "$SRCMD"
}

for r in "${REPOS[@]}"; do
    src="$GITHUB_SRC/$r"
    [ -d "$src" ] || { echo "ERROR: missing clone $src" >&2; exit 1; }
    echo "  repo $r"
    rsync -a --delete "${EXCLUDES[@]}" "$src/" "$BUILD_DIR/sa/repos/$r/"
    record "$r" "$src"
done

echo "  Code/sa"
mkdir -p "$BUILD_DIR/sa/Code"
rsync -a --delete "${EXCLUDES[@]}" "$SA_SRC/" "$BUILD_DIR/sa/Code/sa/"
record "CoralModelReview-Code" "$CODE_SRC"

# [sources] / Manifest paths: "../../../../../../Documents/GitHub/<repo>/..." -> "../../../../repos/<repo>/..."
# (relative to /opt/sa/Code/sa/models/MODEL; a separator written as an escaped backslash is accepted too)
python3 - "$BUILD_DIR/sa/Code/sa/models/$MODEL" "${REPOS[@]}" <<'PY'
import re, sys, pathlib
env, repos = pathlib.Path(sys.argv[1]), set(sys.argv[2:])
SEP = r'(?:/|\\\\)'
pat = re.compile(r'(path\s*=\s*")(?:\.\.' + SEP + r')+Documents' + SEP + r'GitHub' + SEP + r'([^"]*)(")')
for name in ("Project.toml", "Manifest.toml"):
    f = env / name
    s = f.read_text(encoding="utf-8")
    n, missing = 0, set()
    def fix(m):
        global n; n += 1
        rel = re.sub(SEP, "/", m.group(2))
        repo = rel.split("/")[0]
        if repo not in repos:
            missing.add(repo)
        return m.group(1) + "../../../../repos/" + rel + m.group(3)
    s2 = pat.sub(fix, s)
    if missing:
        sys.exit(f"ERROR: {name} path-depends on {sorted(missing)}, which repos.txt does not stage")
    left = re.findall(r'path\s*=\s*"[^"]*Documents[^"]*"', s2)
    if left:
        sys.exit(f"ERROR: {name}: unrewritten paths remain: {left}")
    f.write_text(s2, encoding="utf-8", newline="\n")
    print(f"    {name}: rewrote {n} source paths")
PY

cp -f "$DEF" "$BUILD_DIR/${MODEL}_sa.def"
cp -f "$HERE/build_model_image.sh" "$BUILD_DIR/build_model_image.sh"
rsync -a --delete "$HERE/jobs/" "$BUILD_DIR/upload/jobs/"
find "$BUILD_DIR/${MODEL}_sa.def" "$BUILD_DIR/build_model_image.sh" "$BUILD_DIR/upload/jobs" "$BUILD_DIR/sa/Code/sa/hpc" \
    -type f \( -name '*.sh' -o -name '*.slurm' -o -name '*.def' \) -exec sed -i 's/\r$//' {} +
chmod +x "$BUILD_DIR/build_model_image.sh" "$BUILD_DIR/upload/jobs/"*.sh "$BUILD_DIR/upload/jobs/"*.slurm 2>/dev/null || true

echo
cat "$SRCMD"
echo
du -sh "$BUILD_DIR/sa" "$BUILD_DIR/upload"
echo "Next: bash $BUILD_DIR/build_model_image.sh $MODEL"
