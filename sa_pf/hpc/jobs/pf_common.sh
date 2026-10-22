# pf_common.sh -- sourced by the pick-and-freeze job scripts in this folder (not run on its own).
#
# Same layout and image conventions as Code/sa/hpc/jobs/sa_common.sh (sourced below), plus:
#   ~/CoralModelReview/sa/
#       images/<model>_pf.sif               the model's image with Code/sa_pf inside (stage_pf.sh)
#       pf_inputs/<domain>/dhw_draws.jld2   the DHW generator's fixed random draws (01_dhw_draws.jl)
#       <ModelFolder>/model_run_pf/band<NN>/design_pf[_zero].{jld2,toml}   designs in
#       <ModelFolder>/model_run_pf/band<NN>/<domain>[/pf_zero_dhw]/        chunks and pf_results.jld2 out

SA_BASE="${SA_BASE:-$HOME/CoralModelReview/sa}"
MODEL="${MODEL:?set MODEL (e.g. coralpop)}"
SIF="${SIF:-$SA_BASE/images/${MODEL}_pf.sif}"
source "$SA_BASE/jobs/sa_common.sh"
PF_CODE=/opt/sa/Code/sa_pf                   # the scripts inside the image
export CMR_PF_SUBDIR="${CMR_PF_SUBDIR:-model_run_pf}"
export CMR_PF_INPUTS_DIR="$SA_BASE/pf_inputs"

# pf_band_dir <band> -> $SA_BASE/<ModelFolder>/model_run_pf/band<NN> (band 0.5 -> band50)
pf_band_dir() {
    local nn
    nn=$(awk -v b="$1" 'BEGIN { printf "%02d", b * 100 + 0.5 }')
    echo "$SA_BASE/$(model_folder "$MODEL")/$CMR_PF_SUBDIR/band$nn"
}
