#=
Tests of the C~scape side of the pick-and-freeze design (models/cscape.jl, adapters/cscape.jl). Separate
from runtests.jl (CoralPop): one process loads one model.

    julia +release --project=Code/sa/models/cscape Code/sa_pf/test/runtests_cscape.jl
=#

using Test, Statistics
const MODEL = "cscape"
include(joinpath(@__DIR__, "..", "pf_common.jl"))
include(joinpath(@__DIR__, "..", "..", "sa", "shapley", "groups.jl"))

@testset "C~scape pick-and-freeze" begin
    spec = cscape_native_spec(; band=0.5)
    D = pf_design(0.5, 4)

    @testset "factors as CScape_sampled_parameters.md" begin
        @test nrow(spec) == 62 && !any(spec.is_constant)
        @test count(startswith("draw_"), String.(spec.fieldname)) == 18
        @test ncol(D.A) == 83 && ncol(pf_design(0.5, 4; set="zero").A) == 75
        r(n) = spec[findfirst(==(n), spec.fieldname), :]
        @test r(:draw_surv_acro_table).dist_params == (1.0, 8000.0)
        @test r(:fert_prob_brooder).dist_params == (0.5, 1.0)                       # capped at 1
        @test all(r(:bleaching_slope).dist_params .≈ (0.347 * 0.8, 0.347 * 1.2))   # ±20 %
        @test all(r(:recruitment_multiplier).dist_params .≈ (5.0, 15.0))
        @test all(cscape_native_spec(; band=0.1)[findfirst(==(:bleaching_slope), spec.fieldname), :].dist_params .≈
                  (0.347 * 0.9, 0.347 * 1.1))                                       # min(band, 0.2)
        ## every input in exactly one group of every scheme (errors otherwise)
        G = load_groupings(joinpath(@__DIR__, "..", "..", "sa", "shapley", "groups", "cscape.toml"), names(D.A))
        @test Set(keys(G)) ⊇ Set(["block", "process", "source"])
    end

    dom = load_sa_domain("SingleReef_SobolDHW")
    draws = load_dhw_draws(joinpath(PF_INPUTS, "SingleReef_SobolDHW", "dhw_draws.jld2"))
    row = walk_rows(D, 1)[1, :]

    @testset "inputs reach C~scape" begin
        traj = generate_dhw(row, Int(row[:dhw_weather]), draws)
        dom.dhw_scens.data[:, :, PF_DHW_SLOT] .= traj
        d, _ = apply_domain_multipliers(dom, row)
        inp = domain_to_cscape_inputs(d; dhw_idx=PF_DHW_SLOT, brood_site_value=Float64(row[:brooder_retention]))
        T = size(inp.temporal_override.disturbances, 1)
        @test inp.temporal_override.disturbances[:, :, CScape.Layout.DIST_DHW_SITE] == traj[1:T, :]
        @test inp.connectivity_override.matrices[1][1, 1, end] ≈ 0.9 * row[:brooder_retention]
        @test inp.spatial_override.ub_med ≈ Float64.(dom.loc_data.ub_med) .* row[:m_ub]   # m_ub reaches ub_med
        dem = cscape_row_demography(row, inp.spatial_override)
        i = findfirst(==("acro_corym"), CSCAPE_SA_FTS)
        @test dem[i].fert_prob == row[:fert_prob_acro_corym] && dem[i].juvi_cap == row[:juvi_cap_acro_corym]
        @test dem[i].distri_setter_cont == (row[:size1yr_mean_acro_corym], row[:size1yr_sd_acro_corym])
        cfg = CScape.build_scenario_input(cscape_scenario_table(row), CSCAPE_SA_SCENARIO_ID, cscape_project_dir(),
            CSCAPE_NULL_LOG, 1; use_cached_ipm=false, write_interventions=false, seed=CSCAPE_SA_SEED)
        @test cfg.bleaching_coefs == Tuple(Float64(row[k]) for k in CSCAPE_BLEACHING_FACTORS)
        @test cfg.recruitment_multiplier == row[:recruitment_multiplier] && cfg.plasticity == 0.0
        @test CScape.stsp_num(cfg.heritability) == [row[:heritability], 0.0]
        @test CScape.stsp_num(cfg.heat_init) == [0.0, row[:heat_tolerance_sd]]
    end

    @testset "runs: deterministic, outputs in range" begin
        r1 = pf_run(dom, row, draws)
        r2 = pf_run(dom, row, draws)
        @test isequal(r1.state, r2.state)       # fixed seed, heritability SD 0 (isequal: NaN where no cover)
        St = r1.state
        L = size(St["cover"], 2)
        @test size(St["cover"]) == (26, L) && size(St["class_cover"]) == (26, 3, L)
        @test size(St["thermal_mean"]) == (26, 6) && length(St["weights"]) == L
        @test all(0 .<= St["cover"] .<= 1)
        @test maximum(abs.(dropdims(sum(St["class_cover"]; dims=2); dims=2) .- St["cover"])) < 1e-12
        M, W = St["thermal_mean"], St["thermal_weight"]
        @test all(-5 .<= M[W .> 0] .<= 8)
        Yr = derive_outputs(St; thermal_groups=PF_THERMAL_GROUPS)
        @test all(0 .<= Yr["acropora_share"] .<= 1) && all(0 .<= Yr["massive_share"] .<= 1)
        @test all(k -> haskey(Yr, "thermal_$k"), CSCAPE_SA_FTS) && !haskey(Yr, "spatial_cv")
        @test all(isfinite, Yr["thermal_all"]) && all(-5 .<= Yr["thermal_all"] .<= 8)   # defined while coral is left
        z = pf_run_zero(dom, zero_rows(D, 1)[1, :])
        Yz = derive_outputs(z.state; thermal_groups=PF_THERMAL_GROUPS)
        @test all(isfinite, Yz["relative_cover"]) && all(0 .<= Yz["relative_cover"] .<= 1)
    end
end
