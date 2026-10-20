#=
Tests of the shared Sobol'-run helpers (sa_common.jl).

    julia +release --project=Code/sa Code/sa/test/runtests.jl

Needs the synthetic domains (01_build_domains.jl) for the domain-level checks.
=#

using Test
include(joinpath(@__DIR__, "..", "sa_common.jl"))

@testset "Code/sa helpers" begin
    @testset "specs" begin
        @test nrow(dhw_spec()) == 7 && !any(dhw_spec().is_constant)
        @test nrow(multiplier_spec(MULTIPLIER_BOUNDS)) == 13
        @test all(==(MULTIPLIER_BOUNDS), multiplier_spec(MULTIPLIER_BOUNDS).dist_params)
        s = dhw_selector_spec(2:1025)
        @test s.fieldname == [:dhw_scenario] && s.dist_params[1] == (2.0, 1025.0)
        X = dhw_sample(; n=64)
        @test nrow(X) == 64 && X == dhw_sample(; n=64)  # reproducible
        for (n, lo, hi, _) in DHW_FACTORS
            @test all(lo .<= X[!, n] .<= hi)
        end
        Random.seed!(1)
        D = ADRIA.sample(dhw_selector_spec(2:1025), 256)
        @test all(isinteger, D.dhw_scenario) && all(2 .<= D.dhw_scenario .<= 1025)
        @test length(unique(D.dhw_scenario)) > 200                          # spread over the members
        @test factor_block(:dhw_trend_per_decade) == "DHW" && factor_block(:m_k) == "domain" &&
              factor_block(:seed) == "model"
    end

    @testset "multiplier bands (multiplier_bands.toml)" begin
        shipped = multiplier_bands(joinpath(@__DIR__, "..", "multiplier_bands.toml"))
        @test sort(collect(keys(shipped))) == sort(MULTIPLIER_NAMES)
        @test nrow(multiplier_spec(shipped)) == 13
        mktempdir() do dir
            f = joinpath(dir, "b.toml")
            write(f, "default = [0.8, 1.2]
[bands]
m_k = [0.9, 1.1]
m_depth = [1.0, 1.0]
")
            b = multiplier_bands(f)
            @test b[:m_k] == (0.9, 1.1) && b[:m_area] == (0.8, 1.2) && b[:m_depth] == (1.0, 1.0)
            s = multiplier_spec(b)
            i(n) = findfirst(==(n), s.fieldname)
            @test s.is_constant[i(:m_depth)] && s.val[i(:m_depth)] == 1.0      # lo == hi: fixed
            @test !s.is_constant[i(:m_k)] && s.dist_params[i(:m_k)] == (0.9, 1.1)
            Random.seed!(3)
            M = ADRIA.sample(s, 64)
            @test all(M.m_depth .== 1.0) && all(0.9 .<= M.m_k .<= 1.1) && all(0.8 .<= M.m_area .<= 1.2)
            write(f, "[bands]
m_x = [0.5, 1.5]
");          @test_throws ErrorException multiplier_bands(f)
            write(f, "[bands]
m_k = [1.5, 0.5]
");          @test_throws ErrorException multiplier_bands(f)
            write(f, "[bands]
m_k = [0.0, 1.5]
");          @test_throws ErrorException multiplier_bands(f)
            write(f, "[bands]
m_k = [0.5]
");               @test_throws ErrorException multiplier_bands(f)
            write(f, "defaults = [0.5, 1.5]
");               @test_throws ErrorException multiplier_bands(f)
            write(f, "");                                         @test multiplier_bands(f) == Dict(n => MULTIPLIER_BOUNDS for n in MULTIPLIER_NAMES)
        end
        ## a design's saved spec comes back with its Distributions types (06_collect.jl uses it)
        spec = vcat(dhw_selector_spec(2:1025), multiplier_spec(MULTIPLIER_BOUNDS))
        back = spec_from_disk(spec_for_disk(spec))
        @test back.dist == spec.dist && back.dist_params == spec.dist_params && back.fieldname == spec.fieldname
        @test spec_from_disk(spec) == spec
    end

    @testset "connectivity scaling" begin
        C = [0.2 0.05 0.01; 0.04 0.2 0.03; 0.0 0.06 0.2]
        same, n0 = scale_connectivity(C, 1, 1, 1)
        @test same ≈ C && n0 == 0
        D, _ = scale_connectivity(C, 1.5, 0.5, 1.2)
        @test mean(sum(D; dims=2)) ≈ 1.2 * mean(sum(C; dims=2))           # total set by m_retention
        ratio(M) = sum(M[i, i] for i in 1:3) / (sum(M) - sum(M[i, i] for i in 1:3))
        @test ratio(D) ≈ ratio(C) * 1.5 / 0.5                               # split set by self / strength
        D1, _ = scale_connectivity(C, 2.0, 2.0, 1.0)
        @test D1 ≈ C                                                        # equal scaling is undone
        big, ncap = scale_connectivity([0.9 0.3; 0.2 0.9], 1.5, 1.0, 1.5)
        @test ncap > 0 && all(sum(big; dims=2) .<= 1 + 1e-12)
        single, _ = scale_connectivity(fill(0.2, 1, 1), 1.5, 1.3, 0.7)
        @test single[1, 1] ≈ 0.14                                           # one reef: only m_retention
    end

    @testset "initial cover scaling" begin
        cover = zeros(35, 2)
        for g in 1:5
            cover[((g - 1) * 7 + 1):(g * 7), 1] .= g * 0.01 .* (1:7) ./ 28
        end
        cover[1:7, 2] .= 0.05 / 7                                           # site 2: one group only
        same = scale_initial_cover(cover, 1.0, ones(5))
        @test same ≈ cover
        out = scale_initial_cover(cover, 1.4, [2.0, 1, 1, 1, 0.5])
        @test sum(out[:, 1]) ≈ 1.4 * sum(cover[:, 1])
        tot(M, g, l) = sum(M[((g - 1) * 7 + 1):(g * 7), l])
        w = [tot(cover, g, 1) for g in 1:5] .* [2.0, 1, 1, 1, 0.5]
        @test [tot(out, g, 1) for g in 1:5] ./ sum(out[:, 1]) ≈ w ./ sum(w)
        @test out[1:7, 1] ./ sum(out[1:7, 1]) ≈ cover[1:7, 1] ./ sum(cover[1:7, 1])  # size split kept
        @test sum(out[:, 2]) ≈ 1.4 * 0.05 && all(out[8:end, 2] .== 0)
        @test sum(scale_initial_cover(fill(0.5 / 35, 35, 1), 3.0, ones(5))) ≈ 1.0  # capped at k
    end

    if all(d -> isdir(domain_path(d)), SA_DOMAINS)
        @testset "domain multipliers on $d" for d in SA_DOMAINS
            dom = load_sa_domain(d)
            ones_row = Dict(n => 1.0 for n in MULTIPLIER_NAMES)
            d1, info = apply_domain_multipliers(dom, ones_row)
            @test info.n_conn_capped == 0
            @test d1.loc_data.k ≈ dom.loc_data.k && d1.loc_data.area ≈ dom.loc_data.area
            @test d1.loc_data.depth_med ≈ dom.loc_data.depth_med
            @test Array(d1.init_coral_cover) ≈ Array(dom.init_coral_cover)
            @test Array(d1.conn) ≈ Array(dom.conn) && Array(d1.wave_scens) ≈ Array(dom.wave_scens)
            r = merge(ones_row, Dict(:m_k => 1.2, :m_area => 0.5, :m_depth => 1.5, :m_ub => 0.8, :m_retention => 1.3))
            d2, _ = apply_domain_multipliers(dom, r)
            @test d2.loc_data.k ≈ min.(dom.loc_data.k .* 1.2, 1) && d2.loc_data.area ≈ dom.loc_data.area ./ 2
            @test d2.loc_data.depth_med ≈ dom.loc_data.depth_med .* 1.5
            @test Array(d2.wave_scens) ≈ Array(dom.wave_scens) .* 0.8
            @test mean(sum(Array(d2.conn); dims=2)) ≈ 1.3 * mean(sum(Array(dom.conn); dims=2))
            @test Array(dom.conn) ≈ Array(d1.conn)    # the original is untouched
        end
    end
end
