#=
reefmod_fallback.jl prepare|finish -- the fallback ReefMod route, for when the host MATLAB cannot run
inside the container: Julia stages here (in the container), MATLAB on the host in between.

    prepare   (jobs/reefmod_A.slurm)  sample ReefMod's factors and write one MATLAB input file per
                                      scenario into <project>/reefmod/run/inputs/
    [MATLAB]  (jobs/reefmod_B.slurm)  one Slurm array task per scenario: run_reefmod_scenario.m
    finish    (jobs/reefmod_C.slurm)  read the MATLAB outputs into a ReefModResultSet and write the
                                      SAME <project>/reefmod/ layout 08_RUN_ENGINE.jl writes

The calls are the ones ADRIA's :reefmod engine makes (ADRIA/ext/ADRIAReefModExt.jl) with the stage
boundary in the middle, and the settings come from the same `engine_config()`, so 09 cannot tell the
two routes apart -- except in the benchmarks, where this route's run time is the sum of the array
tasks' wall times (read from their `time_reefmod_B_<task>.log` files) instead of one process's.
=#

include(joinpath(@__DIR__, "..", "engines.jl"))

stage = isempty(ARGS) ? "" : ARGS[1]
stage in ("prepare", "finish") || error("usage: julia reefmod_fallback.jl prepare|finish")

cfg = engine_config()
const OUT = engine_dir(cfg.project, :reefmod)
const JOB_DIR = joinpath(OUT, "run")
dom, dom_fixed, YEARS = load_shared_domain(cfg)

if stage == "prepare"
    t = @elapsed begin
        isdir(JOB_DIR) && rm(JOB_DIR; recursive=true, force=true)
        for f in ("result.jld2", "samples.jld2", "benchmark.csv", "size_structure.jld2", "thermal_tolerance.jld2")
            isfile(joinpath(OUT, f)) && rm(joinpath(OUT, f))
        end
        spec = make_spec(:reefmod, dom_fixed, cfg, YEARS, JOB_DIR)
        Random.seed!(cfg.seed)
        scens = ADRIA.sample(spec, cfg.n_samples)
        inp = domain_to_reefmod_inputs(dom; dhw_idx=cfg.dhw_idx)
        scenarios, _, _ = reefmod_sample_to_scenarios(scens, 1)
        paths = prepare_reefmod_job(inp, scenarios, JOB_DIR; domain_path=cfg.domain_path,
            rcp=parse(Int, cfg.rcp), dhw_idx=cfg.dhw_idx, spec)
        JLD2.save(joinpath(OUT, "samples.jld2"), Dict(
            "scens" => scens, "varying" => String.(spec.fieldname[.!spec.is_constant]),
            "spec" => spec_for_disk(spec), "settings" => settings_dict(cfg, YEARS)))
    end
    write(joinpath(OUT, "prepare_seconds.txt"), string(t))
    @info "ReefMod inputs written" n = length(paths) dir = joinpath(JOB_DIR, "inputs")
else
    s = JLD2.load(joinpath(OUT, "samples.jld2"))
    meta = JLD2.load(joinpath(JOB_DIR, "job_meta.jld2"))
    n_out = length(reefmod_output_paths(JOB_DIR))
    n_out == nrow(meta["scenarios"]) ||
        @warn "Only $n_out of $(nrow(meta["scenarios"])) ReefMod scenarios produced output; continuing with those."
    rs, t_load = timed_call(() -> load_reefmod_results(JOB_DIR; scenarios=meta["scenarios"], dom,
        spec=meta["spec"], name="reefmod"))
    save_engine_result(:reefmod, rs, OUT)
    ss = extract_size_structure(:reefmod, rs, JOB_DIR, dom, cfg, s["scens"])
    ss === nothing || JLD2.save(joinpath(OUT, "size_structure.jld2"), ss)
    tt = extract_thermal_tolerance(:reefmod, rs, JOB_DIR, dom, cfg, s["scens"])
    tt === nothing || JLD2.save(joinpath(OUT, "thermal_tolerance.jld2"), tt)

    ## per-task wall time / peak memory from GNU time's verbose logs written by reefmod_B.slurm
    task_logs = filter(f -> startswith(f, "time_reefmod_B_"), readdir(dirname(OUT)))
    secs, rss = Float64[], Float64[]
    for f in task_logs
        txt = read(joinpath(dirname(OUT), f), String)
        m = match(r"Elapsed \(wall clock\) time.*: (?:(\d+):)?(\d+):([\d.]+)", txt)
        m === nothing || push!(secs, (m[1] === nothing ? 0 : parse(Int, m[1])) * 3600 + parse(Int, m[2]) * 60 + parse(Float64, m[3]))
        r = match(r"Maximum resident set size \(kbytes\): (\d+)", txt)
        r === nothing || push!(rss, parse(Float64, r[1]) / 1e6)
    end
    n = nrow(rs.inputs)
    t_main = isempty(secs) ? NaN : sum(secs)
    bench = (
        engine="reefmod", label="ReefMod", host=gethostname(), date=string(now()),
        n_scenarios=n, n_locations=length(rs.loc_ids), n_years=length(YEARS),
        threads=1, n_workers=length(secs), cpus=1,
        slurm_job_id=get(ENV, "REEFMOD_STAGE_B_JOB", ""), slurm_cpus="1",
        t_load_s=NaN, t_domain_s=parse(Float64, read(joinpath(OUT, "prepare_seconds.txt"), String)),
        t_first_s=isempty(secs) ? NaN : first(secs), t_first_compile_s=NaN, t_first_alloc_gb=NaN,
        t_main_s=t_main, t_main_compile_s=NaN, t_main_gc_s=NaN, t_main_alloc_gb=NaN,
        t_per_scen_s=t_main / n, t_per_scen_site_year_ms=1000 * t_main / (n * length(rs.loc_ids) * length(YEARS)),
        t_size_extract_s=NaN, t_thermal_extract_s=NaN, maxrss_gb=isempty(rss) ? NaN : maximum(rss),
        disk_gb=dir_size_gb(JOB_DIR), disk_per_scen_gb=dir_size_gb(JOB_DIR) / n, t_total_s=NaN,
    )
    write_csv(joinpath(OUT, "benchmark.csv"), [bench])
    @info "ReefMod (fallback route) result written" OUT n t_main
end
