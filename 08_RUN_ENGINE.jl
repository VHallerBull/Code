#=
08_RUN_ENGINE.jl -- run ONE coral engine on the shared ADRIA domain and save everything `09` needs.

    julia +1.12 --project=. 08_RUN_ENGINE.jl <engine>        # engine: coralblox cscape coralpop kora corset reefmod coconet

One engine per invocation, so each can be its own HPC job (see `hpc/RUNNING_ON_HPC.md`): the seven
write to separate `ModelRuns/<project>/<engine>/` folders and never touch each other's files. The
per-engine code (factor specs, the run call, saving, size-structure extraction) is in `engines.jl`
and is lifted from 06/07, so an engine runs here exactly as it does there.

Configuration comes from `CMR_<NAME>` environment variables -- see `engine_config` in `engines.jl`.
The ones you are most likely to change:

    CMR_PROJECT_NAME   output folder under ModelRuns/ (default TestDomain_N8)
    CMR_DOMAIN_PATH    ADRIA domain (default: the ADRIA Test domain, container path first)
    CMR_DHW_IDX / CMR_WAVE_IDX / CMR_CYCLONE_IDX   shared climate member (default 40 / 16 / 1)
    CMR_N_SAMPLES      scenarios per engine; a power of 2 (default 8)
    CMR_SEED           seeds ADRIA's (scrambled, so random) Sobol' sampler (default 1234)
    CMR_N_WORKERS      parallelism, engine-specific meaning (default 0 = serial; see engines.jl)
    CMR_WARMUP         run 1 scenario first to separate compile/start-up time (default true)
    CMR_MODEL_RUNS_DIR where ModelRuns/ is (the HPC jobs point it outside the read-only image)

ReefMod additionally needs `matlab` reachable (`REEFMOD_MATLAB_EXE`, default `matlab` on PATH).
CoCoNet: CMR_COCONET_KM_OFFSHORE (default 50, or `coastline`), CMR_COCONET_OCEAN_ACIDIFICATION and
CMR_COCONET_COTS (both default false); the local CoCoNet clone is found through COCONET_DIR.

Benchmarks (`<engine>/benchmark.csv`, one row):
    t_load_s        loading the packages (from process start)
    t_domain_s      loading the domain and sampling
    t_first_*       the 1-scenario warm-up: time to first result, including compilation and
                    worker/MATLAB start-up (only if CMR_WARMUP=true)
    t_main_*        the timed N-scenario run; t_per_scen_s = t_main_s / N
    maxrss_gb       peak resident memory of THIS Julia process (workers/subprocesses not included --
                    use the Slurm accounting for those, see hpc/jobs/collect_benchmarks.sh)
    disk_gb         size of the engine's output folder
=#

const T_PROCESS_START = time()

include(joinpath(@__DIR__, "engines.jl"))       # sets interpreter env vars, loads every package

const T_LOADED = time()

isempty(ARGS) && error("usage: julia --project=. 08_RUN_ENGINE.jl <engine>, engine one of $(ENGINE_ORDER)")
const KEY = Symbol(lowercase(ARGS[1]))
KEY in ENGINE_ORDER || error("Unknown engine '$(ARGS[1])', expected one of $(ENGINE_ORDER)")
const LABEL = ENGINE_INFO[KEY].label

cfg = engine_config()
@info "Running $LABEL" cfg threads = Threads.nthreads()
ispow2(cfg.n_samples) || @warn "N_SAMPLES=$(cfg.n_samples) is not a power of 2; ADRIA's Sobol' sampler may refuse it."

## C~scape only: check its demography files against brms's own predictions BEFORE anything runs (a file
## whose link/offset/formula is wrong stops the job here), and record the outcome with the results.
const CSCAPE_DEMOGRAPHY = KEY === :cscape ? cscape_demography_summary(cfg) : nothing
KEY === :cscape && @info "C~scape demography: $CSCAPE_DEMOGRAPHY"

# =============================================================================
# Output folder -- a re-run replaces this engine's previous output (and only this engine's)
# =============================================================================

"""
    remove_dir(d; move_aside=false)

`rm(d; recursive=true)`, tolerating NFS (the HPC home): a file some handle still has open survives there
as a `.nfs*` placeholder, so the folder cannot be removed (ENOTEMPTY -- seen right after CoralBlox's
warm-up, whose discarded result set still held its files). The garbage collector runs first so a dropped
result set closes them; if the folder still cannot be removed it is left in place with a warning, or with
`move_aside=true` renamed out of the way so the run starts from an empty folder.
"""
function remove_dir(d::AbstractString; move_aside::Bool=false)
    isdir(d) || return nothing
    GC.gc()
    try
        rm(d; recursive=true, force=true)
    catch err
        if move_aside
            aside = "$(d)_stale_$(Dates.format(now(), "yyyymmdd_HHMMSS"))"
            @warn "Could not remove $d (files still open on NFS?); moved it to $aside" exception = err
            mv(d, aside)
        else
            @warn "Could not remove $d (files still open on NFS?); left in place" exception = err
        end
    end
    return nothing
end

const OUT = engine_dir(cfg.project, KEY)
const RUN_DIR = joinpath(OUT, "run")
const WARMUP_DIR = joinpath(OUT, "warmup")
for d in (RUN_DIR, WARMUP_DIR)
    remove_dir(d; move_aside=true)
end
for f in ("result.jld2", "result_path.txt", "samples.jld2", "benchmark.csv", "size_structure.jld2", "thermal_tolerance.jld2")
    isfile(joinpath(OUT, f)) && rm(joinpath(OUT, f))
end

# =============================================================================
# Domain, factor spec, samples
# =============================================================================

t_domain = @elapsed begin
    dom, dom_fixed, YEARS = load_shared_domain(cfg; require_calibration=(KEY === :coralblox))
    spec = make_spec(KEY, dom_fixed, cfg, YEARS, RUN_DIR)
    Random.seed!(cfg.seed)
    scens = ADRIA.sample(spec, cfg.n_samples)
end
varying = spec.fieldname[.!spec.is_constant]
@info "$LABEL: $(nrow(scens)) scenarios, $(length(varying)) varying factors, years $(first(YEARS))-$(last(YEARS))"

# =============================================================================
# Warm-up (1 scenario, untimed for the comparison) and the timed run
# =============================================================================

first_run = (time_s=NaN, gc_s=NaN, alloc_gb=NaN, compile_s=NaN)
if cfg.warmup
    @info "$LABEL: warm-up run (1 scenario)"
    _, first_run = timed_call(() -> run_engine(KEY, dom, dom_fixed, scens[1:1, :], spec, WARMUP_DIR, cfg, YEARS))
    remove_dir(WARMUP_DIR)   # a leftover is only disk space; the next run of this engine clears it
end

@info "$LABEL: main run ($(nrow(scens)) scenarios)"
rs, main_run = timed_call(() -> run_engine(KEY, dom, dom_fixed, scens, spec, RUN_DIR, cfg, YEARS))

# =============================================================================
# Save: result set, samples + settings, size structure, benchmark
# =============================================================================

save_engine_result(KEY, rs, OUT)
run_settings = settings_dict(cfg, YEARS)
## The calibration file only affects CoralBlox: say so in the other engines' records instead of a path they never used.
KEY === :coralblox || (run_settings["CALIBRATION_PATH"] = "n/a (CoralBlox only)")
CSCAPE_DEMOGRAPHY === nothing || (run_settings["CSCAPE_DEMOGRAPHY"] = CSCAPE_DEMOGRAPHY)
JLD2.save(joinpath(OUT, "samples.jld2"), Dict(
    "scens" => scens, "varying" => String.(varying), "spec" => spec_for_disk(spec),
    "settings" => run_settings))

t_size = @elapsed (ss = extract_size_structure(KEY, rs, RUN_DIR, dom, cfg, scens))
ss === nothing || JLD2.save(joinpath(OUT, "size_structure.jld2"), ss)
t_thermal = @elapsed (tt = extract_thermal_tolerance(KEY, rs, RUN_DIR, dom, cfg, scens))
tt === nothing || JLD2.save(joinpath(OUT, "thermal_tolerance.jld2"), tt)

n = nrow(scens)
n_locs = length(rs.loc_ids)
bench = (
    engine=String(KEY), label=LABEL, host=gethostname(), date=string(now()),
    n_scenarios=n, n_locations=n_locs, n_years=length(YEARS),
    threads=Threads.nthreads(), n_workers=cfg.n_workers, cpus=Sys.CPU_THREADS,
    slurm_job_id=get(ENV, "SLURM_JOB_ID", ""), slurm_cpus=get(ENV, "SLURM_CPUS_PER_TASK", ""),
    t_load_s=T_LOADED - T_PROCESS_START, t_domain_s=t_domain,
    t_first_s=first_run.time_s, t_first_compile_s=first_run.compile_s, t_first_alloc_gb=first_run.alloc_gb,
    t_main_s=main_run.time_s, t_main_compile_s=main_run.compile_s, t_main_gc_s=main_run.gc_s,
    t_main_alloc_gb=main_run.alloc_gb, t_per_scen_s=main_run.time_s / n,
    t_per_scen_site_year_ms=1000 * main_run.time_s / (n * n_locs * length(YEARS)),
    t_size_extract_s=t_size, t_thermal_extract_s=t_thermal, maxrss_gb=Sys.maxrss() / 1e9, disk_gb=dir_size_gb(RUN_DIR),
    disk_per_scen_gb=dir_size_gb(RUN_DIR) / n, t_total_s=time() - T_PROCESS_START,
)
write_csv(joinpath(OUT, "benchmark.csv"), [bench])

@info "$LABEL done" OUT t_main_s = round(main_run.time_s; digits=1) t_per_scen_s = round(bench.t_per_scen_s; digits=2) maxrss_gb = round(bench.maxrss_gb; digits=2) size_structure = ss !== nothing thermal_tolerance = tt !== nothing
