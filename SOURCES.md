# Source clones

The Julia packages this environment uses are local development clones (see `Project.toml` `[sources]`).
Recorded by `record_sources.jl` on 2026-10-06 19:40.

| Repository (local path) | Packages used | Branch | Commit | Uncommitted changes |
|---|---|---|---|---|
| `C:\Users\vhaller\Documents\GitHub\ADRIA.jl` | ADRIA, ADRIAanalysis, ADRIAviz | CscapeIntegration | 97fb44ee :coralpop/:kora engines: forward larvae_via_connectivity | no |
| `C:\Users\vhaller\Documents\GitHub\ADRIAIndicators.jl` | ADRIAIndicators | main | 6d22933 Merge pull request #26 from open-AIMS/juve-vec-case | no |
| `C:\Users\vhaller\Documents\GitHub\CORSET.jl` | CORSET | master | 11f19e3 Fix regional coral recruitment double-counting the local supply | no |
| `C:\Users\vhaller\Documents\GitHub\CORSETInterface.jl` | CORSETInterface | master | e31e24e Store the domain's initial cover as timestep 1; store k as max cover | no |
| `C:\Users\vhaller\Documents\GitHub\CScape.jl` | CScape | ADRIAIntegration | 4fc1dc8 Point the demography export at CscapeInterface; drop data-raw's environment | no |
| `C:\Users\vhaller\Documents\GitHub\CoralPopInterface.jl` | CoralPopInterface | main | 21a9632 feat: optional larvae_via_connectivity (matrix self-retention, own retention off) | no |
| `C:\Users\vhaller\Documents\GitHub\CscapeInterface.jl` | CscapeInterface | CScapeUpdate | f403236 Add ADRIAviz and CairoMakie to the package environment | no |
| `C:\Users\vhaller\Documents\GitHub\Kora.jl` | Kora | integration/combined-fixes | 334f281 Merge branch 'feat/larval-connectivity' into integration/combined-fixes | no |
| `C:\Users\vhaller\Documents\GitHub\KoraInterface.jl` | KoraInterface | main | b1bde13 feat: optional larval connectivity (larvae_via_connectivity) | no |
| `C:\Users\vhaller\Documents\GitHub\ReefModInterface.jl` | ReefModInterface | master | fb5b4dd Run ReefMod via matlab subprocesses as ADRIA's :reefmod engine, not Slurm-only | no |

Other local inputs the scripts read (not packages):
- the ADRIA domain package and CoralBlox calibration file named in each script's CONFIG section
- the CoralPop Python engine, `CoralPop-Engine` (branch `engine-api`), used through CoralPopInterface
  - C:\Users\vhaller\Documents\GitHub\CoralPop-Engine: engine-api @ 032824c run_network: include_self keeps the connectivity diagonal as self-recruitment
