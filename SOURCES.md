# Source clones

The Julia packages this environment uses are local development clones (see `Project.toml` `[sources]`).
Recorded by `record_sources.jl` on 2026-10-06 15:48.

| Repository (local path) | Packages used | Branch | Commit | Uncommitted changes |
|---|---|---|---|---|
| `C:\Users\vhaller\Documents\GitHub\ADRIA.jl` | ADRIA, ADRIAanalysis, ADRIAviz | CscapeIntegration | 28835f5f Add engine=:reefmod to run_scenarios via a weak-dependency extension | no |
| `C:\Users\vhaller\Documents\GitHub\ADRIAIndicators.jl` | ADRIAIndicators | main | 6d22933 Merge pull request #26 from open-AIMS/juve-vec-case | no |
| `C:\Users\vhaller\Documents\GitHub\CORSET.jl` | CORSET | master | 11f19e3 Fix regional coral recruitment double-counting the local supply | no |
| `C:\Users\vhaller\Documents\GitHub\CORSETInterface.jl` | CORSETInterface | master | e31e24e Store the domain's initial cover as timestep 1; store k as max cover | no |
| `C:\Users\vhaller\Documents\GitHub\CScape.jl` | CScape | ADRIAIntegration | 4fc1dc8 Point the demography export at CscapeInterface; drop data-raw's environment | no |
| `C:\Users\vhaller\Documents\GitHub\CoralPopInterface.jl` | CoralPopInterface | main | 64f4d15 Sandbox: default domain-translation check to ADRIA's bundled Test_domain | no |
| `C:\Users\vhaller\Documents\GitHub\CscapeInterface.jl` | CscapeInterface | CScapeUpdate | f403236 Add ADRIAviz and CairoMakie to the package environment | no |
| `C:\Users\vhaller\Documents\GitHub\Kora.jl` | Kora | integration/combined-fixes | 660be58 Merge branch 'fix/test-load-packaged-models' into integration/combined-fixes | no |
| `C:\Users\vhaller\Documents\GitHub\KoraInterface.jl` | KoraInterface | main | 6d7feb8 feat: simulate one 250 m² section per site and extrapolate to the site | no |
| `C:\Users\vhaller\Documents\GitHub\ReefModInterface.jl` | ReefModInterface | master | fb5b4dd Run ReefMod via matlab subprocesses as ADRIA's :reefmod engine, not Slurm-only | no |

Other local inputs the scripts read (not packages):
- the ADRIA domain package and CoralBlox calibration file named in each script's CONFIG section
- the CoralPop Python engine, `CoralPop-Engine` (branch `engine-api`), used through CoralPopInterface
  - C:\Users\vhaller\Documents\GitHub\CoralPop-Engine: engine-api @ a4ef31a Engine: cap recruits at free substrate; floor benthos and settlement at zero
