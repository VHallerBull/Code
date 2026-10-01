# Source clones

The Julia packages this environment uses are local development clones (see `Project.toml` `[sources]`).
Recorded by `record_sources.jl` on 2026-10-02 07:52.

| Repository (local path) | Packages used | Branch | Commit | Uncommitted changes |
|---|---|---|---|---|
| `C:\Users\vhaller\Documents\GitHub\ADRIA.jl` | ADRIA, ADRIAanalysis, ADRIAviz | CscapeIntegration | cbd8d30f Add engine=:kora to run_scenarios via a weak-dependency extension | **yes** |
| `C:\Users\vhaller\Documents\GitHub\ADRIAIndicators.jl` | ADRIAIndicators | main | baf4b68 Merge pull request #23 from open-AIMS/bump-version-0.3 | no |
| `C:\Users\vhaller\Documents\GitHub\CORSET.jl` | CORSET | master |  | no |
| `C:\Users\vhaller\Documents\GitHub\CORSETInterface.jl` | CORSETInterface | master |  | no |
| `C:\Users\vhaller\Documents\GitHub\CScape.jl` | CScape | ADRIAIntegration | ae3d1f8 Cache create_ipm's growth-kernel Normal-CDF pass per (wave_exposure, depth_cat) | **yes** |
| `C:\Users\vhaller\Documents\GitHub\CoralPopInterface.jl` | CoralPopInterface | main | 14ffa59 README: Phase 7 lives in CoralModelReview/Code | **yes** |
| `C:\Users\vhaller\Documents\GitHub\CscapeInterface.jl` | CscapeInterface | CScapeUpdate | 8b63533 Fix relative shelter volume: use ADRIA's 95 cm reference colony, clamp to [0, 1] | **yes** |
| `C:\Users\vhaller\Documents\GitHub\Kora.jl` | Kora | fix/real-domain-area-and-interpolations-deps | b627e99 fix: initialize_reef crashed on real (non-round) per-location area/density | no |
| `C:\Users\vhaller\Documents\GitHub\KoraInterface.jl` | KoraInterface | main | 6f97b86 Add KoraInterface.jl: run Kora as a fourth ADRIA engine | no |
| `C:\Users\vhaller\Documents\GitHub\ReefModInterface.jl` | ReefModInterface | master |  | **yes** |

Other local inputs the scripts read (not packages):
- the ADRIA domain package and CoralBlox calibration file named in each script's CONFIG section
- the CoralPop Python engine, `CoralPop-Engine` (branch `engine-api`), used through CoralPopInterface
  - C:\Users\vhaller\Documents\GitHub\CoralPop-Engine: engine-api @ a4ef31a Engine: cap recruits at free substrate; floor benthos and settlement at zero
