# canopymicroenv

> **Microclimate niche modelling and vertical colonization simulation for epiphytic orchids — NW Ecuador**
> Master's thesis project — University of Bonn, 2026

---

## What this is

This repository contains the full analysis pipeline for my master's thesis on the **vertical stratification of epiphytic Maxillariinae orchids** and their microenvironmental correlates across five cloud forest sites in the Chocó Andino of northwestern Ecuador (Maquipucuna, MiradorMindo, Mashpi, MindoTarabita, Yanayacu).

The core idea: instead of using coarse climate data to describe orchid habitat, this pipeline models the **exact microclimate at the height and location where each individual was observed** — at 0.1 m resolution across the full canopy vertical gradient. These microclimate profiles characterise the realised niche of each species and feed a **3D spatially explicit colonization model** that simulates population dynamics (dispersal, establishment, survival, growth) across the canopy landscape, driven by microclimate-based vital rates.

> ⚠️ **Work in progress.** Microclimate pipeline complete for all 5 sites; colonization model and sensitivity experiments running. See [Status](#status).

---

## Study system

- **Focal group:** Maxillariinae orchids (tribe Maxillarieae, subtribe Maxillariinae)
- **Sites:** 5 cloud forest sites, ~800–1800 m elevation, NW Ecuador
- **Observations:** 141 individuals across all sites, with observed canopy height, elevation, genus, and photo references
- **Approach:** Field observations → microclimate modelling → niche characterisation → 3D vertical colonization simulation → parameter sensitivity experiments

---

## Pipeline overview

```
GeoJSON field exports
        ↓
  scripts/convert_observations.py     # parse field notes → structured CSV
        ↓
  geojson_to_csv/csv/                 # per-site CSVs (one per field site)
        ↓
  data/csv/combinedv3.csv             # manually edited, merged dataset with final species IDs
        ↓
  run_microclimate.R                  # interactive: all sites sequentially
  run_microclimate_site.R             # HPC: single site via Rscript / SLURM
    └── get_climateinputs.R           # ERA5 download + preprocessing
    └── get_microenv.R                # vegetation · soil · runpointmodela()
                                      # height loop: 0.1 m steps, per-height RDS to scratch
  microenv_<site>.rds                 # manifest → data/processed/ (~350 KB)
  /lustre/scratch/.../microenv_<site>_heights/  # per-height spatial arrays (~6 GB each)
        ↓
  run_colonization_onesite.R          # single-site colonization run
  run_colonization_allsites.R         # all-sites loop (or SLURM array)
    └── get_colonization.R            # build_forest() · dispersal · establishment
                                      # survival/growth (IPM-style, equation primitives)
        ↓
  make_params.R                       # build parameter sweep RDS files for experiments
  batch_exp.sh / run_colonization.sh  # SLURM: 5 sites × 6 experiments in parallel
  run_experiments.R                   # parallel sensitivity experiments (thesis model)
  simple_colonization.R               # standalone 3D model without microclimate
  simple_experiments.R                # sensitivity experiments for simple model
  simple_model/run_extinction_heatmap.R  # fine-scale extinction threshold scan
  plot_all.R / plot_functions.R       # all project figures in one pass
```

---

## Colonization model

The colonization model simulates a 3D canopy landscape as a voxel grid `[x, y, z]` where z is height in the canopy. Each run goes through:

1. **Forest structure** (`build_forest`): stochastic tree placement from forest inventory parameters (stem density, crown radius, height distribution). Each tree's trunk and crown are labelled with Johansson (1974) zones (JZ1–JZ5), and the maximum epiphyte carrying capacity of each voxel (`carCap_voxel`) is derived from available bark surface area.

2. **Dispersal** (pass 1): adults produce seeds; each seed travels an exponentially distributed distance in the downwind direction. Dispersal is vectorized using `tabulate()` over a padded array.

3. **Establishment** (pass 2): seeds become seedlings if the voxel is within the canopy, below carrying capacity, and a germination event occurs. Germination probability is scaled by local relative humidity and shortwave radiation from `get_clim()`.

4. **Survival and growth** (pass 3): three stage classes — seedlings (S, 0–1 cm pseudobulb), juveniles (J, 1–7 cm), adults (A, 7–20 cm). Survival is computed via `survival_logit()` — logistic in size and microclimate (temperature, RH, light). Stage transitions use `transition_logit()` driven by annual precipitation and RH. Pseudobulb growth uses `size_increment()` with a cost-of-reproduction penalty for fruiting individuals. Each function names its literature source directly in the code.

All inner loops run at height-slice level (not voxel level), replacing `rbinom(1,...)` per voxel with `rbinom(n,...)` per height slice. Experiments are parallelized across parameter values via `parallel::mclapply`.

---

## HPC usage (Marvin cluster, University of Bonn)

The microclimate model is computationally intensive. Each site is submitted as an independent SLURM job via `run_microenv.sh`, coordinated by `microenv_array.sh`:

```bash
# From /home/s38leste_hpc/canopymicroenv/
sbatch scripts/microenv_array.sh 12          # all 5 sites × 12 months of ERA5
sbatch scripts/run_microenv.sh Maquipucuna 12  # single site

# Sensitivity experiments (after microclimate is done)
sbatch scripts/batch_exp.sh                  # 5 sites × 6 experiments = 30 jobs

# Test environment before a full run
sbatch scripts/hpc_test.sh
```

Each microclimate job runs on the `lm_short` partition (large-memory nodes) with 4 CPUs and 500 GB RAM. It loads R/4.4.2, the Miniforge3 conda environment (`canopy_rgee`) for Earth Engine access, and allocates a Lustre scratch workspace for per-height temp files via `ws_allocate`. Logs are written to `logs/log_<jobid>.out`.

**Lustre scratch workspace:** per-height `.rds` files (~6 GB each) are written to `/lustre/scratch/data/s38leste_hpc-canopymicroenv/microenv_<site>_heights/` during computation. The final `microenv_<site>.rds` in `data/processed/` is a small manifest (~350 KB) that records the height vector, scratch path, and baseline weather — not the full spatial arrays. **Do not release the scratch workspace** (`ws_release`) until downstream analysis is complete. The workspace expires in 90 days and can be extended up to 3 times.

---

## Repository structure

```
canopymicroenv/
├── scripts/
│   │   # ── Microclimate pipeline ──────────────────────────────────────────
│   ├── run_microclimate.R             # data acquisition + point model + grid model (interactive)
│   ├── run_microclimate_site.R        # same, single site — called by SLURM
│   │
│   │   # ── Colonization model ────────────────────────────────────────────
│   ├── allsites.R                     # all-sites run: load models, extract niches, run colonization
│   ├── run_colonization_onesite.R     # colonization run: single site (interactive / SLURM)
│   ├── run_colonization_allsites.R    # colonization loop: all sites
│   ├── run_experiments.R              # sensitivity experiment driver (parallel, thesis model)
│   ├── make_params.R                  # build parameter sweep RDS files for batch_exp.sh
│   │
│   │   # ── Simple / standalone model ─────────────────────────────────────
│   ├── simple_colonization.R          # standalone 3D colonization model (no microclimate)
│   ├── simple_experiments.R           # simple model sensitivity experiments + animations
│   ├── simple_model/                  # additional simple-model analyses
│   │   ├── run_colonization_allsites.R
│   │   ├── run_experiments.R
│   │   └── run_extinction_heatmap.R   # extinction threshold scan (p_est × repro_rate)
│   │
│   │   # ── Plotting ──────────────────────────────────────────────────────
│   ├── plot_all.R                     # entry point: runs all project figures
│   ├── plot_functions.R               # shared plotting helpers
│   │
│   │   # ── Function libraries ────────────────────────────────────────────
│   ├── get_colonization.R             # colonization functions: build_forest(),
│   │                                  #   runcolonization(), pass1/2/3 sub-models,
│   │                                  #   survival_logit(), transition_logit(), size_increment(),
│   │                                  #   load_height() (manifest-aware height loader)
│   ├── get_microenv.R                 # microenvironment functions: niche extraction,
│   │                                  #   canopy grid, climate helpers, get_clim()
│   ├── get_climateinputs.R            # ERA5 download + preprocessing, concat_era5_nc()
│   ├── microenv_helpers.R             # additional microenvironment utilities
│   ├── helper_functions.R             # logging, coordinate parsing, DMS normalisation
│   ├── patches.R                      # runtime monkey-patches for ecmwfr / microclimdata
│   ├── paths.R                        # directory constants (BASE_DIR, RAW_DIR, PARAMS_DIR, …)
│   │
│   │   # ── SLURM scripts ─────────────────────────────────────────────────
│   ├── microenv_array.sh              # submit one microclimate job per site
│   ├── run_microenv.sh                # single-site microclimate job (lm_short, 4 CPUs, 500G)
│   ├── batch_exp.sh                   # submit 5 sites × 6 experiments = 30 colonization jobs
│   ├── run_colonization.sh            # single site × experiment colonization job (lm_short, 500G)
│   ├── run_plots.sh                   # run plot_all.R (local or interactive node)
│   ├── hpc_test.sh                    # connectivity / environment smoke test
│   ├── test_env.R                     # R package + Python environment check
│   │
│   │   # ── Data conversion (GeoJSON → CSV) ──────────────────────────────
│   ├── convert_observations.py        # parse iNaturalist / field GeoJSON → structured CSV
│   ├── build_photo_lookup.py          # build photo-lookup table from field exports
│   └── geojson-csv-sql-conversion-tools/  # Node.js + Python toolkit for
│                                          #   GeoJSON ↔ CSV ↔ SQL round-trips;
│                                          #   used for QA and manual corrections
│                                          #   during the observation cleaning step.
│                                          #   © Rudo Kemper, GPL-3.0
│                                          #   github.com/rudokemper/geojson-csv-sql-conversion-tools
├── geojson_to_csv/
│   ├── raw/                           # original GeoJSON exports from iNaturalist / OrganicMaps app
│   └── csv/                           # per-site CSVs produced by convert_observations.py
├── data/
│   ├── csv/                           # combined*.csv, Processed*.csv — observation datasets
│   ├── raw/                           # ERA5, DTM, LAI, albedo, etc. (not tracked)
│   ├── params/                        # parameter sweep RDS files built by make_params.R
│   └── processed/                     # microenv_*.rds (manifests, ~350 KB each),
│                                      #   pointmodel_*.rds, colonization outputs (not tracked)
│                                      #   NB: full per-height arrays (~6 GB each) live in Lustre scratch
├── output/                            # figures and animation outputs (not tracked)
└── logs/                              # timestamped run logs (not tracked)
```

---

## Status

| Step | Status |
|---|---|
| Field data collection (5 sites) | ✅ Complete |
| GeoJSON → CSV conversion | ✅ Complete |
| Observation cleaning + `hObs` extraction | ✅ Complete |
| ERA5 climate data download | ✅ Complete |
| DTM, landcover, vegetation, soil parameters | ✅ Complete |
| Point model height loop (`runpointmodela`, all 5 sites) | ✅ Complete |
| Grid microclimate model (`runmicro`, all 5 sites) | ✅ Complete (height files in Lustre scratch) |
| Colonization model — single site | 🔄 In progress |
| Colonization model — all sites | 🔄 In progress |
| Forest structure (`build_forest`, Myster 2017 params) | ✅ Implemented |
| IPM equation primitives (`survival_logit` etc.) | ✅ Implemented |
| Parameter sensitivity experiments | ✅ Implemented |
| Extinction threshold heatmap | ✅ Complete |
| Species identification | 🔄 In progress |

---

## Key dependencies

- [`microclimf`](https://github.com/ilyamaclean/microclimf) — mechanistic microclimate model (Maclean 2026)
- [`microclimdata`](https://github.com/ilyamaclean/microclimdata) — automated input data acquisition
- [`mcera5`](https://github.com/dklinges9/mcera5) — ERA5 climate data download
- [`rgee`](https://github.com/r-spatial/rgee) — Google Earth Engine interface from R
- [`plotly`](https://plotly.com/r/) — interactive 3D visualisation
- [`gganimate`](https://gganimate.com/), [`gifski`](https://gif.ski/) — dispersal animations
- [`patchwork`](https://patchwork.data-imaginist.com/) — multi-panel plots
- [`sf`](https://r-spatial.github.io/sf/), `ggplot2`, `ggspatial`, `rnaturalearth` — spatial visualisation
- `terra`, `parallel`, `dplyr`, `readr`

---

## Notes

- **Johansson (1974) zones** (JZ1–JZ5) are used to assign habitat suitability and bark surface area within `build_forest()`. Zone boundaries are proportional canopy height: JZ1 < 10 %, JZ2 < 30 %, JZ3 < 50 %, JZ4 < 80 %, JZ5 = emergent crown. Carrying capacity per voxel is derived from trunk or effective crown surface area divided by mean epiphyte footprint.
- **Forest inventory baseline:** Myster (2017), Maquipucuna primary cloud forest, 1400 m: mean dsh 22.7 cm (trunk radius 0.114 m), 272–324 stems/ha (≥10 cm dsh). Used as default `forestparams` in one-site and all-sites runs.
- **Incremental microclimate saves:** `run_microclimate_site.R` saves each height to its own `.rds` in Lustre scratch before proceeding. If the job is killed (e.g. SLURM timeout or OOM), resubmitting resumes from the last completed height automatically.
- **microenv manifest format:** `microenv_<site>.rds` is a small list with `.heights` (numeric vector), `.height_dir` (path to scratch), and `.weather` (ERA5-cell baseline). The full 6 GB per-height spatial arrays stay in scratch. Load a specific height with `readRDS(file.path(microenv$.height_dir, sprintf("h%.2f.rds", h)))`.
- **HPC path:** `paths.R` sets `BASE_DIR` to the HPC home directory. The `CANOPY_PYTHON` environment variable controls which Python interpreter is used for Earth Engine; it defaults to the `canopy_rgee` conda environment.
- **Runtime patches:** two monkey-patches are applied in `patches.R` to fix known bugs in `ecmwfr` and `microclimdata` without modifying package source.
- **Credentials:** `credentials.rds` (CDS API, NASA Earthdata, Google credentials) is excluded from version control. You will need your own.
- **Array dimensions:** the colonization model tracks `[x, y, z, timestep, species]` per life stage. Heights are scoped to the current site's observed range to keep the z-dimension manageable.

---

*Lizeth Estévez Tobar — Universität Bonn, 2026*  
*Supervisor: Juliano Sarmento Cabral*
