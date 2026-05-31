# canopymicroenv

> **Microclimate niche modelling and vertical colonization simulation for epiphytic orchids — NW Ecuador**
> Master's thesis project — University of Bonn, 2026

---

## What this is

This repository contains the full analysis pipeline for my master's thesis on the **vertical stratification of epiphytic Maxillariinae orchids** and their microenvironmental correlates across five cloud forest sites in the Chocó Andino of northwestern Ecuador (Maquipucuna, MiradorMindo, Mashpi, MindoTarabita, Yanayacu).

The core idea: instead of using coarse climate data to describe orchid habitat, this pipeline models the **exact microclimate at the height and location where each individual was observed** — at 0.1 m resolution across the full canopy vertical gradient. These microclimate profiles characterise the realised niche of each species and feed a **3D spatially explicit colonization model** that simulates population dynamics (dispersal, establishment, survival, growth) across the canopy landscape, driven by microclimate-based vital rates.

> ⚠️ **Work in progress.** The pipeline runs end-to-end but statistical analysis is not yet started. See [Status](#status).

---

## Study system

- **Focal group:** Maxillariinae orchids (tribe Maxillarieae, subtribe Maxillariinae)
- **Sites:** 5 cloud forest sites, ~800–1800 m elevation, NW Ecuador
- **Observations:** 141 individuals across all sites, with observed canopy height, elevation, genus, and photo references
- **Approach:** Field observations → microclimate modelling → niche characterisation → vertical colonization simulation → statistical analysis

---

## Pipeline overview

```
GeoJSON field exports
        ↓
  convert_observations.py       # parse field notes → structured CSV
        ↓
  combined.csv                  # merged observation dataset (all sites)
        ↓
  combinedv3.csv                # manually edited dataset to add final ID for each observation
        ↓
  runmicroenv.R                 # ERA5 climate · DTM · landcover
    └── get_microenv.R          # vegetation · soil · runpointmodela()
        ↓                       # height loop: 0.1 m steps per site
  pointmodel_<site>.rds         # micropoint objects saved to data/processed/
        ↓
  main.R  ── or ──  allsites.R  # load models → extract niches → run colonization
    └── get_colonization.R      # dispersal · establishment · survival/growth
                                # plot_abundance() · plot_3d_abundance()
        ↓
  experiments.R                 # parameter sensitivity runs
  plotmap.R                     # field site map from GeoJSON
```

---

## Repository structure

```
canopymicroenv/
├── scripts/
│   ├── runmicroenv.R              # executable: data acquisition + point model loop
│   ├── main.R                     # executable: per-site colonization run
│   ├── allsites.R                 # executable: all-sites colonization loop
│   ├── experiments.R              # parameter sensitivity experiments
│   ├── plotmap.R                  # field site map (sf + ggplot2)
│   ├── get_microenv.R             # microenvironment functions — observation
│   │                              #   prep, niche extraction, canopy grid,
│   │                              #   climate helpers, data acquisition wrappers
│   ├── get_colonization.R         # colonization sub-models (dispersal,
│   │                              #   establishment, survival, growth) + plotting
│   ├── get_climateinputs.R        # ERA5 / climate input helpers
│   ├── helper_functions.R         # utilities: logging, coordinate parsing,
│   │                              #   DMS normalisation, network simulation
│   ├── paths.R                    # directory constants (BASE_DIR, RAW_DIR, …)
│   └── config_processing/
│       ├── convert_observations.py  # GeoJSON → structured CSV (Python)
│       └── csv_processing.R         # observation cleaning + hObs extraction
├── geojson_to_csv/
│   ├── raw/                       # raw GeoJSON exports from field app
│   └── csv/                       # per-site processed CSVs
├── data/
│   ├── csv/                       # combined observation dataset (combinedv3.csv)
│   ├── raw/                       # ERA5, DTM, LAI, albedo, etc. (not tracked)
│   └── processed/                 # model outputs — pointmodel_*.rds, niches.rds
│                                  #   colonization_*.rds (not tracked)
├── logs/                          # timestamped run logs (not tracked)
└── output/                        # figures and results (not tracked)
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
| Niche extraction (microclimate per observation) | ✅ Complete |
| Colonization model — per-site (`main.R`) | ✅ Running |
| Colonization model — all sites (`allsites.R`) | 🔄 In progress |
| Parameter sensitivity experiments (`experiments.R`) | 🔄 In progress |
| Species identification | 🔄 In progress |
| Statistical analysis | ⏳ Not yet started |

---

## Key dependencies

- [`microclimf`](https://github.com/ilyamaclean/microclimf) — mechanistic microclimate point model (Maclean 2026)
- [`microclimdata`](https://github.com/ilyamaclean/microclimdata) — automated input data acquisition
- [`mcera5`](https://github.com/dklinges9/mcera5) — ERA5 climate data download
- [`rgee`](https://github.com/r-spatial/rgee) — Google Earth Engine interface from R
- [`scico`](https://github.com/thomasp85/scico) — perceptually uniform colour palettes
- [`plotly`](https://plotly.com/r/) — interactive 3D visualisation
- [`sf`](https://r-spatial.github.io/sf/), `ggplot2`, `ggspatial`, `rnaturalearth` — spatial visualisation
- `terra`, `dplyr`, `readr`

---

## Notes

- Johansson canopy zones (JZ1–JZ5) were considered and dropped in favour of directly observed individual heights (`hObs`). All analysis uses actual measurement data.
- Two monkey-patches are applied at runtime in `runmicroenv.R` to fix known bugs in `ecmwfr` and `microclimdata` without modifying package source.
- `credentials.rds` (CDS API, NASA Earthdata, Google credentials) is excluded from version control. You will need your own.
- The colonization model uses a 5D array `[x, y, z, timestep, species]` per life stage (S/J/A). Heights are scoped to the current site to avoid inflating the z-dimension across sites. Use `allsites = TRUE` in `runcolonization()` to run with the full multi-site height range.

---

*Lizeth Estévez Tobar — Universidad de Bonn, 2026*
*Supervisor: Juliano Sarmento Cabral*
