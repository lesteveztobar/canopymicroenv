# main.R
# Canopy colonization model — vertical niche partitioning of epiphytic Maxillariinae
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
library(dplyr)
library(readr)
library(rgee)
library(rgl)
library(plotly)
source("scripts/helper_functions.R")
source("scripts/get_climateinputs.R")
source("scripts/get_microenv.R")
source("scripts/get_colonization.R")
source("scripts/paths.R")

# ── Setup virtual environment ─────────────────────────────────────────────────
# reticulate::use_virtualenv("/Users/lizethestevezt/.virtualenvs/rgee", required = TRUE)

# ── Logging ───────────────────────────────────────────────────────────────────
dir.create(LOGS_DIR, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(
  LOGS_DIR,
  sprintf("main_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
log_msg <- function(msg) {
  stamped <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)
  message(stamped)
  cat(stamped, "\n", file = log_file, append = TRUE)
}
log_msg("main.R started")
# ── STEP 1: Load microclimate outputs from getmicroenv.R ─────────────────────
# Models are saved per site — load all and merge into one list
model_files <- list.files(PROCESSED_DIR, pattern = "^pointmodel_.*\\.rds$", full.names = TRUE)

models <- list()
for (f in model_files) {
  site_models <- readRDS(f)
  models <- c(models, site_models)
}

log_msg(sprintf(
  "Loaded %d height models across %d sites",
  length(models), length(model_files)
))

# Check — any valid cell per site
sapply(c("Maquipucuna", "Mashpi", "MindoTarabita", "MiradorMindo", "Yanayacu"), function(s) {
  site_keys <- names(models)[grepl(s, names(models))]
  sum(sapply(site_keys, function(k) {
    any(sapply(models[[k]], function(cell) inherits(cell, "micropoint")))
  }))
})

# ── STEP 2: Extract microclimate niche ─────────────
prep <- prepare_observations("data/csv/combinedv3.csv", models)
obs <- prep$obs
niches <- extract_niches(obs, models, prep$valid_per_model)
saveRDS(niches, file.path(PROCESSED_DIR, "niches.rds"))
# ── STEP 3: Run colonization model ───────────────────────────────────────────
sites_df <- make_sites("data/csv/combinedv3.csv", pad = 0.15)
site <- sites_df[sites_df$Site == "Maquipucuna", ]
res <- 10  # set once, use everywhere
# need xDim/yDim first — compute them the same way runcolonization() does
site_obs    <- niches[niches$Area_or_Site == site$Site, ]
lat_range_m <- (max(site_obs$lat) - min(site_obs$lat)) * 111000
lon_range_m <- (max(site_obs$lon) - min(site_obs$lon)) * 111000 *
  cos(mean(site_obs$lat) * pi / 180)
xDim <- max(round(lon_range_m / res), 10) + 3
yDim <- max(round(lat_range_m / res), 10) + 3

# this uses earth engine so authenticate and initialize
ee_Authenticate()
ee_Initialize(project = "ee-lizethestevezt")
canopy_grid <- get_canopy_grid(site, xDim, yDim,
  resolution = res,
  out_dir = file.path(RAW_DIR, site$Site)
)

# setup parameters
params <- list(
  # ── Growth parameters (Mondragón et al. 2007 — pseudobulb size increments) ──
  # size thresholds for stage transitions — not growth rates per se
  beta0GrowthS = 0.22,  # mean size of seedlings (pseudobulbs)
  beta0GrowthJ = 0.45,  # mean size of juveniles (pseudobulbs)
  beta0GrowthA = 1.40,  # mean size of adults (pseudobulbs)
  
  beta1Growth  = 0.8,   # size autocorrelation — how much current size predicts next size
  # prototype: 0.3 (slower transitions); literature: 0.8
  sigma        = 0.5,   # residual growth variance (stochastic noise around mean growth)
  # prototype: 0.2; literature: 0.5
  
  # ── Survival parameters (Raventós et al. 2015; Mondragón et al. 2007) ──
  # logistic regression intercepts — higher = better survival
  # s = 1 / (1 + exp(-(beta0_stage + beta1 * size - climate_penalty)))
  beta0Seedling = -1.5,  # literature value; prototype forgiving: 1.5
  beta0Juvenile = -0.5,  # literature value; prototype forgiving: 1.0
  beta0Adult    =  1.0,  # literature value; prototype forgiving: 2.5
  # adults have high stasis (Zotz 1998)
  
  beta1 = 0.105,  # size effect on survival — larger individuals survive better
  
  # ── Reproduction parameters ──
  p_flower = 0.80,   # proportion of adults flowering per year (Winkler et al. 2009)
  p_poll   = 0.05,   # pollination success / fruit set probability
  # prototype relaxed: 0.50; literature: 0.05
  p_germ   = 0.001,  # mycorrhizal germination probability (Taylor & Bruns 1999)
  # prototype relaxed: 0.50; literature: 0.001
  p_s1     = 0.83,   # first-year seedling survival (Zotz 1998)
  # prototype relaxed: 0.90
  
  # ── Dispersal parameters (Murren & Ellison 1998) ──
  Ut       = 1,   # terminal seed fall velocity (m/s) — affects mean dispersal distance
  lambda   = 1,   # kernel spread factor — higher = wider dispersal
  canopy_z = mean(canopy_grid, na.rm = TRUE)  # mean canopy height across site (m)
)

forgivingparams <- list(
  # ── Growth parameters ──
  beta0GrowthS = 0.22,  # mean size of seedlings (pseudobulbs)
  beta0GrowthJ = 0.45,  # mean size of juveniles (pseudobulbs)
  beta0GrowthA = 1.40,  # mean size of adults (pseudobulbs)
  
  beta1Growth  = 0.4,   # size autocorrelation — literature: 0.8; forgiving: 0.5
  sigma        = 0.2,   # residual growth variance — literature: 0.5; forgiving: 0.3
  
  # ── Survival parameters ──
  beta0Seedling =  0.8,  # literature: -1.5; forgiving: 0.5
  # gives ~57% monthly survival at relhum=85%
  beta0Juvenile =  1,  # literature: -0.5; forgiving: 0.8
  # gives ~68% monthly survival
  beta0Adult    =  1.5,  # literature: 1.0; forgiving: 1.5
  # gives ~82% monthly survival; stasis dominates (Zotz 1998)
  
  beta1 = 0.105,  # size effect on survival — larger individuals survive better
  
  # ── Reproduction parameters ──
  p_flower = 0.80,   # proportion of adults flowering (Winkler et al. 2009)
  p_poll   = 0.20,   # pollination success — literature: 0.05; forgiving: 0.20
  p_germ   = 0.1,   # mycorrhizal germination — literature: 0.001; forgiving: 0.05
  p_s1     = 0.83,   # first-year seedling survival (Zotz 1998)
  # repRate = 0.80 * 0.20 * 0.05 * 0.83 = 0.00664
  # reproduce(6, 0.00664) → ~4% chance of 1 seed per adult per year
  # with forgiving survival adults accumulate → more seeds over time
  
  # ── Dispersal parameters (Murren & Ellison 1998) ──
  Ut       = 1,   # terminal seed fall velocity (m/s)
  lambda   = 1,   # kernel spread factor
  canopy_z = mean(canopy_grid, na.rm = TRUE)
)

superforgivingparams <- list(
  # ── Growth parameters ──
  beta0GrowthS = 0.22,  # mean size of seedlings (pseudobulbs)
  beta0GrowthJ = 0.45,  # mean size of juveniles (pseudobulbs)
  beta0GrowthA = 1.40,  # mean size of adults (pseudobulbs)
  
  beta1Growth  = 0.3,   # slow transitions — literature: 0.8
  sigma        = 0.8,   # very low variance — predictable slow growth
  
  # ── Survival parameters ──
  # at relhum=85%, dry_penalty=0.15; at temp=17°C, temp_penalty=0
  beta0Seedling =  2.0,  # ~88% monthly survival → ~21% annual
  beta0Juvenile =  2.5,  # ~92% monthly survival → ~37% annual
  beta0Adult    =  3.0,  # ~95% monthly survival → ~54% annual
  
  beta1 = 0.105,  # size effect on survival
  
  # ── Reproduction parameters ──
  # repRate = 0.80 * 0.50 * 0.30 * 0.83 = 0.0996 ≈ 0.10
  # reproduce(7, 0.10) → guaranteed 1 seed most calls
  p_flower = 0.80,  # proportion of adults flowering
  p_poll   = 0.50,  # pollination success
  p_germ   = 0.30,  # mycorrhizal germination
  p_s1     = 0.83,  # first-year seedling survival (Zotz 1998)
  
  # ── Dispersal parameters ──
  Ut       = 1,
  lambda   = 1,
  canopy_z = mean(canopy_grid, na.rm = TRUE)
)

# run log 
dir.create(LOGS_DIR, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(
  LOGS_DIR,
  sprintf("runcolonization_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
log_msg <- function(msg) {
  stamped <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)
  message(stamped)
  cat(stamped, "\n", file = log_file, append = TRUE)
}
log_msg("colonization run started")

# for viewing

result <- runcolonization(
  site = site,
  niches = niches,
  canopy_grid = canopy_grid,
  models = models,
  valid_per_model = prep$valid_per_model,
  timesteps = 40, # short run for debugging
  resolution = res,
  carCap = 1,
  maxDisp = 5,
  spinup = 15,
  stochastic = FALSE,
  Visualize = TRUE,
  sleeptime = 0.2,
  visualize_dispersion = FALSE,
  parameters = superforgivingparams
)

rgl::rglwidget()
plot_3d_abundance(
  result
)

plot_total_abundance(result, t = 40)

# ── STEP 4: Analysis + figures ───────────────────────────────────────────────
