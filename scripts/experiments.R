# experiments.R
# experimenting with the model running
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────
library(dplyr)
library(readr)
library(rgee)
source("scripts/helper_functions.R")
source("scripts/get_climateinputs.R")
source("scripts/get_microenv.R")
source("scripts/get_colonization.R")
source("scripts/paths.R")

# run steps in main.R first

# ── Logging ───────────────────────────────────────────────────────────────────
dir.create(LOGS_DIR, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(
  LOGS_DIR,
  sprintf("experiment_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
log_msg <- function(msg) {
  stamped <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)
  message(stamped)
  cat(stamped, "\n", file = log_file, append = TRUE)
}
log_msg("experiments.R started")

# sourcing
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

sapply(c("Maquipucuna", "Mashpi", "MindoTarabita", "MiradorMindo", "Yanayacu"), function(s) {
  site_keys <- names(models)[grepl(s, names(models))]
  sum(sapply(site_keys, function(k) {
    any(sapply(models[[k]], function(cell) inherits(cell, "micropoint")))
  }))
})
prep <- prepare_observations("data/csv/combinedv3.csv", models)
obs <- prep$obs
niches <- extract_niches(obs, models, prep$valid_per_model)
sites_df <- make_sites("data/csv/combinedv3.csv", pad = 0.15)
site <- sites_df[sites_df$Site == "Maquipucuna", ]
resolution <- 10 # metres per cell
site_obs <- niches[niches$Area_or_Site == site$Site, ]
xDim <- max(round((max(site_obs$lon) - min(site_obs$lon)) * 111000 / resolution), 10)
yDim <- max(round((max(site_obs$lat) - min(site_obs$lat)) * 111000 / resolution), 10)
canopy_grid <- get_canopy_grid(site, xDim, yDim,
                               resolution = 5,
                               out_dir = file.path(RAW_DIR, site$Site)
)


# EXPERIMENT

# first I wanna play with beta0Seedling from -1.5 to 1.5 by 0.1 
seedlingSurvival <- seq(from = -1.5, to = 1.5, by = 0.1)

results <- list()

for (i in seedlingSurvival){
  # setup parameters
  params <- list(
    # growth parameters
    beta0GrowthS = 0.22, # seedling growth rate derived from bulb growth per year
    beta0GrowthJ = 0.45, # juvenile growth rate derived from bulb growth per year
    beta0GrowthA = 1.40, # adult growth rate derived from bulb growth per year
    
    beta1Growth = 0.8, # I believe this is mean growth
    sigma = 0.5, # residual variance (assumed independent of size)
    
    # survival parameters
    beta0Seedling = i, # survival rate seedlings
    beta0Juvenile = -0.5, # survival rate juveniles
    beta0Adult = 1.0, #survival rate adults
    
    beta1 = 0.105, # mean survival I guess? 
    
    # reproduction parameters
    p_flower = 0.80, # proportion of adults that flower
    p_poll = 0.5, # pollination probability / fruit set
    p_germ = 0.5, # mycorrhizal germination probability
    p_s1 = 0.9, # first-year seedling survival
    
    # dispersal parameters
    Ut = 1, # Terminal seed velocity
    lambda = 1, # Spread factor
    canopy_z = mean(canopy_grid, na.rm = TRUE) # Mean canopy height
  )
  
  results[[i]] <- runcolonization(
    site = site,
    niches = niches,
    canopy_grid = canopy_grid,
    models = models,
    valid_per_model = prep$valid_per_model,
    timesteps = 50, # short run for debugging
    resolution = 10,
    carCap = 1,
    maxDisp = 5,
    spinup = 10,
    stochastic = FALSE,
    Visualize = TRUE,
    sleeptime = 0.2,
    visualize_dispersion = FALSE,
    parameters = params
  )
}

