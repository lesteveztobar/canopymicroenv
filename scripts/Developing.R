# Developing.R
# Canopy colonization model — vertical niche partitioning of epiphytic Maxillariinae
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────

library(fields)      # image.plot for visualization
library(viridis)     # colour palettes
library(dplyr)
library(readr)
library(microclimdata)
library(microclimf)

source("scripts/functions.R")
source("scripts/helper_functions.R")
source("scripts/paths.R")

# ── Load microclimate outputs ─────────────────────────────────────────────────
models <- readRDS(file.path(BASE_DIR, "data/processed/pointmodelv1.rds"))

# ── Actual development space  ─────────────────────────────────────────────────
# nothing works (:

mout <- runpointmodel(weather = climdata, dtm = dtmcaerth, vegp = vegp, soilc = soilc)
mouttt <- runmicro(micropoint = mout, reqhgt = 3, vegp = vegp, soilc = soilc, dtm = dtmcaerth)
