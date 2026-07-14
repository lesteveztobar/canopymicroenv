# archive.R
# Dead code retired from get_colonization.R, kept here for reference only —
# NOT sourced by any part of the active pipeline. Each function below was
# superseded by a vectorized, per-height-slice reimplementation inline in
# run_pass1_disperse() / run_pass2_establish() / run_pass3_survive_grow()
# (see get_colonization.R), which replace these scalar-per-voxel /
# scalar-per-individual versions with a single rbinom()/vectorized call over
# an entire [xDim, yDim] slice. Confirmed unreferenced anywhere else in the
# codebase before archiving (grep audit, 2026-07-09).
#
# survive(), grow(), and size_increment() are superseded by the inlined
# logic in run_pass3_survive_grow(); establish() is superseded by the
# inlined logic in run_pass2_establish() (which also no longer applies
# p_germ — see get_colonization.R's Pass 2 for the corrected, single-
# application-of-p_germ design). ricker()/stochRicker() were never wired
# into the main loop at all (self-documented as "kept for future use").
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────

# Annual pseudobulb size increment for adults, climate-scaled and with
# optional cost-of-reproduction penalty.
# Δz = delta_z_base × (P/P_ref) × (RH/RH_ref) + ε  [ε ~ N(0, σ/2)]
# If fruited: Δz × cost_repro (Zotz 1998 — "reduced growth after fruiting")
size_increment <- function(z_A, precip_annual, relhum_mean,
                           delta_z_base, precip_ref, rh_ref, sigma,
                           cost_repro, fruited, z_A_min, z_A_max) {
  delta_z <- delta_z_base *
    (precip_annual / precip_ref) *
    (relhum_mean   / rh_ref) +
    rnorm(1, 0, sigma * 0.5)
  if (fruited) delta_z <- delta_z * cost_repro
  return(pmin(z_A_max, pmax(z_A_min, z_A + delta_z)))
}

# ── s(z, e): Survival probability ─────────────────────────────────────────────
# Stage-structured logistic survival with stage-specific intercepts and climate
# covariates. Size z within each stage is drawn from Uniform(stage_min, stage_max)
# as an approximation (no individual size data; disclosed in methods).
# Form: logit(s) = β₀_stage + β₁·z + climate_terms
# β₀ intercepts calibrated to match annual survival rates per stage:
#   S: ~0.44 (Mondragón et al. 2007, Zotz 1998)
#   J: ~0.60 (Winkler et al. 2009 — mean across species)
#   A: ~0.85 (Zotz & Schmidt 2006, Borrero et al. 2023)
# Climate effects:
#   Seedlings: RH-dependent (Zotz 1998 — "most deaths during dry season")
#   Adults:    temperature-dependent — heat stress above 23°C (Olaya-Arenas 2011)
# Sources: Zotz (1998), Olaya-Arenas et al. (2011), Izuddin et al. (2018),
#          Mondragón et al. (2007), Winkler et al. (2009), Zotz & Schmidt (2006)
survive <- function(N, stage, clim_month, mean_swdown_site,
                    beta0S = -0.24,   # logit intercept for S → annual s ≈ 0.44
                    beta0J =  0.41,   # logit intercept for J → annual s ≈ 0.60
                    beta0A =  1.73,   # logit intercept for A → annual s ≈ 0.85
                    beta1  =  0.10,   # size effect (positive; larger = higher survival)
                    z_S_min = 0,  z_S_max = 1,   # cm — seedling size range
                    z_J_min = 1,  z_J_max = 7,   # cm — juvenile size range
                    z_A_min = 7,  z_A_max = 20,  # cm — adult size range (min ~7cm: Zotz 2006)
                    stochastic = FALSE) {
  if (is.na(N) || N == 0) return(0L)
  temp        <- mean(clim_month$temp,   na.rm = TRUE)
  relhum      <- mean(clim_month$relhum, na.rm = TRUE)
  swdown_mean <- mean(clim_month$swdown[clim_month$swdown > 0], na.rm = TRUE)
  swdown_rel  <- if (!is.na(swdown_mean) && mean_swdown_site > 0)
                   swdown_mean / mean_swdown_site else 1.0
  z <- switch(stage,
    S = runif(1, z_S_min, z_S_max),
    J = runif(1, z_J_min, z_J_max),
    A = runif(1, z_A_min, z_A_max))
  s <- survival_logit(stage, z, temp, relhum, swdown_rel, beta0S, beta0J, beta0A, beta1)
  return(if (stochastic) rbinom(1, N, s) else round(N * s))
}

# ── g(z'|z, e) applied: Advance individuals and update adult size ──────────────
# Wraps growth_prob() for stage transitions (S→J, J→A).
# Also computes the annual size increment for adults (Δz) and applies the
# cost-of-reproduction penalty directly to Δz for cells that fruited.
#
# Annual pseudobulb growth:
#   Δz = delta_z_base × (precip / precip_ref) × (rh / rh_ref) + ε
#   where delta_z_base ≈ 0.8 cm/yr (mean from Zotz 1998 Fig. 3, adults).
#   Fruiting individuals: Δz ← Δz × cost_repro (= ×0.5), so ~0.4 cm/yr.
#   "After reproduction, plants showed reduced vegetative growth" — Zotz (1998).
# Sources: Zotz (1998), Zotz & Schmidt (2006), Mondragón et al. (2009)
grow <- function(nS, nJ, nA, z_A, clim_year, fruited = FALSE,
                 psi0S        = -3.30,
                 psi0J        = -2.70,
                 beta_precip  =  3e-4,
                 beta_rh      =  0.010,
                 sigma        =  0.10,   # individual stochasticity (Raventós 2015)
                 delta_z_base =  0.80,   # mean adult annual growth cm/yr (Zotz 1998)
                 precip_ref   =  2500,   # reference precipitation (mm/yr)
                 rh_ref       =  85,     # reference RH (%)
                 cost_repro   =  0.50,   # Δz multiplier after fruiting (Zotz 1998)
                 z_A_min      =  7.0,
                 z_A_max      = 20.0) {
  nS  <- if (is.na(nS))  0L   else nS
  nJ  <- if (is.na(nJ))  0L   else nJ
  nA  <- if (is.na(nA))  0L   else nA
  z_A <- if (is.na(z_A) || z_A < z_A_min) z_A_min else z_A

  # clim_year has one row per hourly climate record (see get_clim() in
  # get_colonization.R — 576 rows in the current format: 12 months x 24
  # hours x 2 thermal extremes). precip is the mean hourly ERA5 value;
  # scale to annual by × 8760 hours.
  precip_annual <- mean(clim_year$precip, na.rm = TRUE) * 8760
  relhum_mean   <- mean(clim_year$relhum, na.rm = TRUE)

  p_StoJ <- growth_prob("S", precip_annual, relhum_mean, psi0S, psi0J, beta_precip, beta_rh)
  p_JtoA <- growth_prob("J", precip_annual, relhum_mean, psi0S, psi0J, beta_precip, beta_rh)

  epsilon <- rnorm(1, 0, sigma)
  p_StoJ  <- pmin(1, pmax(0, p_StoJ + epsilon))
  p_JtoA  <- pmin(1, pmax(0, p_JtoA + epsilon))

  n_StoJ <- rbinom(1, nS, p_StoJ)
  n_JtoA <- rbinom(1, nJ, p_JtoA)

  z_A_new <- size_increment(z_A, precip_annual, relhum_mean,
                            delta_z_base, precip_ref, rh_ref, sigma,
                            cost_repro, fruited, z_A_min, z_A_max)

  list(nS  = nS - n_StoJ,
       nJ  = nJ + n_StoJ - n_JtoA,
       nA  = nA + n_JtoA,
       z_A = z_A_new)
}

# ── p_est(e): Establishment probability ───────────────────────────────────────
# Probability that a dispersed seed germinates and survives to the seedling stage,
# as a function of local microclimate. Humidity and light are the key covariates.
# "Relative humidity influenced survival of 4/11 species" — Izuddin et al. (2018).
# "Humus presence and microsite (fork) drove survival and growth" — Izuddin (2018).
# Light (swdown) modulates establishment: too little → mycorrhizal fungus absent;
# too much → desiccation. Ratio to site mean captures relative openness.
#
# NOTE (archived version, do not reuse as-is): this version still applies
# p_germ a second time on top of Pass 1's fecundity kernel, and uses a
# deterministic (new_seeds * p_establish) >= 0.1 threshold in the
# non-stochastic branch — both superseded by run_pass2_establish()'s
# corrected, purely-stochastic, single-p_germ-application design.
# Sources: McCormick & Jacquemyn (2014), Izuddin et al. (2018), Zotz (1998)
establish <- function(new_seeds, clim, mean_swdown_site,
                      p_germ = 0.001, stochastic = FALSE) {
  relhum_mean <- mean(clim$relhum, na.rm = TRUE)
  swdown_mean <- mean(clim$swdown[clim$swdown > 0], na.rm = TRUE)
  p_establish <- (relhum_mean / 100) *
    min(swdown_mean / mean_swdown_site, 1) *
    p_germ
  established <- if (stochastic) rbinom(1, new_seeds, p_establish) > 0
  else (new_seeds * p_establish) >= 0.1
  if (established) 1L else 0L
}

# ── Density-dependent growth (unused in main loop — kept for future use) ──────
# Ricker model for within-voxel density regulation if carrying capacity needed.
ricker      <- function(N, r, K) N * exp(r * (1 - N / K))
stochRicker <- function(N, r, K) rpois(1, lambda = ricker(N, r, K))
