# get_colonization.R
# Population dynamics simulation functions for canopy colonization.
# Vital rate functions follow IPM convention (Merow et al. 2013) and are named
# after the quantity they compute. Each cites the source for its form and
# parameter values. Literature parameters are used throughout — individual size
# within each stage class is drawn from a Uniform distribution over the stage
# range (disclosed approximation; we have no individual size census data).
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────

# ── Core equation primitives ──────────────────────────────────────────────────
# Each named after the quantity it computes. Used inside the IPM vital rate
# functions below. Extracting them makes the biology explicit and testable.

# Logistic survival probability for one stage given size and monthly climate.
# logit(s) = β₀_stage + β₁·z + climate_penalties
#
# Light effect on S and J (Izuddin et al. 2018; Zotz 1998):
#   swdown_rel = local swdown / site-mean swdown (computed in run_pass3).
#   Seedlings benefit from increasing light up to site mean (mycorrhizal fungi
#   require some radiation; deep shade reduces germination success).
#   Above 1.5× site mean, desiccation stress reduces seedling survival.
#   Juveniles experience the same light dependence at half the magnitude.
#   Adults are physiologically buffered and retain only the heat-stress term.
# Form:
#   light_S = beta_light_S * min(swdown_rel, 1) - beta_stress_S * max(0, swdown_rel - 1.5)
#   light_J = 0.5 * light_S
# Other penalties:
#   Seedlings: RH-dependent (Zotz 1998 — "most deaths during dry season")
#   Adults:    temperature-dependent — heat stress above 23°C (Olaya-Arenas 2011)
# Sources: Zotz (1998), Izuddin et al. (2018), Olaya-Arenas et al. (2011),
#          Mondragón et al. (2007), Winkler et al. (2009), Zotz & Schmidt (2006)
survival_logit <- function(stage, z, temp, relhum, swdown_rel,
                           beta0S, beta0J, beta0A, beta1,
                           beta_light_S  = 0.30,   # light benefit for seedlings
                           beta_stress_S = 0.20) { # desiccation penalty above 1.5x mean
  light_S <- beta_light_S  * min(swdown_rel, 1.0) -
             beta_stress_S * max(0, swdown_rel - 1.5)
  eta <- switch(stage,
    S = beta0S + beta1 * z - (100 - relhum) / 100 + light_S,
    J = beta0J + beta1 * z + 0.5 * light_S,
    A = beta0A + beta1 * z - max(0, (temp - 23) * 0.05),
    stop("survival_logit: stage must be 'S', 'J', or 'A'"))
  return(1 / (1 + exp(-eta)))
}

# Logistic stage-transition probability given annual precipitation and RH.
# logit(p) = ψ₀_stage + β_precip·P + β_rh·RH
# Sources: Zotz & Schmidt (2006), Mondragón et al. (2009), Izuddin et al. (2018)
transition_logit <- function(psi0, precip_annual_mm, relhum_mean,
                             beta_precip, beta_rh) {
  return(1 / (1 + exp(-(psi0 + beta_precip * precip_annual_mm +
                         beta_rh * relhum_mean))))
}

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

# ── Vital rate functions ──────────────────────────────────────────────────────
# Named following IPM notation: s(z,e), g(z'|z,e), p_r(z), f_s(z), p_est(e)
# z  = pseudobulb length (cm) — state variable
# e  = environmental vector [T, RH, precip, radiation] from microclimf

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

# ── g(z'|z, e): Growth / stage transition probability ─────────────────────────
# Logistic in annual precipitation and mean relative humidity.
# "Annual rainfall significantly affected growth of smaller orchid individuals"
#   — Zotz & Schmidt (2006). "Precipitation increase → faster juvenile growth"
#   — Mondragón et al. (2009).
# At reference conditions (2500 mm/yr, 85% RH) this gives:
#   P(S→J) ≈ 0.18  →  mean time in seedling stage: ~5-6 yr
#   P(J→A) ≈ 0.25  →  mean time to maturity: ~4 yr
# Consistent with "germination to maturity estimated to reach or exceed a decade"
#   — Schmidt & Zotz (2002), cited in Zotz & Schmidt (2006).
# Sources: Zotz & Schmidt (2006), Mondragón et al. (2009), Izuddin et al. (2018)
growth_prob <- function(stage, precip_annual_mm, relhum_mean,
                        psi0S       = -3.30,  # S→J intercept
                        psi0J       = -2.70,  # J→A intercept
                        beta_precip =  3e-4,  # precipitation slope (Zotz 2006)
                        beta_rh     =  0.010) { # RH slope (Izuddin 2018)
  psi0 <- switch(stage,
    S = psi0S, J = psi0J,
    stop("growth_prob: stage must be 'S' or 'J'"))
  return(transition_logit(psi0, precip_annual_mm, relhum_mean, beta_precip, beta_rh))
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

  # clim_year has 48 rows (representative warm + cold hours across the year).
  # precip is the mean hourly ERA5 value; scale to annual by × 8760 hours.
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

# ── p_r(z): Flowering probability ─────────────────────────────────────────────
# Logistic in pseudobulb size z (cm).
# Reproductive minimum PsbL ~7 cm → P(flower | z=7) ≈ 0.50.
# "Reproductive effort increased strongly with plant size" — Zotz (1998).
# Form follows Raventós (2015) and Jacquemyn et al. (2010): logistic on size.
# At z=7 cm: P ≈ 0.50 | z=10 cm: P ≈ 0.82 | z=15 cm: P ≈ 0.98
# Sources: Raventós (2015), Zotz & Schmidt (2006), Zotz (1998), Jacquemyn (2010)
flowering_prob <- function(z,
                           alpha0 = -3.5,  # intercept: threshold near 7 cm
                           alpha1 =  0.5) { # size slope
  1 / (1 + exp(-(alpha0 + alpha1 * z)))
}

# ── f_s(z): Fruit production (conditional on flowering) ───────────────────────
# Poisson mean as function of pseudobulb size z (cm).
# "Larger individuals produced fruits in larger numbers" — Zotz (1998).
# Form: E[fruits | z] = exp(a + b·z) — Poisson regression on size (Raventós 2015).
# At z=7 cm: ~0.9 fruits | z=12 cm: ~1.7 fruits | z=20 cm: ~4.1 fruits
# Sources: Raventós (2015) — Poisson regression of fruit number on size alone
#          Zotz (1998) — positive size effect on fruit number and fruit size
fruit_number <- function(z, a = -1.0, b = 0.12, stochastic = FALSE) {
  mu <- exp(a + b * z)
  if (stochastic) rpois(1, lambda = max(mu, 0)) else mu
}

# ── F(z', z | e): Fecundity kernel ────────────────────────────────────────────
# Full annual seed output: flowering × fruits × pollination × germination × yr-1 survival.
# p_poll: probability a flower is pollinated. Epiphytic orchids are strongly
#   pollinator-limited. "Complete pollination would raise λ to persistence threshold"
#   — Zotz & Schmidt (2006). Default 0.30 (30% of flowers pollinated).
# p_germ: probability a seed germinates. Dust seeds require mycorrhizal partner.
#   "Dispersal not limiting at landscape scale; microsite (mycorrhizal) is"
#   — McCormick & Jacquemyn (2014). Default 0.001 (McCormick & Jacquemyn 2014).
# p_s1: first-year seedling survival. "Fewer than 50% of seedlings survived
#   the first dry season" — Zotz (1998). Default 0.45.
# z is the tracked mean pseudobulb size (cm) for the cell, updated annually by grow().
# Sources: Zotz (1998), Zotz & Schmidt (2006), Raventós (2015),
#          McCormick & Jacquemyn (2014)
reproduce <- function(N, z,
                      p_poll = 0.30,
                      p_germ = 0.001,
                      p_s1   = 0.45,
                      stochastic = FALSE) {
  if (is.na(N) || N == 0) return(0L)
  p_flower <- flowering_prob(z)
  n_fruits <- fruit_number(z, stochastic = stochastic)
  seeds    <- N * p_flower * n_fruits * p_poll * p_germ * p_s1
  if (stochastic) rpois(1, lambda = max(seeds, 0))
  else floor(seeds) + rbinom(1, 1, seeds - floor(seeds))
}

# ── p_est(e): Establishment probability ───────────────────────────────────────
# Probability that a dispersed seed germinates and survives to the seedling stage,
# as a function of local microclimate. Humidity and light are the key covariates.
# "Relative humidity influenced survival of 4/11 species" — Izuddin et al. (2018).
# "Humus presence and microsite (fork) drove survival and growth" — Izuddin (2018).
# Light (swdown) modulates establishment: too little → mycorrhizal fungus absent;
# too much → desiccation. Ratio to site mean captures relative openness.
# p_germ is the baseline mycorrhizal-gated germination probability.
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

# ── d(x'|x): Dispersal kernel ─────────────────────────────────────────────────
# Wind-mediated exponential dispersal with canopy attenuation.
# Mean dispersal distance decreases exponentially with height within canopy
# (Winkler et al. 2009 — "most seeds land near source").
# The canopy attenuation coefficient a = 23(1 - difrac) is derived from the
# fraction of diffuse radiation (difrac), calibrated by microclimf.
# Sources: Winkler et al. (2009), McCormick & Jacquemyn (2014)
# Draw all indN seeds at once — vectorized over seeds, no per-seed loop.
.ind_disperse <- function(x, y, z, indN, winddir, meanDisp, Disp, pad,
                          maxDispZ = 5) {
  wind_rad <- (winddir + 180) %% 360 * pi / 180
  dist  <- pmin(round(rexp(indN, rate = 1 / max(meanDisp, 0.1))), pad)
  angle <- wind_rad + runif(indN, -pi / 4, pi / 4)
  tx <- x + round(dist * sin(angle)) + pad
  ty <- y + round(dist * cos(angle)) + pad
  tz <- z + sample(-maxDispZ:maxDispZ, indN, replace = TRUE) + pad
  D  <- dim(Disp)
  ok <- tx >= 1L & tx <= D[1] & ty >= 1L & ty <= D[2] & tz >= 1L & tz <= D[3]
  if (any(ok)) {
    idx  <- (tx[ok] - 1L) * D[2] * D[3] + (ty[ok] - 1L) * D[3] + tz[ok]
    Disp <- Disp + array(tabulate(idx, nbins = prod(D)), dim = D)
  }
  Disp
}

disperse <- function(x, y, z, seeds, clim, height, canopy_z, a,
                     lambda = 1, Ut = 1, maxDisp = 5, maxDispZ = 5, Disp,
                     stochastic = FALSE) {
  wind     <- mean(clim$windspeed, na.rm = TRUE)
  winddir  <- mean(clim$winddir,   na.rm = TRUE)
  meanDisp <- max(1, min(round((wind * exp(-(a * height / canopy_z))) /
                                 (lambda * Ut)), maxDisp))
  .ind_disperse(x, y, z, indN = seeds, winddir = winddir,
                meanDisp = meanDisp, Disp = Disp, pad = maxDisp,
                maxDispZ = maxDispZ)
}

# ── Density-dependent growth (unused in main loop — kept for future use) ──────
# Ricker model for within-voxel density regulation if carrying capacity needed.
ricker      <- function(N, r, K) N * exp(r * (1 - N / K))
stochRicker <- function(N, r, K) rpois(1, lambda = ricker(N, r, K))

# ── Climate helpers ───────────────────────────────────────────────────────────

# Build the climate data frame for one height tier.
# Spatial mean across the raster at that height → one value per representative
# timestep (24 hours of warmest day + 24 hours of coldest day = 48 rows).
# Macro variables (precip, winddir) come from microenv$.weather (ERA5 hourly),
# summarised to match the 48-row structure by taking the overall mean.
# The result is the single climate object used for all voxels at this height.
get_clim <- function(height, microenv) {
  h_keys <- names(microenv)[names(microenv) != ".spatial" & names(microenv) != ".weather"]
  avail  <- as.numeric(sub("h", "", h_keys))
  h_key  <- h_keys[which.min(abs(avail - height))]
  h      <- microenv[[h_key]]
  if (is.null(h)) return(NULL)

  .smean <- function(arr) {
    if (length(dim(arr)) == 3) apply(arr, 3, mean, na.rm = TRUE)
    else rep(mean(arr, na.rm = TRUE), 24)
  }

  make_df <- function(slot) {
    data.frame(
      temp      = .smean(slot$Tz),
      relhum    = .smean(slot$relhum),
      windspeed = .smean(slot$windspeed),
      swdown    = .smean(slot$Rdirdown) + .smean(slot$Rdifdown),
      difrad    = .smean(slot$Rdifdown)
    )
  }

  df <- rbind(make_df(h$tmax), make_df(h$tmin))  # 48 rows: 24 warm + 24 cold hours

  # append macro variables from ERA5 weather — same value replicated across rows
  w <- microenv$.weather
  if (!is.null(w)) {
    df$precip  <- mean(w$precip,  na.rm = TRUE)
    df$winddir <- mean(w$winddir, na.rm = TRUE)
  } else {
    df$precip  <- NA_real_
    df$winddir <- NA_real_
  }
  df
}

# Return the rows of clim corresponding to a given month.
# tmax rows: 1–24 (hours of warmest day), tmin rows: 25–48 (coldest day).
# Each month maps to two rows: one tmax hour + one tmin hour.
get_clim_month <- function(clim, month) {
  if (is.null(clim) || nrow(clim) == 0) return(NULL)
  clim[c(month, month + 24), ]
}

# ── Forest structure ──────────────────────────────────────────────────────────

# Populate the landscape array with randomly placed trees and derive per-voxel
# carrying capacity from bark surface area (Myster 2017, Johansson 1974).
#
# forestparams must contain:
#   stems_per_ha         tree density for trees ≥10 cm dsh (Myster 2017: ~298/ha)
#   mean_hgt / sd_hgt    tree height distribution (m)
#   mean_crown_r / sd_crown_r  crown radius distribution (m)
#   trunk_r              mean trunk radius in m (Myster 2017: mean dsh 22.7 cm → 0.114 m)
#   branch_density       m² branch surface per m² projected crown area (literature: 2–5)
#   epiphyte_footprint_m2  bark area per individual Maxillariinae (~0.02 m²)
#
# Crown shape: bell curve peaking at 75% of tree height (widest) tapering to
# point at top — approximates tropical montane cloud forest crown architecture.
# Johansson zones 1–2 = trunk only; zones 3–5 = expanding horizontal crown.
build_forest <- function(landscape, heights, forestparams, site_obs, resolution) {
  dims <- dim(landscape)
  xDim <- dims[1]; yDim <- dims[2]; zDim <- dims[3]

  lat_range_m <- (max(site_obs$lat) - min(site_obs$lat)) * 111000
  lon_range_m <- (max(site_obs$lon) - min(site_obs$lon)) * 111000 *
                  cos(mean(site_obs$lat) * pi / 180)
  area_ha <- (lat_range_m * lon_range_m) / 10000
  nTree   <- max(1L, round(area_ha * forestparams$stems_per_ha))

  trees <- data.frame(
    x       = sample(1:xDim, nTree, replace = TRUE),
    y       = sample(1:yDim, nTree, replace = TRUE),
    height  = pmax(1.0, rnorm(nTree, forestparams$mean_hgt,    forestparams$sd_hgt)),
    crown_r = pmax(0.5, rnorm(nTree, forestparams$mean_crown_r, forestparams$sd_crown_r))
  )
  trees$crown_r_cells <- trees$crown_r / resolution

  # voxel height thickness (m) — used for bark surface area calculation
  vox_heights <- diff(c(0, heights))

  landscape[]  <- FALSE
  zone         <- array(0L, dim = c(xDim, yDim, zDim))
  carCap_voxel <- array(1L, dim = c(xDim, yDim, zDim))

  for (ti in seq_len(nTree)) {
    tx <- trees$x[ti];  ty <- trees$y[ti]
    th <- trees$height[ti];  cr <- trees$crown_r_cells[ti]

    x_range <- max(1, tx - ceiling(cr)):min(xDim, tx + ceiling(cr))
    y_range <- max(1, ty - ceiling(cr)):min(yDim, ty + ceiling(cr))

    for (x in x_range) {
      for (y in y_range) {
        horiz_dist <- sqrt((x - tx)^2 + (y - ty)^2)

        for (z in seq_len(zDim)) {
          h <- heights[z]
          if (h > th || h < 0.5) next

          rel_h <- h / th

          jzone <- if      (rel_h < 0.10) 1L
                   else if (rel_h < 0.30) 2L
                   else if (rel_h < 0.50) 3L
                   else if (rel_h < 0.80) 4L
                   else                   5L

          # Trunk zones: only the single column cell; crown zones: bell-shaped radius
          in_tree <- if (jzone <= 2) {
            horiz_dist == 0
          } else {
            crown_fraction <- (rel_h - 0.5) / 0.5       # 0 at zone 3 base, 1 at top
            effective_r    <- cr * sin(crown_fraction * pi)  # peaks at 75% height
            horiz_dist <= effective_r
          }
          if (!in_tree) next

          landscape[x, y, z] <- TRUE
          if (jzone > zone[x, y, z]) zone[x, y, z] <- jzone

          # Bark surface area (m²) per voxel → carrying capacity
          vox_h_m <- vox_heights[z]
          bark_area <- if (jzone <= 2) {
            2 * pi * forestparams$trunk_r * vox_h_m
          } else {
            crown_fraction <- (rel_h - 0.5) / 0.5
            eff_r_m        <- cr * resolution * sin(crown_fraction * pi)  # grid cells → m
            pi * eff_r_m^2 * forestparams$branch_density * vox_h_m
          }
          cap <- max(1L, as.integer(floor(bark_area / forestparams$epiphyte_footprint_m2)))
          if (cap > carCap_voxel[x, y, z]) carCap_voxel[x, y, z] <- cap
        }
      }
    }
  }

  log_msg(sprintf("build_forest: %d trees | %.1f ha | valid voxels: %d | mean carCap: %.1f",
                  nTree, area_ha, sum(landscape), mean(carCap_voxel[landscape])))
  list(landscape = landscape, zone = zone, carCap_voxel = carCap_voxel,
       trees = trees, n_trees = nTree, area_ha = area_ha)
}

# ── Simulation setup ──────────────────────────────────────────────────────────

init_colonization <- function(site, niches, canopy_grid, microenv,
                              resolution = 10, carCap = 1, maxDisp = 5, params,
                              forestparams = NULL, allsites = FALSE) {
  site_name <- site$Site
  h_keys    <- names(microenv)[!names(microenv) %in% c(".spatial", ".weather")]
  heights   <- sort(as.numeric(sub("h", "", h_keys)))
  site_obs  <- if (allsites) niches else niches[niches$Area_or_Site == site_name, ]
  zDim <- length(heights)
  lat_range_m <- (max(site_obs$lat) - min(site_obs$lat)) * 111000
  lon_range_m <- (max(site_obs$lon) - min(site_obs$lon)) * 111000 *
    cos(mean(site_obs$lat) * pi / 180)
  xDim <- max(round(lon_range_m / resolution), 10) + 4
  yDim <- max(round(lat_range_m / resolution), 10) + 4

  species_ids <- sort(unique(site_obs$FinalID))
  n_species   <- length(species_ids)
  sp_index    <- setNames(seq_along(species_ids), species_ids)

  log_msg(sprintf(
    "Landscape: %d x %d x %d cells (%.0f x %.0f m, %.1f-%.1fm) | %d species",
    xDim, yDim, zDim,
    xDim * resolution, yDim * resolution,
    min(heights), max(heights), n_species))

  resolution_deg <- resolution / 111000
  lon_min <- min(site_obs$lon) - resolution_deg
  lat_min <- min(site_obs$lat) - resolution_deg

  coord_to_idx <- function(lon, lat, h) {
    c(x = min(xDim, max(1, round((lon - lon_min) / resolution_deg) + 1)),
      y = min(yDim, max(1, round((lat - lat_min) / resolution_deg) + 1)),
      z = which.min(abs(heights - h)))
  }

  landscape <- array(FALSE, dim = c(xDim, yDim, zDim))

  if (!is.null(forestparams)) {
    # Tree-based landscape: place stems randomly, derive zone and carCap from geometry
    forest       <- build_forest(landscape, heights, forestparams, site_obs, resolution)
    landscape    <- forest$landscape
    zone         <- forest$zone
    carCap_voxel <- forest$carCap_voxel
  } else {
    # Fallback: flat canopy grid ceiling, uniform carrying capacity
    cg_xDim <- nrow(canopy_grid); cg_yDim <- ncol(canopy_grid)
    for (x in 1:xDim) for (y in 1:yDim) for (z in 1:zDim) {
      cx <- min(x, cg_xDim); cy <- min(y, cg_yDim)
      landscape[x, y, z] <- heights[z] <= canopy_grid[cx, cy] && heights[z] >= 0.5
    }
    zone         <- array(0L, dim = c(xDim, yDim, zDim))
    carCap_voxel <- array(as.integer(carCap), dim = c(xDim, yDim, zDim))
  }

  log_msg("Pre-computing climate lookup table from microenv...")
  clim_by_height       <- lapply(heights, function(h) get_clim(h, microenv))
  clim_month_by_height <- lapply(clim_by_height, function(clim)
    lapply(1:12, function(m) get_clim_month(clim, m)))
  log_msg("Climate lookup ready.")

  valid_clim <- which(!sapply(clim_by_height, is.null))
  mid_zi     <- valid_clim[which.min(abs(heights[valid_clim] -
                              max(min(heights), mean(canopy_grid, na.rm=TRUE)/2)))]
  if (length(mid_zi) == 0) mid_zi <- valid_clim[1]
  clim_mid <- clim_by_height[[mid_zi]]
  difrac   <- mean(clim_mid$difrad / (clim_mid$swdown + 0.001), na.rm = TRUE)
  a        <- 23 * (1 - difrac)

  med_zi           <- valid_clim[which.min(abs(heights[valid_clim] -
                                               median(heights[valid_clim])))]
  mean_swdown_site <- mean(clim_by_height[[med_zi]]$swdown, na.rm = TRUE)

  log_msg(sprintf("Valid climate heights: %d/%d | mean_swdown: %.1f",
                  length(valid_clim), zDim, mean_swdown_site))

  params$a                <- a
  params$mean_swdown_site <- mean_swdown_site
  # ensure size-tracking defaults exist if caller omitted them
  if (is.null(params$z_A_min))      params$z_A_min      <- 7.0
  if (is.null(params$z_A_max))      params$z_A_max      <- 20.0
  if (is.null(params$delta_z_base)) params$delta_z_base <- 0.80
  if (is.null(params$cost_repro))   params$cost_repro   <- 0.50

  list(
    site_name            = site_name,
    site_obs             = site_obs,
    heights              = heights,
    xDim = xDim, yDim = yDim, zDim = zDim,
    n_species            = n_species,
    species_ids          = species_ids,
    sp_index             = sp_index,
    landscape            = landscape,
    zone                 = zone,
    carCap_voxel         = carCap_voxel,
    coord_to_idx         = coord_to_idx,
    clim_by_height       = clim_by_height,
    clim_month_by_height = clim_month_by_height,
    params               = params,
    carCap               = carCap,
    maxDisp              = maxDisp,
    dispersalmatrix      = array(0L, dim = c(xDim + 2*maxDisp,
                                             yDim + 2*maxDisp,
                                             zDim + 2*maxDisp,
                                             n_species))
  )
}

# ── Simulation passes ─────────────────────────────────────────────────────────

# Pass 1: Adults reproduce and seeds are dispersed through the 3D landscape.
# Uses tracked pseudobulb size (size_A[x,y,z,t,sp]) for the F kernel.
# Returns both the dispersal matrix and a fruited boolean array so pass 3
# can apply the cost-of-reproduction penalty to the right cells.
run_pass1_disperse <- function(state, abundanceA, size_A, t) {
  p       <- state$params
  Disp    <- array(0L, dim = c(
    state$xDim + 2 * state$maxDisp,
    state$yDim + 2 * state$maxDisp,
    state$zDim + 2 * state$maxDisp,
    state$n_species))
  fruited <- array(FALSE, dim = c(state$xDim, state$yDim, state$zDim, state$n_species))
  total_seeds <- 0L

  for (sp in 1:state$n_species) {
    nonzero <- which(abundanceA[,,,t,sp] > 0, arr.ind = TRUE)
    for (i in seq_len(nrow(nonzero))) {
      x <- nonzero[i,1]; y <- nonzero[i,2]; z <- nonzero[i,3]
      if (!state$landscape[x, y, z]) next
      N <- abundanceA[x, y, z, t, sp]
      if (is.na(N) || N == 0) next
      psb_z <- size_A[x, y, z, t, sp]  # tracked pseudobulb size (cm)
      seeds <- reproduce(N, z = psb_z,
                         p_poll = p$p_poll, p_germ = p$p_germ, p_s1 = p$p_s1)
      total_seeds <- total_seeds + seeds
      if (seeds == 0) next
      fruited[x, y, z, sp] <- TRUE
      Disp[,,,sp] <- disperse(
        x = x, y = y, z = z, seeds = seeds,
        clim    = state$clim_by_height[[z]],
        height  = state$heights[z], canopy_z = p$canopy_z,
        a = p$a, lambda = p$lambda, Ut = p$Ut,
        maxDisp = state$maxDisp, Disp = Disp[,,,sp])
    }
  }
  message("Pass1: seeds produced=", total_seeds, " dispersed=", sum(Disp))
  list(Disp = Disp, fruited = fruited)
}

# Pass 2: Dispersed seeds try to establish in new cells.
# Vectorized per height: p_establish is a scalar per height tier (climate is
# spatially uniform at each z), so one rbinom() call handles the full [xDim,yDim]
# slice instead of looping voxel-by-voxel.
run_pass2_establish <- function(state, abundanceS, abundanceJ, abundanceA,
                                dispersalmatrix, t) {
  p   <- state$params
  pad <- state$maxDisp
  tnext <- t + 1L
  if (tnext > dim(abundanceS)[4]) return(abundanceS)

  # Total occupancy across all species — needed for carCap check
  total_occ <- apply(abundanceS[,,,t,,drop=FALSE], 1:3, sum) +
               apply(abundanceJ[,,,t,,drop=FALSE], 1:3, sum) +
               apply(abundanceA[,,,t,,drop=FALSE], 1:3, sum)
  total_seeds_seen <- 0L; total_established <- 0L

  for (sp in seq_len(state$n_species)) {
    for (zi in seq_len(state$zDim)) {
      clim <- state$clim_by_height[[zi]]
      if (is.null(clim)) next

      # Extract the [xDim, yDim] seed-rain slice (trim padding)
      seeds_slice <- dispersalmatrix[(pad+1L):(pad+state$xDim),
                                     (pad+1L):(pad+state$yDim),
                                     zi + pad, sp]
      if (sum(seeds_slice) == 0L) next

      # p_establish: scalar for this height (spatial mean already in clim)
      relhum_mean <- mean(clim$relhum, na.rm = TRUE)
      swdown_mean <- mean(clim$swdown[clim$swdown > 0], na.rm = TRUE)
      p_est <- (relhum_mean / 100) *
               min(swdown_mean / p$mean_swdown_site, 1) *
               p$p_germ

      # Mask: valid landscape, not at capacity, has seeds
      can_est <- state$landscape[,,zi] &
                 total_occ[,,zi] < state$carCap_voxel[,,zi] &
                 seeds_slice > 0L

      if (!any(can_est)) next

      total_seeds_seen <- total_seeds_seen + sum(seeds_slice[can_est])
      n  <- sum(can_est)
      established <- rbinom(n, as.integer(seeds_slice[can_est]), p_est)
      # Clamp to remaining capacity
      space <- pmax(0L, state$carCap_voxel[,,zi][can_est] - total_occ[,,zi][can_est])
      established <- pmin(as.integer(established), space)
      abundanceS[,,zi,tnext,sp][can_est] <-
        abundanceS[,,zi,tnext,sp][can_est] + established
      total_established <- total_established + sum(established)
    }
  }
  message("Pass2: seeds seen=", total_seeds_seen, " established=", total_established)
  abundanceS
}

# Pass 3: Survival and stage transitions, vectorized per height tier.
#
# Key insight: climate is spatially uniform at each height (get_clim returns
# a spatial mean), so survival probability is a *scalar* per (height, month,
# stage). One rbinom() call over the full [xDim × yDim] slice replaces the
# old loop over every occupied voxel — same biology, far fewer R calls.
#
# size_A is also updated in bulk: delta_z is drawn for every cell in the slice
# and masked post-hoc to cells that actually have adults.
run_pass3_survive_grow <- function(state, abundanceS, abundanceJ, abundanceA,
                                   size_A, fruited, t) {
  p    <- state$params
  xDim <- state$xDim; yDim <- state$yDim

  # Helper: apply rbinom to a 2D slice with a scalar probability.
  # Cells outside the landscape are 0 already — no masking needed.
  .surv_slice <- function(sl, prob) {
    if (prob <= 0 || sum(sl) == 0L) return(sl * 0L)
    if (prob >= 1) return(sl)
    array(rbinom(length(sl), as.integer(sl), prob), dim = dim(sl))
  }

  for (sp in seq_len(state$n_species)) {

    # ── s(z, e): monthly survival ─────────────────────────────────────────────
    # Loop over heights × months — draw one z per stage (same approximation as
    # before; individual size data not available) and apply vectorised rbinom.
    for (zi in seq_len(state$zDim)) {
      if (sum(abundanceS[,,zi,t,sp]) + sum(abundanceJ[,,zi,t,sp]) +
          sum(abundanceA[,,zi,t,sp]) == 0L) next

      for (month in 1:12) {
        cm <- state$clim_month_by_height[[zi]][[month]]
        if (is.null(cm) || nrow(cm) == 0) next
        temp        <- mean(cm$temp,   na.rm = TRUE)
        relhum      <- mean(cm$relhum, na.rm = TRUE)
        swdown_mean <- mean(cm$swdown[cm$swdown > 0], na.rm = TRUE)
        swdown_rel  <- if (!is.na(swdown_mean) && p$mean_swdown_site > 0)
                         swdown_mean / p$mean_swdown_site else 1.0

        z_S <- runif(1, p$z_S_min, p$z_S_max)
        z_J <- runif(1, p$z_J_min, p$z_J_max)
        # Use mean tracked adult size for the height slice (better than fixed draw)
        z_A_mean <- mean(size_A[,,zi,t,sp][state$landscape[,,zi]], na.rm = TRUE)
        if (is.na(z_A_mean) || z_A_mean < p$z_A_min) z_A_mean <- p$z_A_min

        s_S <- survival_logit("S", z_S,      temp, relhum, swdown_rel, p$beta0S, p$beta0J, p$beta0A, p$beta1)
        s_J <- survival_logit("J", z_J,      temp, relhum, swdown_rel, p$beta0S, p$beta0J, p$beta0A, p$beta1)
        s_A <- survival_logit("A", z_A_mean, temp, relhum, swdown_rel, p$beta0S, p$beta0J, p$beta0A, p$beta1)

        abundanceS[,,zi,t,sp] <- .surv_slice(abundanceS[,,zi,t,sp], s_S)
        abundanceJ[,,zi,t,sp] <- .surv_slice(abundanceJ[,,zi,t,sp], s_J)
        abundanceA[,,zi,t,sp] <- .surv_slice(abundanceA[,,zi,t,sp], s_A)
      }
    }

    # ── g(z'|z, e): annual stage transitions + adult size update ──────────────
    for (zi in seq_len(state$zDim)) {
      slS <- abundanceS[,,zi,t,sp]
      slJ <- abundanceJ[,,zi,t,sp]
      slA <- abundanceA[,,zi,t,sp]
      if (sum(slS) + sum(slJ) + sum(slA) == 0L) next

      clim_year <- state$clim_by_height[[zi]]
      if (is.null(clim_year)) next

      precip_annual <- mean(clim_year$precip, na.rm = TRUE) * 8760
      relhum_mean   <- mean(clim_year$relhum, na.rm = TRUE)

      p_StoJ <- growth_prob("S", precip_annual, relhum_mean,
                            p$psi0S, p$psi0J, p$beta_precip, p$beta_rh)
      p_JtoA <- growth_prob("J", precip_annual, relhum_mean,
                            p$psi0S, p$psi0J, p$beta_precip, p$beta_rh)
      # Add shared stochastic perturbation (same epsilon for all cells at this height)
      epsilon <- rnorm(1, 0, p$sigma)
      p_StoJ  <- pmin(1, pmax(0, p_StoJ + epsilon))
      p_JtoA  <- pmin(1, pmax(0, p_JtoA + epsilon))

      n_StoJ <- .surv_slice(slS, p_StoJ)
      n_JtoA <- .surv_slice(slJ, p_JtoA)

      abundanceS[,,zi,t,sp] <- slS - n_StoJ
      abundanceJ[,,zi,t,sp] <- slJ + n_StoJ - n_JtoA
      abundanceA[,,zi,t,sp] <- slA + n_JtoA

      # ── Adult size update (vectorised over xy) ────────────────────────────────
      za_slice <- size_A[,,zi,t,sp]
      # Draw delta_z for every cell; mask to adult-occupied cells afterwards
      delta_z <- p$delta_z_base *
                 (precip_annual / 2500) *
                 (relhum_mean   / 85) +
                 rnorm(xDim * yDim, 0, p$sigma * 0.5)
      dim(delta_z) <- c(xDim, yDim)
      # Cost of reproduction: halve delta_z where the cell fruited this year
      fr_slice <- fruited[,,zi,sp]
      delta_z[fr_slice] <- delta_z[fr_slice] * p$cost_repro
      za_new <- pmin(p$z_A_max, pmax(p$z_A_min, za_slice + delta_z))
      # Cells newly promoted from J (no prior adults) start at z_A_min
      fresh <- n_JtoA > 0L & slA == 0L
      za_new[fresh] <- p$z_A_min
      # Only write back where adults exist (preserve z_A_min in empty cells)
      has_adult <- abundanceA[,,zi,t,sp] > 0L
      za_slice[has_adult] <- za_new[has_adult]
      size_A[,,zi,t,sp] <- za_slice
    }
  }
  list(S = abundanceS, J = abundanceJ, A = abundanceA, size_A = size_A)
}

# ── Internal live visualisation ───────────────────────────────────────────────

.plot_live <- function(state, abundanceS, abundanceJ, abundanceA,
                       totalS, totalJ, totalA, t, carCap, sleeptime = 0.2) {
  stage_cols <- scico::scico(3, palette = "lipari", begin = 0.2, end = 0.8)
  sp_cols    <- if (state$n_species == 1)
    scico::scico(3, palette = "lipari", begin = 0.3, end = 0.7)[2]
  else
    scico::scico(state$n_species, palette = "lipari", begin = 0.2, end = 0.8)
  total  <- totalS + totalJ + totalA
  n_sp   <- state$n_species
  par(mfrow = c(1, n_sp + 1), mar = c(4, 4, 3, 2))
  for (sp in 1:n_sp) {
    ts_S <- sapply(1:t, function(i) sum(abundanceS[,,,i,sp]))
    ts_J <- sapply(1:t, function(i) sum(abundanceJ[,,,i,sp]))
    ts_A <- sapply(1:t, function(i) sum(abundanceA[,,,i,sp]))
    ts_total <- ts_S + ts_J + ts_A
    plot(ts_total, type="b", col=sp_cols[sp], lwd=2,
         ylim=c(0, max(ts_total, 1)), xlab="Year", ylab="Abundance",
         main=paste0(state$species_ids[sp], " (t=", t, ")"), las=1)
    lines(ts_S, type="b", col=stage_cols[1], pch=16, lty=2)
    lines(ts_J, type="b", col=stage_cols[2], pch=17, lty=2)
    lines(ts_A, type="b", col=stage_cols[3], pch=15, lty=2)
    legend("topleft", legend=c("Total","S","J","A"),
           col=c(sp_cols[sp], stage_cols), lty=c(1,2,2,2),
           pch=c(NA,16,17,15), cex=0.6)
  }
  plot(total[1:t], type="b", col="black", lwd=2,
       ylim=c(0, max(total, 1)), xlab="Year", ylab="Abundance",
       main=paste0("All species (t=", t, ")"), las=1)
  lines(totalS[1:t], type="b", col=stage_cols[1], pch=16)
  lines(totalJ[1:t], type="b", col=stage_cols[2], pch=17)
  lines(totalA[1:t], type="b", col=stage_cols[3], pch=15)
  abline(h=carCap * state$xDim * state$yDim * state$zDim, col="red", lty=2)
  legend("topleft", legend=c("Total","S","J","A"),
         col=c("black", stage_cols), lty=1, pch=c(NA,16,17,15), cex=0.6)
  dev.flush(); Sys.sleep(sleeptime)
}

# ── Post-hoc abundance plot ───────────────────────────────────────────────────

plot_abundance <- function(result, t = NULL, species_specific = TRUE) {
  state      <- result$state
  t_max      <- if (is.null(t)) length(result$totalabundanceA) else t
  stage_cols <- scico::scico(3, palette = "lipari", begin = 0.2, end = 0.8)
  sp_cols    <- if (state$n_species == 1)
    scico::scico(3, palette = "lipari", begin = 0.3, end = 0.7)[2]
  else
    scico::scico(state$n_species, palette = "lipari", begin = 0.2, end = 0.8)
  totalS <- result$totalabundanceS; totalJ <- result$totalabundanceJ
  totalA <- result$totalabundanceA; total  <- totalS + totalJ + totalA
  n_panels <- if (species_specific) state$n_species + 1L else 1L
  par(mfrow = c(1, n_panels), mar = c(4, 4, 3, 2))
  if (species_specific) {
    for (sp in 1:state$n_species) {
      ts_S <- sapply(1:t_max, function(i) sum(result$abundanceS[,,,i,sp]))
      ts_J <- sapply(1:t_max, function(i) sum(result$abundanceJ[,,,i,sp]))
      ts_A <- sapply(1:t_max, function(i) sum(result$abundanceA[,,,i,sp]))
      ts_total <- ts_S + ts_J + ts_A
      plot(ts_total, type="b", col=sp_cols[sp], lwd=2,
           ylim=c(0, max(ts_total, 1)), xlab="Year", ylab="Abundance",
           main=paste0(state$species_ids[sp], " (t=", t_max, ")"), las=1)
      lines(ts_S, type="b", col=stage_cols[1], pch=16, lty=2)
      lines(ts_J, type="b", col=stage_cols[2], pch=17, lty=2)
      lines(ts_A, type="b", col=stage_cols[3], pch=15, lty=2)
      legend("topleft", legend=c("Total","S","J","A"),
             col=c(sp_cols[sp], stage_cols), lty=c(1,2,2,2),
             pch=c(NA,16,17,15), cex=0.6)
    }
  }
  plot(total[1:t_max], type="b", col="black", lwd=2,
       ylim=c(0, max(total[1:t_max], 1)), xlab="Year", ylab="Abundance",
       main=paste0(state$site_name, " — All species (t=", t_max, ")"), las=1)
  lines(totalS[1:t_max], type="b", col=stage_cols[1], pch=16)
  lines(totalJ[1:t_max], type="b", col=stage_cols[2], pch=17)
  lines(totalA[1:t_max], type="b", col=stage_cols[3], pch=15)
  abline(h=state$carCap * state$xDim * state$yDim * state$zDim, col="red", lty=2)
  legend("topleft", legend=c("Total","S","J","A"),
         col=c("black", stage_cols), lty=1, pch=c(NA,16,17,15), cex=0.6)
}

# ── 3D post-hoc visualisation ─────────────────────────────────────────────────

plot_3d_abundance <- function(result, t = NULL) {
  state <- result$state
  if (is.null(t)) t <- dim(result$abundanceA)[4]
  sp_cols <- if (state$n_species == 1)
    scico::scico(3, palette="lipari", begin=0.3, end=0.7)[2]
  else
    scico::scico(state$n_species, palette="lipari", begin=0.2, end=0.8)
  rows <- list()
  for (sp in 1:state$n_species) {
    sp_name <- state$species_ids[sp]; col <- sp_cols[sp]
    for (stg in list(list(arr=result$abundanceS, nm="S", sym="circle"),
                     list(arr=result$abundanceJ, nm="J", sym="diamond"),
                     list(arr=result$abundanceA, nm="A", sym="square"))) {
      idx <- which(stg$arr[,,,t,sp] > 0, arr.ind=TRUE)
      if (nrow(idx) > 0)
        rows[[length(rows)+1]] <- data.frame(
          x=idx[,1], y=idx[,2], z=idx[,3],
          species=sp_name, stage=stg$nm, color=col, symbol=stg$sym,
          stringsAsFactors=FALSE)
    }
  }
  if (length(rows) == 0) { message("No individuals to plot at t=", t); return(invisible(NULL)) }
  df     <- do.call(rbind, rows)
  traces <- split(df, paste0(df$species, "_", df$stage))
  fig    <- plotly::plot_ly()
  for (tr in traces)
    fig <- plotly::add_trace(fig, data=tr, x=~x, y=~y, z=~z,
                             type="scatter3d", mode="markers",
                             name=paste0(tr$species[1], " ", tr$stage[1]),
                             marker=list(symbol=tr$symbol[1], color=tr$color[1],
                                         size=6, opacity=0.85))
  fig <- plotly::layout(fig, title=paste0("Abundance (t=", t, ")"),
                        scene=list(xaxis=list(title="x"), yaxis=list(title="y"),
                                   zaxis=list(title="height tier")))
  print(fig); invisible(fig)
}

# ── Spin-up ───────────────────────────────────────────────────────────────────

run_spinup <- function(state, n_gens = 5, Visualize = TRUE,
                       carCap = 1, sleeptime = 0.2, visualize_dispersion = FALSE) {
  p <- state$params
  xDim <- state$xDim; yDim <- state$yDim; zDim <- state$zDim
  n_species <- state$n_species
  spinupS  <- array(0L,  dim=c(xDim, yDim, zDim, n_gens, n_species))
  spinupJ  <- array(0L,  dim=c(xDim, yDim, zDim, n_gens, n_species))
  spinupA  <- array(0L,  dim=c(xDim, yDim, zDim, n_gens, n_species))
  # size_A: mean pseudobulb length (cm) per adult cell; initialised at z_A_min
  size_A   <- array(p$z_A_min, dim=c(xDim, yDim, zDim, n_gens, n_species))
  totalS   <- numeric(n_gens); totalJ <- numeric(n_gens); totalA <- numeric(n_gens)

  for (i in 1:nrow(state$site_obs)) {
    idx <- state$coord_to_idx(state$site_obs$lon[i], state$site_obs$lat[i],
                              state$site_obs$Height_m[i])
    x <- idx[1]; y <- idx[2]; z <- idx[3]
    sp <- state$sp_index[state$site_obs$FinalID[i]]
    if (x >= 1 && x <= xDim && y >= 1 && y <= yDim && state$landscape[x, y, z]) {
      spinupA[x, y, z, 1, sp] <- 1L
      size_A[x, y, z, 1, sp]  <- p$z_A_min
    }
  }
  log_msg(sprintf("Spin-up: placed %d observed individuals (%d species)",
                  sum(spinupA[,,,1,]), n_species))

  fruited <- array(FALSE, dim=c(xDim, yDim, zDim, n_species))
  Disp    <- NULL

  for (gen in 1:n_gens) {
    log_msg(sprintf("Spin-up generation %d/%d", gen, n_gens))
    pass1   <- run_pass1_disperse(state, spinupA, size_A, gen)
    Disp    <- pass1$Disp; fruited <- pass1$fruited
    if (visualize_dispersion) {
      fields::image.plot(apply(Disp, c(1,2), sum),
                         col=scico::scico(25, palette="lajolla"),
                         main=paste0("Dispersal (gen ", gen, ") | seeds: ", sum(Disp)),
                         xlab="x", ylab="y")
      dev.flush(); Sys.sleep(sleeptime)
    }
    spinupS <- run_pass2_establish(state, spinupS, spinupJ, spinupA, Disp, gen)
    result  <- run_pass3_survive_grow(state, spinupS, spinupJ, spinupA,
                                      size_A, fruited, gen)
    spinupS <- result$S; spinupJ <- result$J
    spinupA <- result$A; size_A  <- result$size_A
    if (gen < n_gens) {
      spinupS[,,,gen+1,] <- spinupS[,,,gen+1,] + spinupS[,,,gen,]
      spinupJ[,,,gen+1,] <- spinupJ[,,,gen+1,] + spinupJ[,,,gen,]
      spinupA[,,,gen+1,] <- spinupA[,,,gen+1,] + spinupA[,,,gen,]
      size_A[,,,gen+1,]  <- size_A[,,,gen,]   # carry size forward
    }
    totalS[gen] <- sum(spinupS[,,,gen,]); totalJ[gen] <- sum(spinupJ[,,,gen,])
    totalA[gen] <- sum(spinupA[,,,gen,])
    if (Visualize)
      .plot_live(state, spinupS, spinupJ, spinupA,
                 totalS, totalJ, totalA, gen, carCap, sleeptime)
    log_msg(sprintf("Gen %d: S=%d J=%d A=%d total=%d", gen,
                    totalS[gen], totalJ[gen], totalA[gen],
                    totalS[gen]+totalJ[gen]+totalA[gen]))
  }
  log_msg(sprintf("Spin-up complete: S=%d J=%d A=%d",
                  sum(spinupS[,,,n_gens,]), sum(spinupJ[,,,n_gens,]),
                  sum(spinupA[,,,n_gens,])))
  list(S=spinupS[,,,n_gens,], J=spinupJ[,,,n_gens,], A=spinupA[,,,n_gens,],
       size_A=size_A[,,,n_gens,], last_disp=Disp)
}

# ── Main wrapper ──────────────────────────────────────────────────────────────

runcolonization <- function(site, niches, canopy_grid, microenv,
                            timesteps=50, resolution=10, carCap=1,
                            maxDisp=5, stochastic=FALSE,
                            Visualize=TRUE, sleeptime=0.2,
                            visualize_dispersion=FALSE, spinup=5,
                            parameters, forestparams=NULL, allsites=FALSE,
                            train_frac=0.70, seed=42) {
  set.seed(seed)

  # Per-species train/validation split — train_frac of observations per species
  # are used for spin-up; the remainder are held out for validation.
  site_obs_all <- if (allsites) niches else niches[niches$Area_or_Site == site, ]
  train_idx <- unlist(lapply(split(seq_len(nrow(site_obs_all)),
                                   site_obs_all$FinalID),
                             function(idx) sample(idx, max(1L, round(length(idx) * train_frac)))))
  niches_train <- site_obs_all[ train_idx, ]
  niches_val   <- site_obs_all[-train_idx, ]
  log_msg(sprintf("Train/val split: %d train | %d val observations (train_frac=%.2f)",
                  nrow(niches_train), nrow(niches_val), train_frac))

  state     <- init_colonization(site, niches_train, canopy_grid, microenv,
                                 resolution, carCap, maxDisp, params=parameters,
                                 forestparams=forestparams, allsites=allsites)
  xDim      <- state$xDim; yDim <- state$yDim; zDim <- state$zDim
  n_species <- state$n_species
  abundanceS      <- array(0L,          dim=c(xDim, yDim, zDim, timesteps, n_species))
  abundanceJ      <- array(0L,          dim=c(xDim, yDim, zDim, timesteps, n_species))
  abundanceA      <- array(0L,          dim=c(xDim, yDim, zDim, timesteps, n_species))
  size_A          <- array(parameters$z_A_min, dim=c(xDim, yDim, zDim, timesteps, n_species))
  totalabundanceS <- numeric(timesteps)
  totalabundanceJ <- numeric(timesteps)
  totalabundanceA <- numeric(timesteps)

  sp <- run_spinup(state, Visualize=Visualize, n_gens=spinup,
                   carCap=carCap, sleeptime=sleeptime,
                   visualize_dispersion=visualize_dispersion)
  abundanceS[,,,1,] <- sp$S; abundanceJ[,,,1,] <- sp$J
  abundanceA[,,,1,] <- sp$A; size_A[,,,1,]     <- sp$size_A
  totalabundanceS[1] <- sum(abundanceS[,,,1,])
  totalabundanceJ[1] <- sum(abundanceJ[,,,1,])
  totalabundanceA[1] <- sum(abundanceA[,,,1,])
  log_msg(sprintf("Starting population: %d S, %d J, %d A",
                  totalabundanceS[1], totalabundanceJ[1], totalabundanceA[1]))

  fruited <- array(FALSE, dim=c(xDim, yDim, zDim, n_species))
  Disp    <- sp$last_disp

  for (t in 1:(timesteps-1)) {
    pass1   <- run_pass1_disperse(state, abundanceA, size_A, t)
    Disp    <- pass1$Disp; fruited <- pass1$fruited
    abundanceS <- run_pass2_establish(state, abundanceS, abundanceJ, abundanceA, Disp, t)
    result  <- run_pass3_survive_grow(state, abundanceS, abundanceJ, abundanceA,
                                      size_A, fruited, t)
    abundanceS <- result$S; abundanceJ <- result$J
    abundanceA <- result$A; size_A     <- result$size_A
    abundanceS[,,,t+1,] <- abundanceS[,,,t+1,] + abundanceS[,,,t,]
    abundanceJ[,,,t+1,] <- abundanceJ[,,,t+1,] + abundanceJ[,,,t,]
    abundanceA[,,,t+1,] <- abundanceA[,,,t+1,] + abundanceA[,,,t,]
    size_A[,,,t+1,]     <- size_A[,,,t,]        # carry size to next timestep
    totalabundanceS[t+1] <- sum(abundanceS[,,,t+1,])
    totalabundanceJ[t+1] <- sum(abundanceJ[,,,t+1,])
    totalabundanceA[t+1] <- sum(abundanceA[,,,t+1,])
    if (Visualize)
      .plot_live(state, abundanceS, abundanceJ, abundanceA,
                 totalabundanceS, totalabundanceJ, totalabundanceA,
                 t+1, carCap, sleeptime)
    log_msg(sprintf("t=%d | S=%d J=%d A=%d total=%d", t+1,
                    totalabundanceS[t+1], totalabundanceJ[t+1], totalabundanceA[t+1],
                    totalabundanceS[t+1]+totalabundanceJ[t+1]+totalabundanceA[t+1]))
  }
  list(landscape=state$landscape,
       abundanceS=abundanceS, abundanceJ=abundanceJ, abundanceA=abundanceA,
       size_A=size_A,
       totalabundanceS=totalabundanceS, totalabundanceJ=totalabundanceJ,
       totalabundanceA=totalabundanceA,
       heights=state$heights, species_ids=state$species_ids,
       xDim=xDim, yDim=yDim, zDim=zDim, n_species=n_species,
       state=state, last_disp=Disp,
       obs_train=niches_train, obs_val=niches_val)
}
