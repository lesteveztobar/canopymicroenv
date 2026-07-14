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
# Applied once per month and compounded over 12 months (run_pass3_survive_
# grow()) — survival must succeed every month, so beta0S/J/A are calibrated
# such that p_month = p_annual^(1/12), i.e. 12 reference-condition months
# compound back to the literature annual target (S: 0.44, J: 0.60, A: 0.85
# at default beta0A). Passing an annual-target-calibrated intercept here
# directly (logit^-1(beta0) == p_annual) would silently compound it 12x too
# many times and crash the population almost every year — see methods.tex
# for the derivation and the corresponding shifted parameter table.
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

# size_increment() archived to archive.R — superseded by the vectorized
# adult-size update inlined in run_pass3_survive_grow().

# ── Vital rate functions ──────────────────────────────────────────────────────
# Named following IPM notation: s(z,e), g(z'|z,e), p_r(z), f_s(z), p_est(e)
# z  = pseudobulb length (cm) — state variable
# e  = environmental vector [T, RH, precip, radiation] from microclimf

# survive() archived to archive.R — superseded by the vectorized per-height
# survival draw (.surv_slice()) inlined in run_pass3_survive_grow(). Its
# core logic (survival_logit()) is still live and called directly there.

# ── g(z'|z, e): Growth / stage transition probability ─────────────────────────
# Logistic in annual precipitation and mean relative humidity, evaluated
# monthly and compounded across the year (run_pass3_survive_grow()) — a
# transition is an "at least once this year" event, so the intercepts below
# are calibrated so the monthly complement (1 - q_month) compounds to the
# annual target: q_month = 1 - (1 - p_annual)^(1/12). At reference conditions
# (2500 mm/yr, 85% RH) this gives:
#   P(S→J) ≈ 0.18/yr (q_month ≈ 0.0164)  →  mean time in seedling stage: ~5-6 yr
#   P(J→A) ≈ 0.25/yr (q_month ≈ 0.0237)  →  mean time to maturity: ~4 yr
# Consistent with "germination to maturity estimated to reach or exceed a decade"
#   — Schmidt & Zotz (2002), cited in Zotz & Schmidt (2006).
# "Annual rainfall significantly affected growth of smaller orchid individuals"
#   — Zotz & Schmidt (2006). "Precipitation increase → faster juvenile growth"
#   — Mondragón et al. (2009).
# Sources: Zotz & Schmidt (2006), Mondragón et al. (2009), Izuddin et al. (2018)
growth_prob <- function(stage, precip_annual_mm, relhum_mean,
                        psi0S       = -5.877,  # S→J intercept (monthly-compounded)
                        psi0J       = -5.319,  # J→A intercept (monthly-compounded)
                        beta_precip =  3e-4,  # precipitation slope (Zotz 2006)
                        beta_rh     =  0.010) { # RH slope (Izuddin 2018)
  psi0 <- switch(stage,
    S = psi0S, J = psi0J,
    stop("growth_prob: stage must be 'S' or 'J'"))
  return(transition_logit(psi0, precip_annual_mm, relhum_mean, beta_precip, beta_rh))
}

# grow() archived to archive.R — superseded by the vectorized stage-
# transition + adult-size update inlined in run_pass3_survive_grow().
# growth_prob() above is still live and called directly there.

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
# z is the tracked mean pseudobulb size (cm) for the cell, updated annually
# by the adult size-update step in run_pass3_survive_grow().
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

# establish() archived to archive.R — superseded by the vectorized
# per-height establishment draw inlined in run_pass2_establish(), which also
# corrects the double-application of p_germ this archived version still has
# (see run_pass2_establish() below for the current, corrected design).

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

# ricker()/stochRicker() archived to archive.R — never wired into the main
# loop (Ricker model for within-voxel density regulation, kept for future use).

# ── Climate helpers ───────────────────────────────────────────────────────────

# List of usable height tiers for a microenv object, regardless of storage
# format (see load_height() below for the two formats).
microenv_heights <- function(microenv) {
  if (!is.null(microenv$.height_dir)) return(sort(as.numeric(microenv$.heights)))
  h_keys <- names(microenv)[!names(microenv) %in% c(".spatial", ".weather")]
  sort(as.numeric(sub("h", "", h_keys)))
}

# Load one height tier's raw tmax/tmin rasters.
# New-format microenv objects (run_microclimate_site.R) only carry a manifest
# — heights (.heights) and a directory (.height_dir) of per-height RDS files
# on scratch, saved this way because loading every height into memory at once
# needs hundreds of GB. Old-format microenv objects still embed each height's
# data directly as a list element named "h<value>".
load_height <- function(microenv, height) {
  if (!is.null(microenv$.height_dir)) {
    avail  <- microenv$.heights
    h_near <- avail[which.min(abs(avail - height))]
    h_path <- file.path(microenv$.height_dir, sprintf("h%.2f.rds", h_near))
    if (!file.exists(h_path)) return(NULL)
    return(readRDS(h_path))
  }
  h_keys <- names(microenv)[!names(microenv) %in% c(".spatial", ".weather")]
  avail  <- as.numeric(sub("h", "", h_keys))
  h_key  <- h_keys[which.min(abs(avail - height))]
  microenv[[h_key]]
}

# Build the climate data frame for one height tier.
# Spatial mean across the raster at that height → one value per representative
# timestep. microclimf::subsetpointmodela(tstep="month") + runmicro() produce
# one representative day per calendar month, 24 hourly values each, so a full
# year's tmax (or tmin) array holds 12 x 24 = 288 hourly steps concatenated in
# month order (Jan..Dec). Each row is tagged with its calendar `month` so
# get_clim_month() can select the right block regardless of the exact hours-
# per-month count. Older microenv objects with just one representative day
# (24 steps, no month structure) fall back to reusing that day for every month.
# Macro variables (precip, winddir) come from microenv$.weather (ERA5 hourly),
# replicated across every row (annual mean, not month-resolved).
get_clim <- function(height, microenv) {
  h <- load_height(microenv, height)
  if (is.null(h)) return(NULL)

  .smean <- function(arr) {
    if (length(dim(arr)) == 3) apply(arr, 3, mean, na.rm = TRUE)
    else rep(mean(arr, na.rm = TRUE), 24)
  }

  make_df <- function(slot, day_type) {
    temp      <- .smean(slot$Tz)
    relhum    <- .smean(slot$relhum)
    windspeed <- .smean(slot$windspeed)
    swdown    <- .smean(slot$Rdirdown) + .smean(slot$Rdifdown)
    difrad    <- .smean(slot$Rdifdown)
    n <- length(temp)
    if (n %% 12 == 0 && n > 24) {
      month <- rep(1:12, each = n / 12)
    } else {
      # single representative day (old format) — reuse it for every month
      month     <- rep(1:12, each = n)
      temp      <- rep(temp, 12);      relhum <- rep(relhum, 12)
      windspeed <- rep(windspeed, 12); swdown <- rep(swdown, 12)
      difrad    <- rep(difrad, 12)
    }
    data.frame(day_type = day_type, month = month, temp = temp, relhum = relhum,
               windspeed = windspeed, swdown = swdown, difrad = difrad)
  }

  df <- rbind(make_df(h$tmax, "tmax"), make_df(h$tmin, "tmin"))

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

# Return the rows of clim (both tmax and tmin blocks) for a given calendar month.
get_clim_month <- function(clim, month) {
  if (is.null(clim) || nrow(clim) == 0) return(NULL)
  clim[clim$month == month, ]
}

# Build the full per-height (and per-height-per-month) climate lookup table
# for one microenv object. This is the same for every run against a given
# site/microenv, no matter what biological parameters are being tested — so
# callers that run many simulations against the same microenv (e.g.
# run_experiment()'s parameter sweep) should build this once up front and
# pass it into runcolonization()/init_colonization() via `clim_cache` instead
# of letting each run reload every height's climate raster from disk.
build_clim_cache <- function(microenv) {
  heights              <- microenv_heights(microenv)
  clim_by_height        <- lapply(heights, function(h) get_clim(h, microenv))
  clim_month_by_height <- lapply(clim_by_height, function(clim)
    lapply(1:12, function(m) get_clim_month(clim, m)))
  list(clim_by_height = clim_by_height, clim_month_by_height = clim_month_by_height)
}

# ── Species climate niche ──────────────────────────────────────────────────────
# Realized niche per species: the range of each climate variable at the
# heights where that species was actually observed, centred on the observed
# midpoint and guaranteed at least `min_width_frac` of the SITE's full
# vertical climate range wide — even at pad = 0. Earlier version padded
# proportionally to the raw observed width with a small absolute floor
# (0.5 C / 5% RH); that failed silently because pad=0 means 0 x anything = 0
# regardless of the floor, and observed ranges here are ~0.1 C wide (2-3
# observations per species, often at similar heights) — so *no* pad level
# ever produced a niche wider than the observed range, gating establishment
# to a handful of height tiers out of 55 regardless of every other
# parameter. A floor expressed as a fraction of the site's own vertical
# climate range is a physically meaningful minimum tolerance instead of an
# arbitrary constant, and applies unconditionally so pad=0 still yields a
# usable niche. pad > 0 widens further beyond that floor.
#
# Takes the already-computed heights/clim_by_height (from init_colonization,
# via build_clim_cache) rather than re-deriving climate from microenv, since
# re-reading every height tier's raster from disk here would otherwise
# duplicate the same expensive I/O.
#
# The pad-independent geometry (mid, base_half_width) and the pad-application
# step are factored out into niche_geometry()/apply_niche_pad() below so the
# same math can be reused by characterize_niches.R, which pools observations
# across *all* sites for a species instead of just the one being modeled —
# see scripts/characterize_niches.R and init_colonization()'s cache lookup.

# Pad-independent niche geometry from a matrix of observed climate values
# (rows = observations, cols = variables) and a per-variable minimum width.
# Centred on the observed midpoint; base_half_width is guaranteed at least
# min_width/2 regardless of how narrow the raw observed range is.
niche_geometry <- function(clim_vals, min_width) {
  lo  <- apply(clim_vals, 2, min, na.rm = TRUE)
  hi  <- apply(clim_vals, 2, max, na.rm = TRUE)
  mid <- (lo + hi) / 2
  base_half_width <- pmax((hi - lo) / 2, min_width / 2)
  list(mid = mid, base_half_width = base_half_width)
}

# Apply a tolerance pad to a niche_geometry() result, producing the final
# lo/hi box used by niche_match(). pad=0 keeps the geometry's guaranteed
# minimum width; pad>0 widens further beyond it.
apply_niche_pad <- function(geom, pad) {
  half_width <- geom$base_half_width * (1 + pad)
  list(lo = geom$mid - half_width, hi = geom$mid + half_width)
}

get_niche <- function(site_obs, heights, clim_by_height, vars = c("temp", "relhum"),
                      pad = 0, min_width_frac = 0.10) {
  clim_scalars <- do.call(rbind, lapply(clim_by_height, function(cl) {
    if (is.null(cl)) return(setNames(rep(NA_real_, length(vars)), vars))
    vapply(vars, function(v) mean(cl[[v]], na.rm = TRUE), numeric(1))
  }))
  site_range <- apply(clim_scalars, 2, function(x) diff(range(x, na.rm = TRUE)))
  min_width  <- min_width_frac * site_range

  obs_clim <- function(h) clim_scalars[which.min(abs(heights - h)), , drop = TRUE]

  species_ids <- sort(unique(site_obs$FinalID))
  niches <- lapply(species_ids, function(sp) {
    obs_sp <- site_obs[site_obs$FinalID == sp & !is.na(site_obs$Height_m), ]
    if (nrow(obs_sp) == 0) return(NULL)
    clim_vals <- do.call(rbind, lapply(obs_sp$Height_m, obs_clim))
    apply_niche_pad(niche_geometry(clim_vals, min_width), pad)
  })
  names(niches) <- species_ids
  niches
}

# Fraction of niche variables whose value at a given height/time falls within
# a species' (padded) niche box — 1 = matches on every variable, 0 = outside
# on all of them. A species with no usable observations (niche = NULL) isn't
# gated at all (returns 1), since we have no basis to restrict it.
niche_match <- function(clim_values, niche) {
  if (is.null(niche)) return(1)
  vars <- names(niche$lo)
  hits <- vapply(vars, function(v) {
    val <- clim_values[[v]]
    if (is.null(val) || is.na(val)) return(NA)
    val >= niche$lo[[v]] && val <= niche$hi[[v]]
  }, logical(1))
  if (all(is.na(hits))) return(1)
  mean(hits, na.rm = TRUE)
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

          # Trunk zones: only the single column cell; crown zones (3-5, starting
          # at rel_h=0.30): bell-shaped radius. Domain must start at zone 3's
          # actual lower boundary (0.30), not 0.5 — using 0.5 here previously
          # made crown_fraction negative (and effective_r therefore negative,
          # i.e. never satisfied by any horiz_dist >= 0) for the entire
          # 0.30-0.50 sub-range, silently excluding that whole band from valid
          # canopy habitat on every tree.
          in_tree <- if (jzone <= 2) {
            horiz_dist == 0
          } else {
            crown_fraction <- (rel_h - 0.30) / 0.70     # 0 at zone 3 base, 1 at top
            effective_r    <- cr * sin(crown_fraction * pi)  # peaks at 65% height
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
            crown_fraction <- (rel_h - 0.30) / 0.70
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
                              forestparams = NULL, allsites = FALSE,
                              clim_cache = NULL) {
  site_name <- site$Site
  heights   <- microenv_heights(microenv)
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

  if (is.null(clim_cache)) {
    log_msg("Pre-computing climate lookup table from microenv...")
    clim_cache <- build_clim_cache(microenv)
    log_msg("Climate lookup ready.")
  } else {
    log_msg("Using precomputed climate lookup table.")
  }
  clim_by_height       <- clim_cache$clim_by_height
  clim_month_by_height <- clim_cache$clim_month_by_height

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

  # Per-species realized climate niche, gating establishment in
  # run_pass2_establish() — see get_niche()/niche_match(). niche_pad = 0
  # still yields a usable (site-relative minimum-width) niche, not a
  # zero-width point; make_params.R's niche_pad experiment sweeps pad to
  # widen further beyond that floor.
  #
  # Prefer the pooled cross-site cache from characterize_niches.R (every
  # observation of a species across all 5 sites, not just this one) if it
  # exists and covers every species observed at this site — a species with
  # only 2-3 records at one site may have several more elsewhere. Falls back
  # to this-site-only get_niche() for any species the cache doesn't cover
  # (e.g. added to the field data since the cache was last regenerated), or
  # entirely if the cache doesn't exist at all.
  niche_pad     <- if (!is.null(params$niche_pad)) params$niche_pad else 0
  niche_cache_path <- file.path(PROCESSED_DIR, "species_niches.rds")
  species_ids   <- sort(unique(site_obs$FinalID))
  niches_by_species <- setNames(vector("list", length(species_ids)), species_ids)

  cached_geom <- if (file.exists(niche_cache_path)) readRDS(niche_cache_path) else NULL
  missing_from_cache <- character(0)
  for (sp in species_ids) {
    if (!is.null(cached_geom) && !is.null(cached_geom[[sp]])) {
      niches_by_species[[sp]] <- apply_niche_pad(cached_geom[[sp]], niche_pad)
    } else {
      missing_from_cache <- c(missing_from_cache, sp)
    }
  }
  if (length(missing_from_cache) > 0) {
    fallback <- get_niche(site_obs[site_obs$FinalID %in% missing_from_cache, ],
                          heights, clim_by_height, pad = niche_pad)
    niches_by_species[missing_from_cache] <- fallback[missing_from_cache]
  }

  n_cached <- length(species_ids) - length(missing_from_cache)
  log_msg(sprintf(
    "Niche pad %.3f: %d/%d species have a usable observed niche (%d from cross-site cache, %d this-site-only)",
    niche_pad, sum(!sapply(niches_by_species, is.null)), n_species,
    n_cached, length(missing_from_cache)))

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
    niches_by_species    = niches_by_species,
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

      # p_establish: scalar for this height (spatial mean already in clim).
      # No p_germ term here — mycorrhizal germination potential is already
      # accounted for once, in Pass 1's fecundity kernel (reproduce()), so
      # every dispersed seed has already "passed" that gate. This gate is
      # purely site suitability (humidity/light) for a seed that already has
      # germination potential.
      relhum_mean <- mean(clim$relhum, na.rm = TRUE)
      swdown_mean <- mean(clim$swdown[clim$swdown > 0], na.rm = TRUE)
      p_est <- (relhum_mean / 100) *
               min(swdown_mean / p$mean_swdown_site, 1)

      # Gate by this species' realized climate niche (see get_niche() /
      # niche_match() in init_colonization()) — a height whose climate falls
      # outside where the species was actually observed is less likely to be
      # colonized, tempered by params$niche_pad.
      niche_sp <- state$niches_by_species[[state$species_ids[sp]]]
      p_est <- p_est * niche_match(
        list(temp = mean(clim$temp, na.rm = TRUE), relhum = relhum_mean),
        niche_sp
      )

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

    # ── Monthly loop: survival, transitions, and growth are all evaluated per
    # month and accumulated across the year, interleaved (survive -> maybe
    # transition -> accrue growth, each month in sequence) rather than
    # survival looping monthly while transitions/growth used one annual
    # value — an individual must survive a given month to be eligible for
    # anything else that month.
    #
    # Annual precipitation doesn't vary by month in this dataset (site-wide
    # ERA5 constant, replicated across every row regardless of month tag) —
    # only temperature and RH genuinely vary month to month. Growth's
    # monthly increment therefore still uses the annual precip ratio.
    #
    # All three vital rates are calibrated (Table params_stable/varied) for
    # what a SINGLE evaluation should reproduce ANNUALLY, so each is
    # converted to its monthly-equivalent rate before being applied 12x:
    #   survival: p_month = p_annual^(1/12) — must succeed every month
    #     (multiplicative), so a straight 12th root recovers the annual rate.
    #   transitions: q_month = 1-(1-p_annual)^(1/12) — "at least once this
    #     year" event, so the complement (not-yet-transitioned probability)
    #     is what compounds multiplicatively across months.
    #   growth: delta_z_base/12 per month (additive, not compounded) — see
    #     size_A update below.
    # The corresponding intercepts (beta0S/J/A, psi0S/J) already encode this
    # conversion; see methods.tex for the derivation.
    for (zi in seq_len(state$zDim)) {
      if (sum(abundanceS[,,zi,t,sp]) + sum(abundanceJ[,,zi,t,sp]) +
          sum(abundanceA[,,zi,t,sp]) == 0L) next

      clim_year <- state$clim_by_height[[zi]]
      if (is.null(clim_year)) next
      precip_annual <- mean(clim_year$precip, na.rm = TRUE) * 8760

      slS <- abundanceS[,,zi,t,sp]
      slJ <- abundanceJ[,,zi,t,sp]
      slA <- abundanceA[,,zi,t,sp]
      slA_start     <- slA
      n_JtoA_year   <- array(0L, dim = c(xDim, yDim))
      za_slice      <- size_A[,,zi,t,sp]
      delta_z_total <- array(0, dim = c(xDim, yDim))

      for (month in 1:12) {
        cm <- state$clim_month_by_height[[zi]][[month]]
        if (is.null(cm) || nrow(cm) == 0) next
        temp        <- mean(cm$temp,   na.rm = TRUE)
        relhum      <- mean(cm$relhum, na.rm = TRUE)
        swdown_mean <- mean(cm$swdown[cm$swdown > 0], na.rm = TRUE)
        swdown_rel  <- if (!is.na(swdown_mean) && p$mean_swdown_site > 0)
                         swdown_mean / p$mean_swdown_site else 1.0

        # ── s(z,e): monthly survival ────────────────────────────────────────
        z_S <- runif(1, p$z_S_min, p$z_S_max)
        z_J <- runif(1, p$z_J_min, p$z_J_max)
        z_A_mean <- mean(za_slice[state$landscape[,,zi]], na.rm = TRUE)
        if (is.na(z_A_mean) || z_A_mean < p$z_A_min) z_A_mean <- p$z_A_min

        s_S <- survival_logit("S", z_S,      temp, relhum, swdown_rel, p$beta0S, p$beta0J, p$beta0A, p$beta1)
        s_J <- survival_logit("J", z_J,      temp, relhum, swdown_rel, p$beta0S, p$beta0J, p$beta0A, p$beta1)
        s_A <- survival_logit("A", z_A_mean, temp, relhum, swdown_rel, p$beta0S, p$beta0J, p$beta0A, p$beta1)

        slS <- .surv_slice(slS, s_S)
        slJ <- .surv_slice(slJ, s_J)
        slA <- .surv_slice(slA, s_A)

        # ── g(z'|z,e): monthly stage transitions ────────────────────────────
        p_StoJ <- growth_prob("S", precip_annual, relhum, p$psi0S, p$psi0J, p$beta_precip, p$beta_rh)
        p_JtoA <- growth_prob("J", precip_annual, relhum, p$psi0S, p$psi0J, p$beta_precip, p$beta_rh)
        epsilon <- rnorm(1, 0, p$sigma)
        p_StoJ  <- pmin(1, pmax(0, p_StoJ + epsilon))
        p_JtoA  <- pmin(1, pmax(0, p_JtoA + epsilon))

        n_StoJ <- .surv_slice(slS, p_StoJ)
        n_JtoA <- .surv_slice(slJ, p_JtoA)
        slS <- slS - n_StoJ
        slJ <- slJ + n_StoJ - n_JtoA
        slA <- slA + n_JtoA
        n_JtoA_year <- n_JtoA_year + n_JtoA

        # ── Adult size increment: this month's share, own noise draw ────────
        # sigma/sqrt(12) keeps total annual variance matched to the original
        # (pre-monthly) calibration, since variances of independent draws sum.
        dz_month <- (p$delta_z_base / 12) *
                    (precip_annual / 2500) *
                    (relhum        / 85) +
                    rnorm(xDim * yDim, 0, p$sigma * 0.5 / sqrt(12))
        dim(dz_month) <- c(xDim, yDim)
        delta_z_total <- delta_z_total + dz_month
      }

      abundanceS[,,zi,t,sp] <- slS
      abundanceJ[,,zi,t,sp] <- slJ
      abundanceA[,,zi,t,sp] <- slA

      # Cost of reproduction: reduce the year's total growth where the cell fruited
      fr_slice <- fruited[,,zi,sp]
      delta_z_total[fr_slice] <- delta_z_total[fr_slice] * p$cost_repro
      za_new <- pmin(p$z_A_max, pmax(p$z_A_min, za_slice + delta_z_total))
      # Cells newly promoted from J this year (no adults at year start) start at z_A_min
      fresh <- n_JtoA_year > 0L & slA_start == 0L
      za_new[fresh] <- p$z_A_min
      # Only write back where adults exist (preserve prior value in empty cells)
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

# Same idea as plot_3d_abundance() but across every timestep, using plotly's
# built-in frame/animation support (play button + slider) instead of a
# single static scatter. Saved as a self-contained HTML if out_path is
# given. Not yet run against real output — verify once you have a result
# worth animating (e.g. from a best_case replicate that actually persists).
plot_3d_abundance_animated <- function(result, out_path = NULL) {
  state <- result$state
  n_t <- dim(result$abundanceA)[4]
  sp_cols <- if (state$n_species == 1)
    scico::scico(3, palette = "lipari", begin = 0.3, end = 0.7)[2]
  else
    scico::scico(state$n_species, palette = "lipari", begin = 0.2, end = 0.8)

  rows <- list()
  for (t in seq_len(n_t)) {
    for (sp in seq_len(state$n_species)) {
      sp_name <- state$species_ids[sp]; col <- sp_cols[sp]
      for (stg in list(list(arr = result$abundanceS, nm = "S", sym = "circle"),
                       list(arr = result$abundanceJ, nm = "J", sym = "diamond"),
                       list(arr = result$abundanceA, nm = "A", sym = "square"))) {
        idx <- which(stg$arr[, , , t, sp] > 0, arr.ind = TRUE)
        if (nrow(idx) > 0)
          rows[[length(rows) + 1]] <- data.frame(
            x = idx[, 1], y = idx[, 2], z = idx[, 3], t = t,
            species = sp_name, stage = stg$nm, color = col,
            symbol = stg$sym, trace = paste0(sp_name, " ", stg$nm),
            stringsAsFactors = FALSE)
      }
    }
  }
  if (length(rows) == 0) {
    message("No individuals to plot across any timestep")
    return(invisible(NULL))
  }
  df <- do.call(rbind, rows)

  # trace -> color lookup (one row per unique trace, in matching order —
  # safer than pairing two independently-deduplicated vectors)
  trace_lu     <- df[!duplicated(df$trace), c("trace", "color")]
  colors_named <- setNames(trace_lu$color, trace_lu$trace)

  fig <- plotly::plot_ly(
    df, x = ~x, y = ~y, z = ~z, frame = ~t, color = ~trace,
    colors = colors_named,
    symbol = ~symbol, symbols = c(circle = "circle", diamond = "diamond", square = "square"),
    type = "scatter3d", mode = "markers",
    marker = list(size = 6, opacity = 0.85)
  ) |>
    plotly::layout(
      title = "Abundance over time",
      scene = list(xaxis = list(title = "x"), yaxis = list(title = "y"),
                   zaxis = list(title = "height tier"))
    ) |>
    plotly::animation_opts(frame = 400, transition = 200, redraw = TRUE) |>
    plotly::animation_slider(currentvalue = list(prefix = "Year: "))

  if (!is.null(out_path)) {
    # selfcontained=TRUE needs pandoc (not installed on the cluster); FALSE
    # writes a small "<name>_files/" dependency folder alongside the HTML
    # instead -- keep the two together when copying/viewing elsewhere.
    htmlwidgets::saveWidget(fig, out_path, selfcontained = FALSE)
    message("Saved: ", out_path)
  }
  invisible(fig)
}

# ── Spin-up ───────────────────────────────────────────────────────────────────

# Founders: n_founders individuals per species (a fixed, decoupled count —
# not tied to however many field observations that species happens to have),
# placed at random canopy voxels whose local climate matches that species'
# realized niche (see get_niche()/niche_match()). Previously founders were
# placed only at the exact voxel of an observed individual, which failed
# whenever the stochastic forest didn't happen to mark that exact voxel as
# canopy — with a handful of observations per site, that could (and did)
# place zero founders across every replicate, guaranteeing extinction before
# the simulation even started, independent of any vital-rate parameter.
# Falls back to unrestricted canopy placement if no voxel matches the niche
# (species with too few observations to build one, or a niche too narrow for
# this stochastic forest) so founder placement never silently fails.
.niche_matching_canopy <- function(state, niche_sp) {
  zDim <- state$zDim
  height_ok <- vapply(seq_len(zDim), function(zi) {
    clim <- state$clim_by_height[[zi]]
    if (is.null(clim)) return(FALSE)
    niche_match(list(temp   = mean(clim$temp,   na.rm = TRUE),
                     relhum = mean(clim$relhum, na.rm = TRUE),
                     precip = mean(clim$precip, na.rm = TRUE)),
                niche_sp) >= 1
  }, logical(1))

  candidate_list <- lapply(which(height_ok), function(zi) {
    xy <- which(state$landscape[, , zi], arr.ind = TRUE)
    if (nrow(xy) == 0) return(NULL)
    cbind(xy, z = zi)
  })
  do.call(rbind, Filter(Negate(is.null), candidate_list))
}

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

  n_founders <- if (!is.null(p$n_founders)) p$n_founders else 30

  for (sp in seq_len(n_species)) {
    sp_name    <- state$species_ids[sp]
    niche_sp   <- state$niches_by_species[[sp_name]]
    candidates <- .niche_matching_canopy(state, niche_sp)
    if (is.null(candidates) || nrow(candidates) == 0) {
      log_msg(sprintf(
        "Spin-up: no niche-matching canopy for %s, falling back to unrestricted placement", sp_name))
      candidates <- which(state$landscape, arr.ind = TRUE)
    }
    if (nrow(candidates) == 0) {
      log_msg(sprintf("Spin-up: no canopy at all for %s -- 0 founders placed", sp_name))
      next
    }

    chosen <- candidates[sample.int(nrow(candidates), n_founders,
                                    replace = nrow(candidates) < n_founders), , drop = FALSE]
    for (k in seq_len(nrow(chosen))) {
      x <- chosen[k, 1]; y <- chosen[k, 2]; z <- chosen[k, 3]
      spinupA[x, y, z, 1, sp] <- spinupA[x, y, z, 1, sp] + 1L
      size_A[x, y, z, 1, sp]  <- p$z_A_min
    }
  }
  log_msg(sprintf("Spin-up: placed %d founders/species x %d species = %d total",
                  n_founders, n_species, n_founders * n_species))

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
                            train_frac=0.70, seed=42, clim_cache=NULL) {
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
                                 forestparams=forestparams, allsites=allsites,
                                 clim_cache=clim_cache)
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

run_one <- function(params, tag = "run", timesteps = 20, spinup = 3, clim_cache = NULL, seed = 42) {
  tryCatch(
    runcolonization(
      site         = site,
      niches       = niches,
      canopy_grid  = canopy_grid,
      microenv     = microenv,
      timesteps    = timesteps,
      resolution   = 10,
      carCap       = 5,
      maxDisp      = 10,
      spinup       = spinup,
      Visualize    = FALSE,
      parameters   = params,
      forestparams = forestparams,
      clim_cache   = clim_cache,
      seed         = seed
    ),
    error = function(e) {
      # run_experiment() wraps this call in suppressMessages(), so a plain
      # message() here would vanish with no trace anywhere. Write directly
      # to the shared log file (defined by the caller, e.g.
      # run_colonization_onesite.R) so failed workers are still visible.
      err_line <- sprintf("ERROR [%s]: %s", tag, e$message)
      message(err_line)
      if (exists("log_file", inherits = TRUE)) {
        cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", err_line, "\n"),
            file = log_file, append = TRUE)
      }
      NULL
    }
  )
}

# ── Experiment runner ─────────────────────────────────────────────────────────
# Varies param_name across values; n_reps replicates per value.
# Returns tidy data frame of S, J, A totals over time.

run_experiment <- function(param_name, values, base = base_params,
                           n_reps = 2, timesteps = 20, spinup = 3) {
  jobs <- expand.grid(val = values, rep = seq_len(n_reps),
                      stringsAsFactors = FALSE)
  cat(sprintf("\n====== %s (%d jobs on %d cores) ======\n",
              param_name, nrow(jobs), N_CORES))

  # Climate is identical across every value/rep in this sweep (only the
  # biological parameter differs) — build the per-height lookup table once
  # here, before mclapply forks, so all workers inherit it via copy-on-write
  # instead of every one of them re-reading every height's raster from disk.
  cat("Pre-computing shared climate lookup table for the sweep...\n")
  clim_cache <- build_clim_cache(microenv)

  rows <- mclapply(seq_len(nrow(jobs)), function(i) {
    val <- jobs$val[i]; rep <- jobs$rep[i]
    p   <- base; p[[param_name]] <- val
    # Suppress log_msg inside workers
    suppressMessages(
      r <- run_one(p, tag = sprintf("%s=%s rep%d", param_name, val, rep),
                   timesteps = timesteps, spinup = spinup, clim_cache = clim_cache,
                   seed = i)
    )
    if (is.null(r)) return(NULL)
    T <- timesteps
    data.frame(
      param_value = as.character(val), rep = rep, t = 1:T,
      totalS = r$totalabundanceS, totalJ = r$totalabundanceJ,
      totalA = r$totalabundanceA,
      total  = r$totalabundanceS + r$totalabundanceJ + r$totalabundanceA,
      extinct = all(r$totalabundanceA[(T %/% 2):T] == 0)
    )
  }, mc.cores = N_CORES)

  df <- do.call(rbind, Filter(Negate(is.null), rows))
  df$param_value <- factor(df$param_value, levels = as.character(values))
  cat(sprintf("  Done: %d/%d runs succeeded\n",
              length(Filter(Negate(is.null), rows)), nrow(jobs)))
  df
}

# ── Factorial experiment runner ────────────────────────────────────────────────
# Crosses several parameters at once (unlike run_experiment(), which varies
# only one). `param_values` is a named list, e.g.
#   list(p_poll = c(...), p_germ = c(...), p_s1 = c(...))
# giving length(p_poll) x length(p_germ) x length(p_s1) x n_reps jobs.
# Returns a tidy data frame of S/J/A totals over time, with one column per
# swept parameter recording the value used in that run.
#
# If checkpoint_var/checkpoint_path are given, jobs are processed in blocks
# — one block per unique value of checkpoint_var — and the accumulated
# results so far are saveRDS()'d to checkpoint_path after every block. This
# bounds how much work is lost if the process is killed partway through a
# long factorial (e.g. a SLURM walltime limit): at most one block's worth,
# instead of the entire sweep, since without this the only saveRDS() call
# happens after ALL jobs finish. mclapply's built-in mc.cores parallelism is
# unaffected within each block; blocks just run sequentially relative to
# each other, so total wall-clock time is essentially unchanged.
run_factorial_experiment <- function(param_values, base = base_params,
                                     n_reps = 1, timesteps = 20, spinup = 3,
                                     checkpoint_var = NULL, checkpoint_path = NULL) {
  param_names <- names(param_values)
  jobs <- do.call(expand.grid,
                  c(param_values, list(rep = seq_len(n_reps)), stringsAsFactors = FALSE))
  cat(sprintf("\n====== factorial %s (%d combos x %d reps = %d jobs on %d cores) ======\n",
              paste(param_names, collapse = " x "),
              nrow(jobs) / n_reps, n_reps, nrow(jobs), N_CORES))

  # Climate doesn't depend on any of the swept parameters — build it once,
  # before mclapply forks, same reasoning as run_experiment().
  cat("Pre-computing shared climate lookup table for the factorial...\n")
  clim_cache <- build_clim_cache(microenv)

  run_job <- function(i) {
    p <- base
    for (nm in param_names) p[[nm]] <- jobs[[nm]][i]
    tag <- paste(sprintf("%s=%s", param_names,
                         sapply(param_names, function(nm) jobs[[nm]][i])),
                collapse = " ")
    tag <- paste0(tag, sprintf(" rep%d", jobs$rep[i]))
    suppressMessages(
      r <- run_one(p, tag = tag, timesteps = timesteps, spinup = spinup,
                  clim_cache = clim_cache, seed = i)
    )
    if (is.null(r)) return(NULL)
    T <- timesteps
    df <- data.frame(
      rep = jobs$rep[i], t = 1:T,
      totalS = r$totalabundanceS, totalJ = r$totalabundanceJ,
      totalA = r$totalabundanceA,
      total  = r$totalabundanceS + r$totalabundanceJ + r$totalabundanceA,
      extinct = all(r$totalabundanceA[(T %/% 2):T] == 0)
    )
    for (nm in param_names) df[[nm]] <- jobs[[nm]][i]
    df
  }

  use_checkpoints <- !is.null(checkpoint_var) && !is.null(checkpoint_path) &&
    checkpoint_var %in% param_names

  if (!use_checkpoints) {
    rows <- mclapply(seq_len(nrow(jobs)), run_job, mc.cores = N_CORES)
    df <- do.call(rbind, Filter(Negate(is.null), rows))
    cat(sprintf("  Done: %d/%d runs succeeded\n",
                length(Filter(Negate(is.null), rows)), nrow(jobs)))
    return(df)
  }

  blocks    <- split(seq_len(nrow(jobs)), jobs[[checkpoint_var]])
  all_rows  <- list()
  n_done    <- 0L
  for (block_val in names(blocks)) {
    idx  <- blocks[[block_val]]
    cat(sprintf("  -- block %s=%s: %d jobs --\n", checkpoint_var, block_val, length(idx)))
    rows <- Filter(Negate(is.null), mclapply(idx, run_job, mc.cores = N_CORES))
    n_done   <- n_done + length(rows)
    all_rows <- c(all_rows, rows)
    df_so_far <- do.call(rbind, all_rows)
    saveRDS(df_so_far, checkpoint_path)
    cat(sprintf("  Checkpoint saved (%d/%d jobs done so far): %s\n",
                n_done, nrow(jobs), checkpoint_path))
  }
  cat(sprintf("  Done: %d/%d runs succeeded\n", n_done, nrow(jobs)))
  do.call(rbind, all_rows)
}

# ── Replicated runner for a single (non-swept) parameter set ───────────────────
# Runs n_reps independent replicates of one fixed params set in parallel —
# e.g. to check whether an outcome (like establishment never succeeding) is
# genuinely blocked or just one unlucky stochastic draw. Unlike
# run_experiment()/run_factorial_experiment() (which discard each run's full
# spatial arrays down to S/J/A totals, since they cover hundreds of combos),
# this keeps every replicate's complete runcolonization() output — reasonable
# since n_reps is typically a handful, not hundreds — so results stay usable
# for spatial/animation plotting (see plot_3d_abundance_animated()), not just
# aggregate totals. Returns list(runs = <one runcolonization() output per
# replicate>, summary = <tidy S/J/A-over-time data frame, like
# run_experiment()'s output but without a param_value column>).
run_replicated <- function(params, n_reps = 1, timesteps = 20, spinup = 3,
                           clim_cache = NULL) {
  cat(sprintf("\n====== %d replicate(s) on %d cores ======\n", n_reps, N_CORES))
  if (is.null(clim_cache)) clim_cache <- build_clim_cache(microenv)

  runs <- mclapply(seq_len(n_reps), function(rep) {
    suppressMessages(
      r <- run_one(params, tag = sprintf("rep%d", rep), timesteps = timesteps,
                  spinup = spinup, clim_cache = clim_cache, seed = rep)
    )
    r
  }, mc.cores = min(N_CORES, n_reps))

  ok <- !vapply(runs, is.null, logical(1))
  cat(sprintf("  Done: %d/%d replicates succeeded\n", sum(ok), n_reps))

  summary_df <- do.call(rbind, lapply(seq_along(runs), function(i) {
    r <- runs[[i]]
    if (is.null(r)) return(NULL)
    T <- timesteps
    data.frame(
      rep = i, t = 1:T,
      totalS = r$totalabundanceS, totalJ = r$totalabundanceJ,
      totalA = r$totalabundanceA,
      total  = r$totalabundanceS + r$totalabundanceJ + r$totalabundanceA,
      extinct = all(r$totalabundanceA[(T %/% 2):T] == 0)
    )
  }))

  list(runs = runs, summary = summary_df)
}

# ── Plot: S, J, A panels ──────────────────────────────────────────────────────

plot_experiment <- function(df, param_name, title = NULL) {
  title <- title %||% sprintf("Effect of %s — Maquipucuna", param_name)

  make_panel <- function(y_var, y_lab, col) {
    mean_df <- aggregate(as.formula(paste(y_var, "~ param_value + t")),
                         data = df, FUN = mean)
    ggplot(df, aes(x = t, y = .data[[y_var]], colour = param_value,
                   group = interaction(param_value, rep))) +
      geom_line(alpha = 0.20, linewidth = 0.4) +
      geom_line(data = mean_df,
                aes(x = t, y = .data[[y_var]], colour = param_value,
                    group = param_value),
                linewidth = 1.2, inherit.aes = FALSE) +
      scale_colour_brewer(palette = "RdYlBu", direction = -1, name = param_name) +
      labs(x = "Year", y = y_lab) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "right")
  }

  (make_panel("totalS", "Seedlings", "#4dac26") +
   make_panel("totalJ", "Juveniles", "#f1a340") +
   make_panel("totalA", "Adults",    "#08519c")) +
    plot_annotation(
      title   = title,
      caption = "Thick = mean of replicates; thin = individual runs"
    )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b