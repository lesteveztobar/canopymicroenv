# get_colonization.R
# Population dynamics simulation functions for canopy colonization
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────

# ── Demographic helpers ───────────────────────────────────────────────────────

ricker <- function(N, r, K) {
  N * exp(r * (1 - (N / K)))
}

stochRicker <- function(N, r, K) {
  rpois(1, lambda = ricker(N, r, K))
}

indDisp <- function(x, y, z, indN, winddir, meanDisp, Disp, pad, maxDispZ = 5) {
  wind_rad <- (winddir + 180) %% 360 * pi / 180
  for (ind in 1:indN) {
    dist  <- round(min(rexp(1, rate = 1 / max(meanDisp, 0.1)), pad))
    angle <- wind_rad + runif(1, -pi / 4, pi / 4)
    dx    <- round(dist * sin(angle))
    dy    <- round(dist * cos(angle))
    dz    <- sample(-maxDispZ:maxDispZ, size = 1)
    tx <- x + dx + pad
    ty <- y + dy + pad
    tz <- z + dz + pad
    if (tx >= 1 && tx <= dim(Disp)[1] &&
        ty >= 1 && ty <= dim(Disp)[2] &&
        tz >= 1 && tz <= dim(Disp)[3]) {
      Disp[tx, ty, tz] <- Disp[tx, ty, tz] + 1L
    }
  }
  return(Disp)
}

reproduce <- function(N, repRate, stochastic = FALSE) {
  seeds <- N * repRate
  if (stochastic) rpois(1, lambda = max(seeds, 0))
  else floor(seeds) + rbinom(1, 1, seeds - floor(seeds))
}

disperse <- function(x, y, z, seeds, clim, height, canopy_z, a,
                     lambda = 1, Ut = 1, maxDisp = 5, maxDispZ = 5,
                     Disp, stochastic = FALSE) {
  wind     <- mean(clim$windspeed, na.rm = TRUE)
  winddir  <- mean(clim$winddir,   na.rm = TRUE)
  meanDisp <- max(1, min(round((wind * exp(-(a * height / canopy_z))) / (lambda * Ut)), maxDisp))
  indDisp(x = x, y = y, z = z, indN = seeds, winddir = winddir,
          meanDisp = meanDisp, Disp = Disp, pad = maxDisp, maxDispZ = maxDispZ)
}

establish <- function(new_seeds, clim, mean_swdown_site,
                      p_germ = 0.001, stochastic = FALSE) {
  relhum_mean <- mean(clim$relhum, na.rm = TRUE)
  swdown_mean <- mean(clim$swdown[clim$swdown > 0], na.rm = TRUE)
  p_establish <- (relhum_mean / 100) * min(swdown_mean / mean_swdown_site, 1) * p_germ
  established <- if (stochastic) rbinom(1, new_seeds, p_establish) > 0
  else (new_seeds * p_establish) >= 0.1
  if (established) 1L else 0L
}

survive <- function(N, stage, clim_month,
                    beta0S = -1.5, beta0J = -0.5, beta0A = 1.0,
                    beta1  = 0.105,
                    beta0GrowthS = 0.22, beta0GrowthJ = 0.45, beta0GrowthA = 1.40,
                    stochastic = FALSE) {
  if (is.na(N)) return(0L)
  temp   <- mean(clim_month$temp,   na.rm = TRUE)
  relhum <- mean(clim_month$relhum, na.rm = TRUE)
  s <- if (stage == "S") {
    1 / (1 + exp(-(beta0S + beta1 * runif(1, 0, beta0GrowthS) - (100 - relhum) / 100)))
  } else if (stage == "J") {
    1 / (1 + exp(-(beta0J + beta1 * runif(1, beta0GrowthS, beta0GrowthJ))))
  } else {
    1 / (1 + exp(-(beta0A + beta1 * runif(1, beta0GrowthJ, beta0GrowthA) -
                     max(0, (temp - 23) * 0.05))))
  }
  if (stochastic) rbinom(1, N, s) else round(N * s)
}

grow <- function(nS, nJ, nA,
                 beta0GrowthS = 0.22, beta0GrowthJ = 0.45, beta0GrowthA = 1.40,
                 beta1Growth  = 0.8,  sigma = 0.5) {
  nS <- if (is.na(nS)) 0L else nS
  nJ <- if (is.na(nJ)) 0L else nJ
  nA <- if (is.na(nA)) 0L else nA
  epsilon <- rnorm(1, 0, sigma)
  if (nS > 0 &&
      beta0GrowthS + beta1Growth * runif(1, 0, beta0GrowthS) + epsilon >= beta0GrowthJ) {
    nJ <- nJ + nS; nS <- 0L
  }
  if (nJ > 0 &&
      beta0GrowthJ + beta1Growth * runif(1, beta0GrowthS, beta0GrowthJ) + epsilon >= beta0GrowthA) {
    nA <- nA + nJ; nJ <- 0L
  }
  list(nS = nS, nJ = nJ, nA = nA)
}

# ── Climate helpers ───────────────────────────────────────────────────────────

# Multi-site get_clim: tries each site name until valid data found for height
get_clim <- function(height, models, valid_per_model, site_names) {
  for (site in site_names) {
    h_key  <- sprintf("%s_h%.1f", site, height)
    cell_c <- valid_per_model[[h_key]][1]
    if (!is.null(cell_c) && !is.na(cell_c))
      return(models[[h_key]][[cell_c]]$weather)
  }
  return(NULL)
}

get_clim_month <- function(clim, month) {
  if (is.null(clim) || nrow(clim) == 0) return(NULL)
  hours_per_month <- floor(nrow(clim) / 12)
  month_hours     <- ((month - 1) * hours_per_month + 1):(month * hours_per_month)
  month_hours     <- month_hours[month_hours <= nrow(clim)]
  clim[month_hours, ]
}

# ── Simulation setup ──────────────────────────────────────────────────────────

init_colonization <- function(site, niches, canopy_grid, models, valid_per_model,
                              resolution = 10, carCap = 1, maxDisp = 5, params,
                              allsites = FALSE) {
  site_name  <- site$Site
  
  # all site names present in models (kept for multi-site climate fallback in get_clim)
  site_names <- unique(sub("_h.*", "", names(models)))
  
  if (allsites) {
    # heights: union across all sites — requires more memory but spans full model range
    heights  <- sort(unique(as.numeric(sub(".*_h", "", names(models)))))
    # observations: all sites — xDim/yDim cover the combined spatial extent
    site_obs <- niches
  } else {
    # heights: current site only — avoids inflating zDim with other sites' heights
    site_model_keys <- names(models)[grepl(paste0("^", site_name, "_h"), names(models))]
    heights  <- sort(unique(as.numeric(sub(".*_h", "", site_model_keys))))
    # observations: current site only — keeps xDim/yDim to this site's extent
    site_obs <- niches[niches$Area_or_Site == site_name, ]
  }
  zDim <- length(heights)
  lat_range_m <- (max(site_obs$lat) - min(site_obs$lat)) * 111000
  lon_range_m <- (max(site_obs$lon) - min(site_obs$lon)) * 111000 *
    cos(mean(site_obs$lat) * pi / 180)
  xDim <- max(round(lon_range_m / resolution), 10) + 4
  yDim <- max(round(lat_range_m / resolution), 10) + 4
  
  # species across all sites
  species_ids <- sort(unique(site_obs$FinalID))
  n_species   <- length(species_ids)
  sp_index    <- setNames(seq_along(species_ids), species_ids)
  
  log_msg(sprintf("Landscape: %d x %d x %d cells (%.0f x %.0f m, %.1f-%.1fm height) | %d species",
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
  
  # landscape mask
  landscape <- array(FALSE, dim = c(xDim, yDim, zDim))
  cg_xDim   <- nrow(canopy_grid)
  cg_yDim   <- ncol(canopy_grid)
  for (x in 1:xDim) for (y in 1:yDim) for (z in 1:zDim) {
    cx <- min(x, cg_xDim)
    cy <- min(y, cg_yDim)
    landscape[x, y, z] <- heights[z] <= canopy_grid[cx, cy] && heights[z] >= 0.5
  }
  
  # climate lookup — multi-site: first valid site per height
  log_msg("Pre-computing climate lookup table...")
  clim_by_height       <- lapply(heights, function(h)
    get_clim(h, models, valid_per_model, site_names))
  clim_month_by_height <- lapply(clim_by_height, function(clim)
    lapply(1:12, function(m) get_clim_month(clim, m)))
  log_msg("Climate lookup ready.")
  
  # canopy openness from mid-height climate
  valid_clim <- which(!sapply(clim_by_height, is.null))
  mid_zi     <- valid_clim[which.min(abs(heights[valid_clim] -
                                           max(min(heights), mean(canopy_grid, na.rm=TRUE)/2)))]
  if (length(mid_zi) == 0) mid_zi <- valid_clim[1]
  clim_mid <- clim_by_height[[mid_zi]]
  difrac   <- mean(clim_mid$difrad / (clim_mid$swdown + 0.001), na.rm = TRUE)
  a        <- 23 * (1 - difrac)
  
  med_zi           <- valid_clim[which.min(abs(heights[valid_clim] - median(heights[valid_clim])))]
  mean_swdown_site <- mean(clim_by_height[[med_zi]]$swdown, na.rm = TRUE)
  
  log_msg(sprintf("Valid climate heights: %d/%d | mean_swdown: %.1f",
                  length(valid_clim), zDim, mean_swdown_site))
  
  params$a                <- a
  params$mean_swdown_site <- mean_swdown_site
  params$repRate          <- params$p_flower * params$p_poll * params$p_germ * params$p_s1
  
  list(
    site_name            = site_name,
    site_names           = site_names,
    site_obs             = site_obs,
    heights              = heights,
    xDim = xDim, yDim = yDim, zDim = zDim,
    n_species            = n_species,
    species_ids          = species_ids,
    sp_index             = sp_index,
    landscape            = landscape,
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

run_pass1_disperse <- function(state, abundanceA, t) {
  p    <- state$params
  Disp <- array(0L, dim = c(
    state$xDim + 2 * state$maxDisp,
    state$yDim + 2 * state$maxDisp,
    state$zDim + 2 * state$maxDisp,
    state$n_species))
  total_seeds <- 0L
  
  for (sp in 1:state$n_species) {
    nonzero <- which(abundanceA[,,,t,sp] > 0, arr.ind = TRUE)
    for (i in seq_len(nrow(nonzero))) {
      x <- nonzero[i,1]; y <- nonzero[i,2]; z <- nonzero[i,3]
      if (!state$landscape[x, y, z]) next
      N <- abundanceA[x, y, z, t, sp]
      if (is.na(N) || N == 0) next
      seeds <- reproduce(N, p$repRate)
      total_seeds <- total_seeds + seeds
      if (seeds == 0) next
      Disp[,,,sp] <- disperse(
        x = x, y = y, z = z, seeds = seeds,
        clim    = state$clim_by_height[[z]],
        height  = state$heights[z], canopy_z = p$canopy_z,
        a = p$a, lambda = p$lambda, Ut = p$Ut,
        maxDisp = state$maxDisp, Disp = Disp[,,,sp])
    }
  }
  message("Pass1: total seeds produced=", total_seeds, " dispersed=", sum(Disp))
  Disp
}

run_pass2_establish <- function(state, abundanceS, dispersalmatrix, t) {
  p                 <- state$params
  pad               <- state$maxDisp
  total_seeds_seen  <- 0L
  total_established <- 0L
  
  for (sp in 1:state$n_species) {
    nonzero <- which(dispersalmatrix[,,,sp] > 0, arr.ind = TRUE)
    for (i in seq_len(nrow(nonzero))) {
      x <- nonzero[i,1] - pad
      y <- nonzero[i,2] - pad
      z <- nonzero[i,3] - pad
      if (x < 1 || x > state$xDim ||
          y < 1 || y > state$yDim ||
          z < 1 || z > state$zDim) next
      if (!state$landscape[x, y, z]) next
      if (abundanceS[x, y, z, t, sp] > 0) next
      new_seeds <- dispersalmatrix[x+pad, y+pad, z+pad, sp]
      if (new_seeds == 0) next
      clim <- state$clim_by_height[[z]]
      if (is.null(clim)) next
      total_seeds_seen  <- total_seeds_seen + new_seeds
      result <- establish(new_seeds, clim, p$mean_swdown_site, p$p_germ)
      if (t < dim(abundanceS)[4])
        abundanceS[x, y, z, t+1, sp] <- result
      total_established <- total_established + result
    }
  }
  message("Pass2: seeds seen=", total_seeds_seen, " established=", total_established)
  abundanceS
}

run_pass3_survive_grow <- function(state, abundanceS, abundanceJ, abundanceA, t) {
  p <- state$params
  
  for (sp in 1:state$n_species) {
    occupied <- which(
      abundanceS[,,,t,sp] > 0 |
        abundanceJ[,,,t,sp] > 0 |
        abundanceA[,,,t,sp] > 0,
      arr.ind = TRUE)
    if (nrow(occupied) == 0) next
    
    for (month in 1:12) {
      for (i in seq_len(nrow(occupied))) {
        x <- occupied[i,1]; y <- occupied[i,2]; z <- occupied[i,3]
        if (!state$landscape[x, y, z]) next
        cm <- state$clim_month_by_height[[z]][[month]]
        if (is.null(cm) || nrow(cm) == 0) next
        abundanceS[x,y,z,t,sp] <- survive(abundanceS[x,y,z,t,sp], "S", cm,
                                          beta0S=p$beta0Seedling, beta0J=p$beta0Juvenile, beta0A=p$beta0Adult,
                                          beta1=p$beta1, beta0GrowthS=p$beta0GrowthS,
                                          beta0GrowthJ=p$beta0GrowthJ, beta0GrowthA=p$beta0GrowthA)
        abundanceJ[x,y,z,t,sp] <- survive(abundanceJ[x,y,z,t,sp], "J", cm,
                                          beta0S=p$beta0Seedling, beta0J=p$beta0Juvenile, beta0A=p$beta0Adult,
                                          beta1=p$beta1, beta0GrowthS=p$beta0GrowthS,
                                          beta0GrowthJ=p$beta0GrowthJ, beta0GrowthA=p$beta0GrowthA)
        abundanceA[x,y,z,t,sp] <- survive(abundanceA[x,y,z,t,sp], "A", cm,
                                          beta0S=p$beta0Seedling, beta0J=p$beta0Juvenile, beta0A=p$beta0Adult,
                                          beta1=p$beta1, beta0GrowthS=p$beta0GrowthS,
                                          beta0GrowthJ=p$beta0GrowthJ, beta0GrowthA=p$beta0GrowthA)
      }
    }
    
    for (i in seq_len(nrow(occupied))) {
      x <- occupied[i,1]; y <- occupied[i,2]; z <- occupied[i,3]
      if (!state$landscape[x, y, z]) next
      grown <- grow(abundanceS[x,y,z,t,sp], abundanceJ[x,y,z,t,sp],
                    abundanceA[x,y,z,t,sp],
                    beta0GrowthS=p$beta0GrowthS, beta0GrowthJ=p$beta0GrowthJ,
                    beta0GrowthA=p$beta0GrowthA, beta1Growth=p$beta1Growth,
                    sigma=p$sigma)
      abundanceS[x,y,z,t,sp] <- grown$nS
      abundanceJ[x,y,z,t,sp] <- grown$nJ
      abundanceA[x,y,z,t,sp] <- grown$nA
    }
  }
  list(S = abundanceS, J = abundanceJ, A = abundanceA)
}

# ── Internal live visualisation (called during simulation) ────────────────────

.plot_live <- function(state, abundanceS, abundanceJ, abundanceA,
                       totalS, totalJ, totalA, t, carCap, sleeptime = 0.2) {
  stage_cols <- scico::scico(3, palette = "lipari", begin = 0.2, end = 0.8)
  sp_cols    <- if (state$n_species == 1)
    scico::scico(3, palette = "lipari", begin = 0.3, end = 0.7)[2]
  else
    scico::scico(state$n_species, palette = "lipari", begin = 0.2, end = 0.8)

  total <- totalS + totalJ + totalA
  n_sp  <- state$n_species
  par(mfrow = c(1, n_sp + 1), mar = c(4, 4, 3, 2))

  for (sp in 1:n_sp) {
    ts_S     <- sapply(1:t, function(i) sum(abundanceS[,,,i,sp]))
    ts_J     <- sapply(1:t, function(i) sum(abundanceJ[,,,i,sp]))
    ts_A     <- sapply(1:t, function(i) sum(abundanceA[,,,i,sp]))
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
  dev.flush()
  Sys.sleep(sleeptime)
}

# ── Post-hoc abundance plot ───────────────────────────────────────────────────
# Usage:
#   plot_abundance(result)                           # all species panels + totals, last t
#   plot_abundance(result, t = 20)                   # up to year 20
#   plot_abundance(result, species_specific = FALSE) # totals panel only

plot_abundance <- function(result, t = NULL, species_specific = TRUE) {
  state      <- result$state
  t_max      <- if (is.null(t)) length(result$totalabundanceA) else t
  stage_cols <- scico::scico(3, palette = "lipari", begin = 0.2, end = 0.8)
  sp_cols    <- if (state$n_species == 1)
    scico::scico(3, palette = "lipari", begin = 0.3, end = 0.7)[2]
  else
    scico::scico(state$n_species, palette = "lipari", begin = 0.2, end = 0.8)

  totalS <- result$totalabundanceS
  totalJ <- result$totalabundanceJ
  totalA <- result$totalabundanceA
  total  <- totalS + totalJ + totalA

  n_panels <- if (species_specific) state$n_species + 1L else 1L
  par(mfrow = c(1, n_panels), mar = c(4, 4, 3, 2))

  if (species_specific) {
    for (sp in 1:state$n_species) {
      ts_S     <- sapply(1:t_max, function(i) sum(result$abundanceS[,,,i,sp]))
      ts_J     <- sapply(1:t_max, function(i) sum(result$abundanceJ[,,,i,sp]))
      ts_A     <- sapply(1:t_max, function(i) sum(result$abundanceA[,,,i,sp]))
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
    sp_name <- state$species_ids[sp]
    col     <- sp_cols[sp]
    idx_S <- which(result$abundanceS[,,,t,sp] > 0, arr.ind=TRUE)
    idx_J <- which(result$abundanceJ[,,,t,sp] > 0, arr.ind=TRUE)
    idx_A <- which(result$abundanceA[,,,t,sp] > 0, arr.ind=TRUE)
    if (nrow(idx_S) > 0)
      rows[[length(rows)+1]] <- data.frame(
        x=idx_S[,1], y=idx_S[,2], z=idx_S[,3],
        species=sp_name, stage="S", color=col, symbol="circle",
        stringsAsFactors=FALSE)
    if (nrow(idx_J) > 0)
      rows[[length(rows)+1]] <- data.frame(
        x=idx_J[,1], y=idx_J[,2], z=idx_J[,3],
        species=sp_name, stage="J", color=col, symbol="diamond",
        stringsAsFactors=FALSE)
    if (nrow(idx_A) > 0)
      rows[[length(rows)+1]] <- data.frame(
        x=idx_A[,1], y=idx_A[,2], z=idx_A[,3],
        species=sp_name, stage="A", color=col, symbol="square",
        stringsAsFactors=FALSE)
  }
  
  if (length(rows) == 0) {
    message("No individuals to plot at t=", t)
    return(invisible(NULL))
  }
  
  df     <- do.call(rbind, rows)
  traces <- split(df, paste0(df$species, "_", df$stage))
  
  fig <- plotly::plot_ly()
  for (tr in traces) {
    fig <- plotly::add_trace(fig,
                             data=tr, x=~x, y=~y, z=~z,
                             type="scatter3d", mode="markers",
                             name=paste0(tr$species[1], " ", tr$stage[1]),
                             marker=list(symbol=tr$symbol[1], color=tr$color[1],
                                         size=6, opacity=0.85))
  }
  fig <- plotly::layout(fig,
                        title=paste0("Abundance (t=", t, ")"),
                        scene=list(xaxis=list(title="x"),
                                   yaxis=list(title="y"),
                                   zaxis=list(title="height tier")))
  print(fig)
  invisible(fig)
}

# ── Spin-up ───────────────────────────────────────────────────────────────────

run_spinup <- function(state, n_gens = 5, Visualize = TRUE,
                       carCap = 1, sleeptime = 0.2,
                       visualize_dispersion = FALSE) {
  
  xDim      <- state$xDim
  yDim      <- state$yDim
  zDim      <- state$zDim
  n_species <- state$n_species
  
  spinupS <- array(0L, dim=c(xDim, yDim, zDim, n_gens, n_species))
  spinupJ <- array(0L, dim=c(xDim, yDim, zDim, n_gens, n_species))
  spinupA <- array(0L, dim=c(xDim, yDim, zDim, n_gens, n_species))
  totalS  <- numeric(n_gens)
  totalJ  <- numeric(n_gens)
  totalA  <- numeric(n_gens)
  
  for (i in 1:nrow(state$site_obs)) {
    idx <- state$coord_to_idx(state$site_obs$lon[i],
                              state$site_obs$lat[i],
                              state$site_obs$hSnapped[i])
    x <- idx[1]; y <- idx[2]; z <- idx[3]
    sp <- state$sp_index[state$site_obs$FinalID[i]]
    if (x >= 1 && x <= xDim && y >= 1 && y <= yDim && state$landscape[x, y, z])
      spinupA[x, y, z, 1, sp] <- 1L
  }
  log_msg(sprintf("Spin-up: placed %d observed individuals (%d species)",
                  sum(spinupA[,,,1,]), n_species))
  
  for (gen in 1:n_gens) {
    log_msg(sprintf("Spin-up generation %d/%d", gen, n_gens))
    Disp    <- run_pass1_disperse(state, spinupA, gen)
    
    if (visualize_dispersion) {
      disp_xy <- apply(Disp, c(1,2), sum)
      fields::image.plot(disp_xy,
                         col=scico::scico(25, palette="lajolla"),
                         main=paste0("Dispersal (gen ", gen, ") | seeds: ", sum(Disp)),
                         xlab="x", ylab="y")
      dev.flush(); Sys.sleep(sleeptime)
    }
    
    spinupS <- run_pass2_establish(state, spinupS, Disp, gen)
    result  <- run_pass3_survive_grow(state, spinupS, spinupJ, spinupA, gen)
    spinupS <- result$S; spinupJ <- result$J; spinupA <- result$A
    
    if (gen < n_gens) {
      spinupS[,,,gen+1,] <- spinupS[,,,gen+1,] + spinupS[,,,gen,]
      spinupJ[,,,gen+1,] <- spinupJ[,,,gen+1,] + spinupJ[,,,gen,]
      spinupA[,,,gen+1,] <- spinupA[,,,gen+1,] + spinupA[,,,gen,]
    }
    
    totalS[gen] <- sum(spinupS[,,,gen,])
    totalJ[gen] <- sum(spinupJ[,,,gen,])
    totalA[gen] <- sum(spinupA[,,,gen,])
    
    if (Visualize)
      .plot_live(state, spinupS, spinupJ, spinupA,
                 totalS, totalJ, totalA, gen, carCap, sleeptime)
    
    log_msg(sprintf("Gen %d complete: S=%d J=%d A=%d total=%d",
                    gen, totalS[gen], totalJ[gen], totalA[gen],
                    totalS[gen]+totalJ[gen]+totalA[gen]))
  }
  
  log_msg(sprintf("Spin-up complete: S=%d J=%d A=%d",
                  sum(spinupS[,,,n_gens,]), sum(spinupJ[,,,n_gens,]),
                  sum(spinupA[,,,n_gens,])))
  list(S=spinupS[,,,n_gens,], J=spinupJ[,,,n_gens,], A=spinupA[,,,n_gens,],
       last_disp=Disp)
}

# ── Main wrapper ──────────────────────────────────────────────────────────────

runcolonization <- function(site, niches, canopy_grid, models, valid_per_model,
                            timesteps=50, resolution=10, carCap=1,
                            maxDisp=5, stochastic=FALSE,
                            Visualize=TRUE, sleeptime=0.2,
                            visualize_dispersion=FALSE, spinup=5,
                            parameters, allsites=FALSE) {
  
  state     <- init_colonization(site, niches, canopy_grid, models, valid_per_model,
                                 resolution, carCap, maxDisp, params=parameters,
                                 allsites=allsites)
  xDim      <- state$xDim
  yDim      <- state$yDim
  zDim      <- state$zDim
  n_species <- state$n_species
  
  abundanceS      <- array(0L, dim=c(xDim, yDim, zDim, timesteps, n_species))
  abundanceJ      <- array(0L, dim=c(xDim, yDim, zDim, timesteps, n_species))
  abundanceA      <- array(0L, dim=c(xDim, yDim, zDim, timesteps, n_species))
  totalabundanceS <- numeric(timesteps)
  totalabundanceJ <- numeric(timesteps)
  totalabundanceA <- numeric(timesteps)
  
  sp <- run_spinup(state, Visualize=Visualize, n_gens=spinup,
                   carCap=carCap, sleeptime=sleeptime,
                   visualize_dispersion=visualize_dispersion)
  abundanceS[,,,1,] <- sp$S
  abundanceJ[,,,1,] <- sp$J
  abundanceA[,,,1,] <- sp$A
  totalabundanceS[1] <- sum(abundanceS[,,,1,])
  totalabundanceJ[1] <- sum(abundanceJ[,,,1,])
  totalabundanceA[1] <- sum(abundanceA[,,,1,])
  log_msg(sprintf("Starting population: %d S, %d J, %d A",
                  totalabundanceS[1], totalabundanceJ[1], totalabundanceA[1]))
  
  Disp <- sp$last_disp
  
  for (t in 1:(timesteps-1)) {
    Disp       <- run_pass1_disperse(state, abundanceA, t)
    abundanceS <- run_pass2_establish(state, abundanceS, Disp, t)
    result     <- run_pass3_survive_grow(state, abundanceS, abundanceJ, abundanceA, t)
    abundanceS <- result$S; abundanceJ <- result$J; abundanceA <- result$A
    
    abundanceS[,,,t+1,] <- abundanceS[,,,t+1,] + abundanceS[,,,t,]
    abundanceJ[,,,t+1,] <- abundanceJ[,,,t+1,] + abundanceJ[,,,t,]
    abundanceA[,,,t+1,] <- abundanceA[,,,t+1,] + abundanceA[,,,t,]
    
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
       totalabundanceS=totalabundanceS, totalabundanceJ=totalabundanceJ,
       totalabundanceA=totalabundanceA,
       heights=state$heights, species_ids=state$species_ids,
       xDim=xDim, yDim=yDim, zDim=zDim, n_species=n_species,
       state=state, last_disp=Disp)
}
