# simple_colonization.R
# Simplified 3D epiphyte colonization model for teaching / class report.
#
# The model represents a forest as a 3D voxel grid. Trees are placed randomly
# and give the grid structure (which voxels exist, which Johansson zone they
# belong to). Habitat suitability per voxel is derived from the zone — outer
# crown zones are most suitable for epiphytes, trunk zones least. The population
# dynamic is a two-stage (seedling S / adult A) annual loop with dispersal.
#
# No microclimate data are used. All parameters are plain numbers.
#
# Grid default: 50 × 50 × 20 voxels, each 10 m × 10 m × 2 m
#   → 500 m × 500 m footprint, 0–40 m height = ~25 ha of cloud forest
#
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────

# ── Parameters ────────────────────────────────────────────────────────────────

params <- list(
  # Landscape dimensions
  xDim       = 50,    # voxels W–E
  yDim       = 50,    # voxels S–N
  zDim       = 20,    # height tiers
  resolution = 10,    # m per voxel (horizontal)
  max_height = 40,    # m — top of the tallest tier

  # Forest structure — a single random draw sets n_trees at runtime.
  # Range calibrated to primary TMCF at ~1400 m (Myster 2017: 272–324 trees/ha
  # for stems ≥10 cm dsh; 25 ha grid → 6800–8100 canopy trees).
  # We use a reduced count because not all grid cells map to a tree crown.
  n_trees_min  = 4000,
  n_trees_max  = 6000,
  mean_hgt     = 30,    # m — mean tree height
  sd_hgt       =  4,    # m
  mean_crown_r =  3.5,  # m — mean crown radius
  sd_crown_r   =  1.0,

  # Habitat suitability by Johansson zone (0 = bare air, 1 = best epiphyte habitat)
  # Zones 4–5 (crown) are most suitable; zones 1–2 (trunk base) least.
  suit_by_zone = c(
    `0` = 0.00,   # outside any tree
    `1` = 0.10,   # trunk base
    `2` = 0.20,   # lower trunk
    `3` = 0.35,   # upper trunk / first branches
    `4` = 0.75,   # inner crown
    `5` = 0.85    # outer crown
  ),

  # Population dynamics
  # repro_rate = effective seeds dispersed per adult per year.
  # Already collapses flowering rate, fruit number, and pollination success —
  # so this is NOT total seed output (which would be ~20,000) but the realistic
  # number that leave the parent plant and enter the dispersal kernel.
  # ~100 is conservative for a reproducing Maxillariinae adult.
  repro_rate        = 100,    # effective seeds dispersed per adult per year
  maxDisp           = 6,      # max dispersal distance (voxels, horizontal)
  maxDispZ          = 3,      # max vertical dispersal (tiers)
  # establishment_prob × suitability = P(seed → seedling).
  # With repro_rate=100 and suit~0.7, set this low so net recruitment
  # is ~1–3 seedlings per adult per year at good sites.
  establishment_prob = 0.02,  # P(seed establishes | lands in valid voxel)
  survival_S        = 0.55,   # annual seedling survival (suitability-scaled)
  survival_A        = 0.85,   # annual adult survival   (suitability-scaled)
  maturation_prob   = 0.20,   # probability S → A per year
  carCap            = 8,      # max individuals (all stages) per voxel

  # Simulation
  timesteps  = 40,
  n_founders = 30    # adults placed at random valid voxels to seed the run
)

# ── Forest builder ────────────────────────────────────────────────────────────

# Place trees randomly on the grid and return:
#   landscape [xDim, yDim, zDim]  — TRUE where a voxel is inside a tree
#   zone      [xDim, yDim, zDim]  — Johansson zone 1–5 (0 = not in tree)
#   suitability [xDim, yDim, zDim] — habitat quality 0–1
build_forest_simple <- function(p) {
  xDim <- p$xDim; yDim <- p$yDim; zDim <- p$zDim
  heights <- seq(0, p$max_height, length.out = zDim + 1)[-1]  # midpoints of tiers

  n_trees <- round(runif(1, p$n_trees_min, p$n_trees_max))
  trees <- data.frame(
    x       = sample(1:xDim, n_trees, replace = TRUE),
    y       = sample(1:yDim, n_trees, replace = TRUE),
    height  = pmax(2, rnorm(n_trees, p$mean_hgt,    p$sd_hgt)),
    crown_r = pmax(0.5, rnorm(n_trees, p$mean_crown_r, p$sd_crown_r))
  )
  trees$crown_r_cells <- trees$crown_r / p$resolution

  landscape   <- array(FALSE, dim = c(xDim, yDim, zDim))
  zone        <- array(0L,    dim = c(xDim, yDim, zDim))

  for (ti in seq_len(nrow(trees))) {
    tx <- trees$x[ti]; ty <- trees$y[ti]
    th <- trees$height[ti]; cr <- trees$crown_r_cells[ti]

    x_range <- max(1, tx - ceiling(cr)):min(xDim, tx + ceiling(cr))
    y_range <- max(1, ty - ceiling(cr)):min(yDim, ty + ceiling(cr))

    for (x in x_range) for (y in y_range) {
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

        in_tree <- if (jzone <= 2) horiz_dist == 0
                   else {
                     cf <- (rel_h - 0.5) / 0.5
                     horiz_dist <= cr * sin(cf * pi)
                   }
        if (!in_tree) next

        landscape[x, y, z] <- TRUE
        if (jzone > zone[x, y, z]) zone[x, y, z] <- jzone
      }
    }
  }

  suitability <- array(0, dim = c(xDim, yDim, zDim))
  for (z_val in 0:5)
    suitability[zone == z_val] <- p$suit_by_zone[as.character(z_val)]

  list(landscape = landscape, zone = zone, suitability = suitability,
       heights = heights, n_trees = n_trees)
}

# ── Dispersal ─────────────────────────────────────────────────────────────────

# Spread seeds from one adult voxel. All n_seeds are drawn at once (vectorized)
# and counted into Disp via tabulate — no per-seed loop.
disperse_simple <- function(x, y, z, n_seeds, xDim, yDim, zDim,
                            maxDisp, maxDispZ, Disp) {
  tx <- x + sample(-maxDisp:maxDisp, n_seeds, replace = TRUE)
  ty <- y + sample(-maxDisp:maxDisp, n_seeds, replace = TRUE)
  tz <- z + sample(-maxDispZ:maxDispZ, n_seeds, replace = TRUE)
  ok <- tx >= 1 & tx <= xDim & ty >= 1 & ty <= yDim & tz >= 1 & tz <= zDim
  if (!any(ok)) return(Disp)
  # convert valid (tx,ty,tz) to linear indices and count with tabulate
  idx <- (tx[ok] - 1L) * yDim * zDim + (ty[ok] - 1L) * zDim + tz[ok]
  Disp <- Disp + array(tabulate(idx, nbins = xDim * yDim * zDim),
                        dim = c(xDim, yDim, zDim))
  Disp
}

# ── Main simulation loop ──────────────────────────────────────────────────────

run_simple_colonization <- function(p = params, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  cat("Building forest...\n")
  forest <- build_forest_simple(p)
  landscape   <- forest$landscape
  suitability <- forest$suitability
  heights     <- forest$heights
  xDim <- p$xDim; yDim <- p$yDim; zDim <- p$zDim
  T    <- p$timesteps

  valid_voxels <- which(landscape, arr.ind = TRUE)
  cat(sprintf("  Trees placed: %d | Valid voxels: %d (%.1f%% of grid)\n",
              forest$n_trees, nrow(valid_voxels),
              100 * nrow(valid_voxels) / prod(dim(landscape))))

  # Abundance arrays [xDim, yDim, zDim, T]
  abundS <- array(0L, dim = c(xDim, yDim, zDim, T))
  abundA <- array(0L, dim = c(xDim, yDim, zDim, T))

  # Place founders at random valid voxels
  founder_idx <- valid_voxels[sample(nrow(valid_voxels),
                                     min(p$n_founders, nrow(valid_voxels))), ]
  for (i in seq_len(nrow(founder_idx)))
    abundA[founder_idx[i,1], founder_idx[i,2], founder_idx[i,3], 1] <- 1L

  totalS <- integer(T); totalA <- integer(T)
  totalS[1] <- sum(abundS[,,,1]); totalA[1] <- sum(abundA[,,,1])
  # side-view dispersal history: sum Disp over x → [yDim, zDim, T]
  disp_side <- array(0L, dim = c(yDim, zDim, T))
  cat(sprintf("  Founders: %d adults placed\n", totalA[1]))

  for (t in seq_len(T - 1)) {

    # ── Dispersal: adults produce seeds ────────────────────────────────────────
    Disp <- array(0L, dim = c(xDim, yDim, zDim))
    nonzero_A <- which(abundA[,,,t] > 0, arr.ind = TRUE)
    for (i in seq_len(nrow(nonzero_A))) {
      x <- nonzero_A[i,1]; y <- nonzero_A[i,2]; z <- nonzero_A[i,3]
      n_seeds <- rpois(1, abundA[x,y,z,t] * p$repro_rate)
      if (n_seeds > 0)
        Disp <- disperse_simple(x, y, z, n_seeds, xDim, yDim, zDim,
                                p$maxDisp, p$maxDispZ, Disp)
    }
    disp_side[,,t] <- apply(Disp, c(2, 3), sum)  # collapse x → [yDim, zDim]

    # ── Establishment: vectorized over all cells with seeds ───────────────────
    # Only consider valid landscape cells that have seed rain and aren't full.
    can_establish <- Disp > 0 & landscape &
                     (abundS[,,,t] + abundA[,,,t]) < p$carCap
    if (any(can_establish)) {
      space     <- as.integer(p$carCap - abundS[,,, t] - abundA[,,,t])
      p_est_arr <- p$establishment_prob * suitability   # [xDim,yDim,zDim]
      new_s     <- pmin(
        array(rbinom(prod(dim(Disp)), as.vector(Disp), as.vector(p_est_arr)),
              dim = dim(Disp)),
        space
      )
      new_s[!can_establish] <- 0L
      abundS[,,,t+1] <- abundS[,,,t+1] + new_s
    }

    # ── Survival and maturation: vectorized over all occupied cells ───────────
    has_S <- abundS[,,,t] > 0
    has_A <- abundA[,,,t] > 0

    if (any(has_S)) {
      nS_vec    <- abundS[,,,t]
      surv_S    <- array(rbinom(prod(dim(nS_vec)), as.vector(nS_vec),
                                as.vector(p$survival_S * suitability)),
                         dim = dim(nS_vec))
      matures   <- array(rbinom(prod(dim(surv_S)), as.vector(surv_S),
                                p$maturation_prob),
                         dim = dim(surv_S))
      abundS[,,,t+1] <- abundS[,,,t+1] + surv_S - matures
      abundA[,,,t+1] <- abundA[,,,t+1] + matures
    }

    if (any(has_A)) {
      nA_vec  <- abundA[,,,t]
      surv_A  <- array(rbinom(prod(dim(nA_vec)), as.vector(nA_vec),
                               as.vector(p$survival_A * (0.5 + 0.5 * suitability))),
                        dim = dim(nA_vec))
      abundA[,,,t+1] <- abundA[,,,t+1] + surv_A
    }

    # Enforce carCap — trim excess adults (rare edge case)
    total_next <- abundS[,,,t+1] + abundA[,,,t+1]
    over       <- total_next > p$carCap
    if (any(over))
      abundA[,,,t+1] <- pmax(0L, abundA[,,,t+1] - pmax(0L, total_next - p$carCap)) * 1L

    totalS[t+1] <- sum(abundS[,,,t+1])
    totalA[t+1] <- sum(abundA[,,,t+1])
    cat(sprintf("t=%2d | S=%4d  A=%4d  total=%4d | seeds dispersed=%d\n",
                t+1, totalS[t+1], totalA[t+1], totalS[t+1]+totalA[t+1], sum(Disp)))
  }

  list(abundS = abundS, abundA = abundA,
       totalS = totalS, totalA = totalA,
       disp_side = disp_side,
       landscape = landscape, zone = forest$zone,
       suitability = suitability, heights = heights,
       params = p)
}

# ── Run (skipped when sourced by simple_experiments.R) ───────────────────────

if (!isTRUE(getOption("simple_colonization_no_run"))) {
result <- run_simple_colonization(p = params, seed = 42)

# ── Plots ─────────────────────────────────────────────────────────────────────

library(ggplot2)
T <- params$timesteps

# Plot 1: population trajectory
traj_df <- data.frame(
  t     = rep(1:T, 2),
  N     = c(result$totalS, result$totalA),
  stage = rep(c("Seedling", "Adult"), each = T)
)
p_traj <- ggplot(traj_df, aes(x = t, y = N, colour = stage)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = c(Seedling = "#4dac26", Adult = "#d01c8b")) +
  labs(x = "Year", y = "Abundance", colour = NULL,
       title = "Simple colonization model — population trajectory") +
  theme_minimal(base_size = 12)
print(p_traj)

# Plot 2: side-view heatmap of adult abundance at final timestep
# Collapse x (W–E) by summing → abundance as function of (y, height)
final_A <- result$abundA[,,,T]
side_df <- do.call(rbind, lapply(seq_along(result$heights), function(z) {
  data.frame(
    y      = seq_len(params$yDim),
    height = result$heights[z],
    N      = rowSums(final_A[,,z])   # sum over x
  )
}))
side_df <- side_df[side_df$N > 0 | TRUE, ]   # keep all rows for continuous fill

p_side <- ggplot(side_df, aes(x = y, y = height, fill = N)) +
  geom_raster() +
  scale_fill_gradientn(
    colours  = c("white", "#f7fbff", "#6baed6", "#08519c"),
    name     = "Adults",
    na.value = "grey95"
  ) +
  labs(x = "S → N position (voxels)", y = "Height (m)",
       title = "Adult abundance — vertical cross-section (final timestep)") +
  theme_minimal(base_size = 12)
print(p_side)

# Plot 3: mean adult abundance per Johansson zone over time
zone_arr <- result$zone
zone_time <- do.call(rbind, lapply(1:T, function(t) {
  do.call(rbind, lapply(1:5, function(z_val) {
    mask <- zone_arr == z_val
    data.frame(t = t, zone = z_val,
               N = sum(result$abundA[,,,t][mask], na.rm = TRUE))
  }))
}))
zone_time$zone <- factor(zone_time$zone, labels = c(
  "1 Trunk base", "2 Lower trunk", "3 Upper trunk", "4 Inner crown", "5 Outer crown"))

p_zone <- ggplot(zone_time, aes(x = t, y = N, colour = zone)) +
  geom_line(linewidth = 1) +
  scale_colour_brewer(palette = "RdYlGn", direction = 1) +
  labs(x = "Year", y = "Adult abundance", colour = "Johansson zone",
       title = "Adults by canopy zone over time") +
  theme_minimal(base_size = 12)
print(p_zone)
} # end run guard

