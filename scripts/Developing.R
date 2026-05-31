# Developing.R
# Canopy colonization model — vertical niche partitioning of epiphytic Maxillariinae
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────

library(fields)      # image.plot for visualization
library(viridis)     # colour palettes
library(dplyr)
library(readr)

source("scripts/functions.R")
source("scripts/helper_functions.R")
source("scripts/paths.R")

# ── Load microclimate outputs ─────────────────────────────────────────────────
models <- readRDS(file.path(BASE_DIR, "data/processed/pointmodelv1.rds"))

# ── Actual development space  ─────────────────────────────────────────────────
# nothing works (:

#### testingggg
# ── dummy inputs ──
xDim <- 10; yDim <- 10; zDim <- 5
timesteps <- 20
repRate   <- 0.5
p_germ    <- 0.3
maxDisp   <- 3
carCap    <- 10

lajolla    <- scico::scico(25, palette = "lajolla")
stage_cols <- scico::scico(3, palette = "batlow", begin = 0.2, end = 0.8)

clim <- data.frame(
  temp      = runif(672, 35.3, 36.1),
  relhum    = runif(672, 84, 86),
  swdown    = c(rep(0, 336), rep(200, 336)),
  windspeed = rep(1.5, 672),
  winddir   = rep(180, 672),
  precip    = rep(0.1, 672)
)

clim_month       <- clim[1:56, ]
mean_swdown_site <- mean(clim$swdown[clim$swdown > 0])

abundanceS <- array(0L, dim = c(xDim, yDim, zDim, timesteps))
abundanceJ <- array(0L, dim = c(xDim, yDim, zDim, timesteps))
abundanceA <- array(0L, dim = c(xDim, yDim, zDim, timesteps))
totalS     <- numeric(timesteps)
totalJ     <- numeric(timesteps)
totalA     <- numeric(timesteps)
disp_mat   <- array(0L, dim = c(xDim + 2*maxDisp, yDim + 2*maxDisp, zDim + 2*maxDisp))

# place 5 adults randomly at t=1
for (i in 1:10) {
  x <- sample(1:xDim, 1); y <- sample(1:yDim, 1); z <- sample(1:zDim, 1)
  abundanceA[x, y, z, 1] <- abundanceA[x, y, z, 1] + 1L
}
totalA[1] <- sum(abundanceA[,,,1])
cat("t=1 | A=", totalA[1], "\n")

# Plot initial adult placement
par(mfrow = c(1, 1), mar = c(4, 4, 3, 5))
for (z in 1:zDim) {
  fields::image.plot(abundanceA[,,z,1], col = lajolla,
                     main = paste0("Adults h=", z, " (t=1)"),
                     xlab = "x", ylab = "y",
                     zlim = c(0, 5),
                     axis.args = list(at = 0:5, labels = 0:5))
  idx <- which(abundanceA[,,z,1] > 0, arr.ind = TRUE)
  if (nrow(idx) > 0)
    points((idx[,1] - 0.5) / xDim, (idx[,2] - 0.5) / yDim,
           pch=15, col=stage_cols[3], cex=2)
  legend("topright", legend="A", pch=15, col=stage_cols[3], bg="white", cex=0.8)
}


t <- 1

# ── eco loop ──
for (t in 1:(timesteps - 1)) {
  disp_mat[] <- 0L
  
  # pass 1: reproduce + disperse
  for (x in 1:xDim) for (y in 1:yDim) for (z in 1:zDim) {
    N <- abundanceA[x, y, z, t]
    if (N == 0) next
    seeds <- reproduce(N, repRate)
    if (seeds == 0) next
    disp_mat <- indDisp(x, y, z, seeds, 180, 2, disp_mat, maxDisp)
  }
  cat("t=", t+1, "| seeds dispersed:", sum(disp_mat), "\n")
  
  # pass 2: establish
  for (x in 1:xDim) for (y in 1:yDim) for (z in 1:zDim) {
    if (abundanceS[x, y, z, t] > 0) next
    new_seeds <- disp_mat[x + maxDisp, y + maxDisp, z + maxDisp]
    if (new_seeds == 0) next
    abundanceS[x, y, z, t+1] <- establish(new_seeds, clim, mean_swdown_site, p_germ)
  }
  
  # pass 3: survive + grow
  for (x in 1:xDim) for (y in 1:yDim) for (z in 1:zDim) {
    abundanceS[x,y,z,t] <- survive(abundanceS[x,y,z,t], "S", clim_month)
    abundanceJ[x,y,z,t] <- survive(abundanceJ[x,y,z,t], "J", clim_month)
    abundanceA[x,y,z,t] <- survive(abundanceA[x,y,z,t], "A", clim_month)
    grown <- grow(abundanceS[x,y,z,t], abundanceJ[x,y,z,t], abundanceA[x,y,z,t])
    abundanceS[x,y,z,t] <- grown$nS
    abundanceJ[x,y,z,t] <- grown$nJ
    abundanceA[x,y,z,t] <- grown$nA
  }
  
  # carry forward
  abundanceS[,,,t+1] <- abundanceS[,,,t+1] + abundanceS[,,,t]
  abundanceJ[,,,t+1] <- abundanceJ[,,,t+1] + abundanceJ[,,,t]
  abundanceA[,,,t+1] <- abundanceA[,,,t+1] + abundanceA[,,,t]
  
  totalS[t+1] <- sum(abundanceS[,,,t+1])
  totalJ[t+1] <- sum(abundanceJ[,,,t+1])
  totalA[t+1] <- sum(abundanceA[,,,t+1])
  
  cat(sprintf("t=%d | S=%d J=%d A=%d total=%d\n", t+1,
              totalS[t+1], totalJ[t+1], totalA[t+1],
              totalS[t+1]+totalJ[t+1]+totalA[t+1]))
  
  # visualize
  total <- totalS + totalJ + totalA
  
  for (z in 1:zDim) {
    par(mfrow = c(1, 3), mar = c(4, 4, 3, 5))
    
    # Panel 1: dispersal at this height
    disp_z <- disp_mat[,,z]
    fields::image.plot(disp_z, col = lajolla,
                       main = paste0("Dispersal h=", z, " (t=", t+1, ")"),
                       xlab = "x", ylab = "y",
                       zlim = c(0, max(disp_z, 1)),
                       axis.args = list(at = 0:max(disp_z, 1),
                                        labels = as.integer(0:max(disp_z, 1))))
    
    # Panel 2: abundance at this height with stage points
    fields::image.plot(abundanceA[,,z,t+1], col = lajolla,
                       main = paste0("Abundance h=", z, " (t=", t+1, ")"),
                       xlab = "x", ylab = "y",
                       zlim = c(0, 5),
                       axis.args = list(at = 0:5, labels = 0:5))
    
    idx_S <- which(abundanceS[,,z,t+1] > 0, arr.ind = TRUE)
    idx_J <- which(abundanceJ[,,z,t+1] > 0, arr.ind = TRUE)
    idx_A <- which(abundanceA[,,z,t+1] > 0, arr.ind = TRUE)
    if (nrow(idx_S) > 0) points((idx_S[,1]-0.5)/xDim, (idx_S[,2]-0.5)/yDim,
                                pch=16, col=stage_cols[1], cex=0.5+abundanceS[idx_S[,1],idx_S[,2],z,t+1]*0.3)
    if (nrow(idx_J) > 0) points((idx_J[,1]-0.5)/xDim, (idx_J[,2]-0.5)/yDim,
                                pch=17, col=stage_cols[2], cex=0.5+abundanceJ[idx_J[,1],idx_J[,2],z,t+1]*0.3)
    if (nrow(idx_A) > 0) points((idx_A[,1]-0.5)/xDim, (idx_A[,2]-0.5)/yDim,
                                pch=15, col=stage_cols[3], cex=0.5+abundanceA[idx_A[,1],idx_A[,2],z,t+1]*0.3)
    legend("topright", legend=c("S","J","A"), pch=c(16,17,15),
           col=stage_cols, bg="white", cex=0.7)
    
    # Panel 3: abundance lines (same for all z)
    plot(total[1:(t+1)], type="b", col="black", lwd=2,
         ylim=c(0, max(total, 1)),
         xlab="Year", ylab="Abundance", main="Total abundance", las=1)
    lines(totalS[1:(t+1)], type="b", col=stage_cols[1], pch=16)
    lines(totalJ[1:(t+1)], type="b", col=stage_cols[2], pch=17)
    lines(totalA[1:(t+1)], type="b", col=stage_cols[3], pch=15)
    abline(h=carCap*xDim*yDim*zDim, col="red", lty=2)
    legend("topleft", legend=c("Total","S","J","A"),
           col=c("black",stage_cols), lty=1, pch=c(NA,16,17,15), cex=0.7)
    
    dev.flush()
    Sys.sleep(0.1)
  }
  Sys.sleep(0.3)
}


####
# Test viability across different initial population sizes
n_reps    <- 20      # repetitions per N
n_init_vals <- c(1, 2, 5, 10, 20, 50, 100)
timesteps_test <- 30

results <- data.frame()

for (Ninit in n_init_vals) {
  survived_S <- 0
  survived_J <- 0
  survived_A <- 0
  
  for (rep in 1:n_reps) {
    
    abundanceA_t <- array(0L, dim = c(xDim, yDim, zDim, timesteps_test))
    abundanceS_t <- array(0L, dim = c(xDim, yDim, zDim, timesteps_test))
    abundanceJ_t <- array(0L, dim = c(xDim, yDim, zDim, timesteps_test))
    disp_t       <- array(0L, dim = c(xDim+2*maxDisp, yDim+2*maxDisp, zDim+2*maxDisp))
    
    for (i in 1:Ninit) {
      x <- sample(1:xDim, 1); y <- sample(1:yDim, 1); z <- sample(1:zDim, 1)
      abundanceA_t[x, y, z, 1] <- abundanceA_t[x, y, z, 1] + 1L
    }
    
    for (t in 1:(timesteps_test - 1)) {
      disp_t[] <- 0L
      for (x in 1:xDim) for (y in 1:yDim) for (z in 1:zDim) {
        N <- abundanceA_t[x, y, z, t]
        if (N == 0) next
        seeds <- reproduce(N, repRate)
        if (seeds == 0) next
        disp_t <- indDisp(x, y, z, seeds, 180, 2, disp_t, maxDisp)
      }
      for (x in 1:xDim) for (y in 1:yDim) for (z in 1:zDim) {
        new_seeds <- disp_t[x+maxDisp, y+maxDisp, z+maxDisp]
        if (new_seeds == 0) next
        abundanceS_t[x,y,z,t+1] <- establish(new_seeds, clim, mean_swdown_site, p_germ)
      }
      for (x in 1:xDim) for (y in 1:yDim) for (z in 1:zDim) {
        abundanceS_t[x,y,z,t] <- survive(abundanceS_t[x,y,z,t], "S", clim_month)
        abundanceJ_t[x,y,z,t] <- survive(abundanceJ_t[x,y,z,t], "J", clim_month)
        abundanceA_t[x,y,z,t] <- survive(abundanceA_t[x,y,z,t], "A", clim_month)
        grown <- grow(abundanceS_t[x,y,z,t], abundanceJ_t[x,y,z,t], abundanceA_t[x,y,z,t])
        abundanceS_t[x,y,z,t] <- grown$nS
        abundanceJ_t[x,y,z,t] <- grown$nJ
        abundanceA_t[x,y,z,t] <- grown$nA
      }
      abundanceS_t[,,,t+1] <- abundanceS_t[,,,t+1] + abundanceS_t[,,,t]
      abundanceJ_t[,,,t+1] <- abundanceJ_t[,,,t+1] + abundanceJ_t[,,,t]
      abundanceA_t[,,,t+1] <- abundanceA_t[,,,t+1] + abundanceA_t[,,,t]
    }
    
    if (sum(abundanceS_t[,,,timesteps_test]) > 0) survived_S <- survived_S + 1
    if (sum(abundanceJ_t[,,,timesteps_test]) > 0) survived_J <- survived_J + 1
    if (sum(abundanceA_t[,,,timesteps_test]) > 0) survived_A <- survived_A + 1
  }
  
  results <- rbind(results, data.frame(
    Ninit      = Ninit,
    survival_S = survived_S / n_reps,
    survival_J = survived_J / n_reps,
    survival_A = survived_A / n_reps
  ))
  cat(sprintf("Ninit=%d | S=%.0f%% J=%.0f%% A=%.0f%%\n",
              Ninit,
              survived_S/n_reps*100,
              survived_J/n_reps*100,
              survived_A/n_reps*100))
}

# Plot
plot(results$Ninit, results$survival_A * 100, type="b",
     col=stage_cols[3], pch=15, lwd=2,
     ylim=c(0,100),
     xlab="Initial adults", ylab="Survival rate (%)",
     main="Minimum viable population by stage", las=1)
lines(results$Ninit, results$survival_J * 100, type="b", col=stage_cols[2], pch=17, lwd=2)
lines(results$Ninit, results$survival_S * 100, type="b", col=stage_cols[1], pch=16, lwd=2)
abline(h=50, col="red", lty=2)
legend("bottomright", legend=c("A","J","S"), col=stage_cols[3:1],
       pch=c(15,17,16), lty=1, lwd=2, cex=0.8)
