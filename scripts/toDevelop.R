# To-do (or to-use) list of stuff that are not currently in development
runcolonization <- function(timesteps = 30, xDim = 10, yDim = length(unique(climgrid$hReq)), carCap = 200,
                            stochastic = FALSE, Visualize = TRUE, sleeptime = 0.2,
                            maxDisp = 4, climgrid = NULL, heights = seq(0.5, 25, by = 0.5)) {
  # ── Sanity check ─────────────────────────────────────────────────────────────
  if (!is.null(climgrid) && yDim != length(heights))
    stop("yDim must equal length(heights) when climgrid is provided.")
  
  
  #INITIALIZATION:
  # --- Landscape ---    
  landscape <- array(data = 0, dim = c(yDim, xDim, timesteps))
  for (i in 1:timesteps){
    for (c in 1:(xDim/2)){
      landscape[, c, i] <- dnorm(n = yDim, mean = mean(climgrid$Temp), sd = sd(climgrid$Temp))
    }
  }
  # --- Habitat Preference ---
  happyTsp1 <- runif(n=1, min = min(climgrid$Temp), max = max(climgrid$Temp))
  happyTsp2 <- runif(n=1, min = min(climgrid$Temp), max = max(climgrid$Temp))
  happyTsp3 <- runif(n=1, min = min(climgrid$Temp), max = max(climgrid$Temp))
  
  # --- Abundance per stage ---
  abundanceS <- array (data = 0, dim= c(yDim, xDim, timesteps))
  abundanceJ <- array (data = 0, dim= c(yDim, xDim, timesteps))
  abundanceA <- array (data = 0, dim= c(yDim, xDim, timesteps))
  totalabundanceS <- vector(mode = "numeric", length = timesteps)
  totalabundanceJ <- vector(mode = "numeric", length = timesteps)
  totalabundanceA <- vector(mode = "numeric", length = timesteps)
  
  # --- Growth parameters (IPM-style, Merow et al. 2014) ---
  beta0GrowthS <- 0.22
  beta0GrowthJ <- 0.45
  beta0GrowthA <- 1.40
  beta1Growth   <- 0.8  # placeholder, size-dependent growth dampening
  sigma         <- 0.5  # placeholder residual variance — calibrate later
  
  # --- Reproduction (Winkler et al. 2009; Taylor & Bruns 1999; Zotz 1998) ---
  p_flower <- 0.3
  p_poll <- 0.05
  p_germ <- 0.001 # germination probability, low due to mycorrhizal interactions not possible to be modeled, so accounted per (Taylor & Bruns 1999)
  p_s1 <- 0.83 # first-year seedling survival as per (Zotz 1998)
  
  repInfo <- matrix(0, nrow = 2, ncol = 3,
                    dimnames = list(c("n_seeds", "repRate"),
                                    c("Sp1", "Sp2", "Sp3")))
  repInfo["n_seeds", ] <- round(runif(3, 500, 5000))
  repInfo["repRate", ] <- p_flower * p_poll * repInfo["n_seeds", ] * p_germ * p_s1
  repRate <- mean(repInfo["repRate", ])
  
  # --- Survival parameters (Raventós et al. 2015; Olaya-Arenas et al. 2011) ---
  beta0Seedling <- -1.5  # low baseline survival, to be calibrated
  beta0Juvenile <- -0.5  # medium
  beta0Adult    <-  1.0  # high, from Raventós structure
  beta1         <-  0.105 # from Raventós et al. 2015
  
  # --- Dispersal parameters (Murren & Ellison 1998; Winkler et al. 2009) ---
  Uf_base  <- 2    # wind velocity above canopy (m/s) — overridden by climgrid when available
  Ut       <- 1    # terminal seed velocity (m/s)
  lambda   <- 1    # spread factor
  canopy_z <- 25   # mean canopy height (m)
  a        <- 23   # canopy openness parameter 
  dispersalmatrix <- array(data = 0, dim = c((yDim+2*maxDisp), (xDim+2*maxDisp), timesteps))
  
  # ── Seed initial abundance where suitability > 19 ───────────────────────────
  for (r in 1:yDim) {
    for (c in 1:xDim) {
      if (landscape[r, c, 1] > 19) {
        abundanceS[r, c, 1] <- runif(1, 1, 1 + carCap * landscape[r, c, 1])
        abundanceJ[r, c, 1] <- runif(1, 1, 1 + carCap * landscape[r, c, 1])
        abundanceA[r, c, 1] <- runif(1, 1, 1 + carCap * landscape[r, c, 1])
      }
    }
  }
  totalabundanceS[1] <- sum(abundanceS[, , 1])
  totalabundanceJ[1] <- sum(abundanceJ[, , 1])
  totalabundanceA[1] <- sum(abundanceA[, , 1])
  
  # ── Eco loop ─────────────────────────────────────────────────────────────────
  for (t in 1:(timesteps-1)){
    for (r in 1:yDim){
      for (c in 1:xDim){
        # -- Look up microclimate for this height tier and x position -----------
        if (!is.null(climgrid)) {
          h    <- heights[r]
          clim <- climgrid[climgrid$hReq == h & climgrid$x == c, ]
          temp    <- clim$Temp
          relhum  <- clim$RelHum
          wind    <- clim$WindSpeed
          swdown  <- clim$SWdown
        } else {
          h      <- heights[r]
          temp   <- 20
          relhum <- 85
          wind   <- Uf_base
          swdown <- 200
        }
        
        # -- Growth (size state per stage, drawn fresh each cell/step) ----------
        epsilon = rnorm(1, 0, sigma)
        zSeedling <- array (data = 0, dim= c(yDim, xDim, timesteps))
        zSeedling[1,1,1] <- runif(1, min = 0, max = 0.22)
        zJuvenile <- array (data = 0, dim= c(yDim, xDim, timesteps))
        zJuvenile[1,1,1] <- runif(1, min = 0.22, max = 0.45)
        zAdult <- array (data = 0, dim= c(yDim, xDim, timesteps))
        zAdult[1,1,1] <- runif(1, min = 0.45, max = 1.4)
        
        if (t<timesteps-1){
          zSeedling[t+1] <- beta0GrowthS + beta1Growth * zSeedling[1,1,1] + epsilon
          zJuvenile[t+1] <- beta0GrowthJ + beta1Growth * zJuvenile[1,1,1] + epsilon
          zAdult[t+1]    <- beta0GrowthA + beta1Growth * zAdult[1,1,1] + epsilon
        }
        
        # -- Survival (logistic, size + climate effects) ------------------------
        # relhum penalty on seedlings (Olaya-Arenas et al. 2011: dry_penalty ~ low RH)
        dry_penalty <- (100 - relhum) / 100
        sSeedling <- 1 / (1 + exp(-(beta0Seedling + beta1 * zSeedling[1,1,1] - dry_penalty)))
        # no effect on juveniles
        sJuvenile <- 1 / (1 + exp(-(beta0Juvenile + beta1 * zJuvenile[1,1,1])))
        # temperature lag effect on adults (Olaya-Arenas 2011: temp_penalty)
        temp_penalty <- ifelse(temp > 23, (temp - 23) * 0.05, 0)
        sAdult    <- 1 / (1 + exp(-(beta0Adult    + beta1 * zAdult[1,1,1] - temp_penalty)))
        
        # -- Total N and carrying capacity for this cell ------------------------
        N <- abundanceS[r, c, t] * sSeedling +
          abundanceJ[r, c, t] * sJuvenile +
          abundanceA[r, c, t] * sAdult
        
        # K modulated by radiation (more light = more resources)
        K <- carCap * landscape[r, c, t] * (swdown / 400)
        
        # -- Offspring via Ricker -----------------------------------------------
        offspring <- if (!stochastic) {
          ricker(N = N, r = repRate, K = K)
        } else {
          stochRicker(N = N, r = repRate, K = K)
        }
        
        # -- Dispersal 
        #first I kill them before disperse them 
        offspring <- offspring * p_s1
        if (offspring>0){
          # wind velocity regulates dispersal 
          meanDispersal <- round(climgrid$WindSpeed[r] * exp(-(a * climgrid$hReq[c] / canopy_z))) / (lambda * Ut)
          meanDispersal <- max(1, min(meanDispersal, maxDisp))
          
          dispersalmatrix[, , t] <- indDisp(
            maxDispersal = meanDispersal,
            r    = r,
            c    = c,
            indN = offspring,
            Disp = dispersalmatrix[, , t]
          )
        }
      }
    }
    # ── Advance abundance: recruits enter as seedlings ────────────────────────
    for (r in 1:yDim) {
      for (c in 1:xDim) {
        new_recruits <- dispersalmatrix[r + maxDisp, c + maxDisp, t]
        abundanceS[r, c, t + 1] <- new_recruits
        abundanceJ[r, c, t + 1] <- abundanceS[r, c, t] * plogis(beta0Seedling)
        abundanceA[r, c, t + 1] <- abundanceJ[r, c, t] * plogis(beta0Juvenile)
      }
    }
    
    totalabundanceS[t + 1] <- sum(abundanceS[, , t + 1])
    totalabundanceJ[t + 1] <- sum(abundanceJ[, , t + 1])
    totalabundanceA[t + 1] <- sum(abundanceA[, , t + 1])
    totalabundance <- totalabundanceS + totalabundanceJ + totalabundanceA
    
    # ── Visualisation ─────────────────────────────────────────────────────────
    if (Visualize) {
      par(mfrow = c(1, 3))
      image.plot(landscape[, , t], col = pink_viridis(25),
                 xlab = "x position", ylab = "Canopy height tier")
      mtext(paste0("Landscape (t = ", t + 1, ")"))
      image.plot(abundanceA[, , t], col = pink_viridis(25),
                 xlab = "x position", ylab = "Canopy height tier")
      mtext(paste0("Adults (t = ", t + 1, ")"))
      plot(totalabundance[1:(t + 1)], type = "b", col = "pink",
           xlab = "Time step", ylab = "Total abundance",
           main = "Total abundance", las = 1)
      abline(h = carCap * xDim * yDim, col = "red")
      Sys.sleep(sleeptime)
    }
  }
  return(list(
    landscape       = landscape,
    abundanceS      = abundanceS,
    abundanceJ      = abundanceJ,
    abundanceA      = abundanceA,
    totalabundanceS = totalabundanceS,
    totalabundanceJ = totalabundanceJ,
    totalabundanceA = totalabundanceA
  ))
}