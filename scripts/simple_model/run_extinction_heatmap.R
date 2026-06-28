# run_extinction_heatmap.R
# Fine-scale extinction heatmap across establishment prob × repro rate.
# Holds all other parameters at baseline and scans the threshold region
# where populations transition from extinction to persistence.
#
# Grid: p_est in 0.005–0.05 (step 0.005), repro in 10–150 (step 10)
# → 10 × 15 = 150 combos × N_REPS replicates
#
# Lizeth Estévez Tobar — University of Bonn, 2026
# ─────────────────────────────────────────────────────────────────────────────

source("scripts/simple_colonization.R")
library(parallel)
library(ggplot2)

N_REPS     <- 5
TIMESTEPS  <- 50

# ── Parameter grid ────────────────────────────────────────────────────────────
p_est_vals <- seq(0.005, 0.05, by = 0.005)   # 10 levels
repro_vals <- seq(10,   150,   by = 10)       # 15 levels

grid <- expand.grid(
  p_est     = p_est_vals,
  repro     = repro_vals,
  rep       = seq_len(N_REPS)
)

# ── Run one combination ───────────────────────────────────────────────────────
run_one <- function(i) {
  p <- params
  p$establishment_prob <- grid$p_est[i]
  p$repro_rate         <- grid$repro[i]
  p$timesteps          <- TIMESTEPS

  # suppress per-timestep cat() output from run_simple_colonization
  result <- suppressMessages(
    capture.output(
      out <- run_simple_colonization(p),
      type = "output"
    )
  )

  # extinct if no adults remain in final 10 years
  final_A <- out$abundA[, , , (TIMESTEPS - 9):TIMESTEPS]
  extinct  <- all(final_A == 0)
  total_A  <- sum(out$abundA[, , , TIMESTEPS])

  list(p_est = grid$p_est[i], repro = grid$repro[i],
       extinct = extinct, final_A = total_A)
}

cat("Running", nrow(grid), "simulations on", detectCores(), "cores...\n")
results <- mclapply(seq_len(nrow(grid)), run_one,
                    mc.cores = max(1, detectCores() - 1))

# ── Summarise ─────────────────────────────────────────────────────────────────
res_df <- do.call(rbind, lapply(results, as.data.frame))

heatmap_df <- aggregate(
  cbind(ext_rate = extinct, mean_A = final_A) ~ p_est + repro,
  data = res_df,
  FUN  = mean
)

# ── Plot ──────────────────────────────────────────────────────────────────────
p_heat <- ggplot(heatmap_df, aes(x = factor(repro), y = factor(p_est),
                                  fill = ext_rate)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = paste0(round(ext_rate * 100), "%")),
            size = 3, colour = "black") +
  scale_fill_gradient2(
    low      = "#2c7bb6",
    mid      = "#ffffbf",
    high     = "#d7191c",
    midpoint = 0.5,
    limits   = c(0, 1),
    name     = "Extinction\nrate"
  ) +
  scale_x_discrete(breaks = as.character(seq(10, 150, by = 20))) +
  labs(
    title = "Extinction rate — establishment prob × repro rate",
    subtitle = paste0("Each cell = mean of ", N_REPS,
                      " runs × 50 yr; all other params at baseline"),
    x = expression(paste("Repro rate  (seeds adult"^{-1}, " yr"^{-1}, ")")),
    y     = "Establishment probability"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid  = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

out_path <- "output/extinction_heatmap_fine.png"
ggsave(out_path, p_heat, width = 9, height = 7, dpi = 180)
cat("Saved:", out_path, "\n")

print(p_heat)
