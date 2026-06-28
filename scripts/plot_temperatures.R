library(plotly)
library(ggplot2)
library(abind)

env <- readRDS("data/processed/microenv_Maquipucuna.rds")

# env[[height]]$tmax$Tz is a numeric array [nrow, ncol, ntimesteps].
# Collapse timesteps by mean → [nrow, ncol] per height, then stack heights.

heights_m <- as.numeric(sub("h", "", names(env)))

build_volume <- function(env, day_type = "tmax", var = "Tz") {
  layers <- lapply(env, function(h) {
    r <- h[[day_type]][[var]]
    if (length(dim(r)) == 3) apply(r, c(1, 2), mean, na.rm = TRUE)
    else r
  })
  abind::abind(layers, along = 3)  # [nrow, ncol, nheight]
}

vol_tmax <- build_volume(env, "tmax", "Tz")

nrow_r  <- dim(vol_tmax)[1]
ncol_r  <- dim(vol_tmax)[2]
nheight <- dim(vol_tmax)[3]

# Use real lon/lat if spatial reference was saved alongside the microenv,
# otherwise fall back to pixel indices.
if (!is.null(env$.spatial)) {
  e        <- env$.spatial$ext
  x_coords <- seq(e$xmin, e$xmax, length.out = ncol_r)
  y_coords <- seq(e$ymax, e$ymin, length.out = nrow_r)  # north→south
  x_lab <- "Longitude"; y_lab <- "Latitude"
} else {
  x_coords <- seq_len(ncol_r)
  y_coords <- seq_len(nrow_r)
  x_lab <- "col (W→E)"; y_lab <- "row (S→N)"
}

# ── Option A: interactive 3D volume (plotly) ──────────────────────────────────
grid <- expand.grid(xi = x_coords, yi = y_coords, zi = seq_len(nheight))
grid$z  <- heights_m[grid$zi]
grid$Tz <- vol_tmax[cbind(nrow_r + 1 - grid$yi, grid$xi, grid$zi)]
grid    <- grid[!is.na(grid$Tz), ]

fig <- plot_ly(
  data = grid,
  x = ~xi, y = ~yi, z = ~z,
  value = ~Tz,
  type  = "volume",
  isomin    = quantile(grid$Tz, 0.05),
  isomax    = quantile(grid$Tz, 0.95),
  opacity   = 0.15,
  surface   = list(count = 10),
  colorscale = list(
    c(0,    "#313695"),
    c(0.25, "#74add1"),
    c(0.5,  "#ffffbf"),
    c(0.75, "#f46d43"),
    c(1,    "#a50026")
  ),
  colorbar = list(title = "T (°C)")
) |>
  layout(
    title = "Air temperature 3D volume — Maquipucuna (warmest representative day)",
    scene = list(
      xaxis = list(title = x_lab),
      yaxis = list(title = y_lab),
      zaxis = list(title = "Height (m)")
    )
  )

print(fig)

# ── Option B: static faceted heatmap per height ───────────────────────────────
long_df <- do.call(rbind, lapply(seq_along(env), function(i) {
  r <- env[[i]]$tmax$Tz
  m <- if (length(dim(r)) == 3) apply(r, c(1, 2), mean, na.rm = TRUE) else r
  df <- as.data.frame(as.table(m))
  names(df) <- c("row", "col", "Tz")
  df$row    <- if (!is.null(env$.spatial)) y_coords[as.integer(df$row)] else as.integer(df$row)
  df$col    <- if (!is.null(env$.spatial)) x_coords[as.integer(df$col)] else as.integer(df$col)
  df$height <- heights_m[i]
  df
}))
long_df <- long_df[!is.na(long_df$Tz), ]

# ── Side-view heatmap: space (col, W→E) × height ─────────────────────────────
# Collapse rows by mean so each column becomes one horizontal position.
# Build a parallel all-pixels data frame for the statistical test.
side_df <- do.call(rbind, lapply(seq_along(env), function(i) {
  r <- env[[i]]$tmax$Tz
  m <- if (length(dim(r)) == 3) apply(r, c(1, 2), mean, na.rm = TRUE) else r
  col_means <- colMeans(m, na.rm = TRUE)
  data.frame(col = seq_along(col_means), height = heights_m[i], Tz = col_means)
}))
side_df <- side_df[!is.na(side_df$Tz), ]

# All pixel values per height — more power for the test than col means alone
all_px <- do.call(rbind, lapply(seq_along(env), function(i) {
  r <- env[[i]]$tmax$Tz
  m <- if (length(dim(r)) == 3) apply(r, c(1, 2), mean, na.rm = TRUE) else r
  data.frame(height = factor(heights_m[i]), Tz = as.vector(m))
}))
all_px <- all_px[!is.na(all_px$Tz), ]

# ── Kruskal-Wallis: are any height levels different? ─────────────────────────
kw <- kruskal.test(Tz ~ height, data = all_px)
kw_label <- sprintf("Kruskal-Wallis: χ²(%.0f) = %.2f, p %s",
                    kw$parameter, kw$statistic,
                    ifelse(kw$p.value < 0.001, "< 0.001",
                           sprintf("= %.3f", kw$p.value)))
message(kw_label)

# ── Dunn post-hoc pairwise comparisons (Holm correction) ─────────────────────
if (!requireNamespace("dunn.test", quietly = TRUE)) install.packages("dunn.test")
dunn_res <- dunn.test::dunn.test(all_px$Tz, all_px$height,
                                  method = "holm", kw = FALSE, label = FALSE)
dunn_df <- data.frame(
  comparison = dunn_res$comparisons,
  p_adj      = dunn_res$P.adjusted,
  sig        = ifelse(dunn_res$P.adjusted < 0.001, "***",
               ifelse(dunn_res$P.adjusted < 0.01,  "**",
               ifelse(dunn_res$P.adjusted < 0.05,  "*", "ns")))
)
message("\nPairwise Dunn tests (Holm-adjusted):")
print(dunn_df[order(dunn_df$p_adj), ], row.names = FALSE)

# ── Compact letter display for the boxplot annotation ────────────────────────
if (!requireNamespace("multcompView", quietly = TRUE)) install.packages("multcompView")
p_mat <- setNames(dunn_res$P.adjusted, dunn_res$comparisons)
cld   <- multcompView::multcompLetters(p_mat, threshold = 0.05)$Letters
cld_df <- data.frame(height = as.numeric(names(cld)), letter = cld)

# ── Plot A: side-view heatmap ─────────────────────────────────────────────────
p_heat <- ggplot(side_df, aes(x = col, y = height, fill = Tz)) +
  geom_raster() +
  scale_fill_gradientn(
    colours  = c("#313695", "#74add1", "#ffffbf", "#f46d43", "#a50026"),
    name     = "T (°C)", na.value = "grey90"
  ) +
  labs(x = "Horizontal position (col, W→E)", y = "Height (m)",
       title = "Vertical cross-section — Maquipucuna (warmest representative day)",
       subtitle = kw_label) +
  theme_minimal(base_size = 11)

# ── Plot B: boxplot per height with CLD letters ───────────────────────────────
p_box <- ggplot(all_px, aes(x = Tz, y = height, group = height)) +
  geom_boxplot(aes(fill = after_stat(middle)), width = 0.06, outlier.size = 0.5) +
  geom_text(data = cld_df,
            aes(x = max(all_px$Tz, na.rm = TRUE), y = height, label = letter),
            hjust = -0.2, size = 3.5, inherit.aes = FALSE) +
  scale_fill_gradientn(
    colours = c("#313695", "#74add1", "#ffffbf", "#f46d43", "#a50026"),
    guide   = "none"
  ) +
  labs(x = "T (°C)", y = NULL,
       caption = "Letters = Dunn post-hoc (Holm); shared letter → no significant difference") +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_blank())

# ── Combine side-by-side ──────────────────────────────────────────────────────
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)
p_heat + p_box + plot_layout(widths = c(3, 1))

