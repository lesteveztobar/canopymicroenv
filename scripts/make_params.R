source("scripts/paths.R")

params_ppoll <- list(
    # ── s(z, e): survival ────────────────────────────────────────────────────
    beta0S  = -0.24,  beta0J  =  0.41,  beta0A  =  1.73,
    beta1   =  0.10,
    z_S_min =  0.0,  z_S_max =  1.0,
    z_J_min =  1.0,  z_J_max =  7.0,
    z_A_min =  7.0,  z_A_max = 20.0,
    # ── g(z'|z, e): growth / stage transitions ────────────────────────────
    psi0S        = -3.30,  psi0J        = -2.70,
    beta_precip  =  3e-4,  beta_rh      =  0.010,
    sigma        =  0.10,  delta_z_base =  0.80,
    cost_repro   =  0.50,
    # ── p_r(z) × f_s(z): fecundity ──────────────────────────────────────
    p_poll  = c(0.05, 0.15, 0.30, 0.50, 0.70),  p_germ  = 0.001,  p_s1 = 0.45,
    # ── d(x'|x): dispersal ────────────────────────────────────────────────
    # canopy_z omitted — it's site-specific and run_colonization_onesite.R
    # always overwrites it with the per-site mean canopy height.
    lambda = 1,  Ut = 1
)

out_path <- file.path(PARAMS_DIR, "p_poll.rds")
saveRDS(params_ppoll, out_path)


params_beta0A <- list(
    # ── s(z, e): survival ────────────────────────────────────────────────────
    beta0S  = -0.24,  beta0J  =  0.41,  beta0A  =  c(0.50, 1.00, 1.73, 2.50, 3.50),
    beta1   =  0.10,
    z_S_min =  0.0,  z_S_max =  1.0,
    z_J_min =  1.0,  z_J_max =  7.0,
    z_A_min =  7.0,  z_A_max = 20.0,
    # ── g(z'|z, e): growth / stage transitions ────────────────────────────
    psi0S        = -3.30,  psi0J        = -2.70,
    beta_precip  =  3e-4,  beta_rh      =  0.010,
    sigma        =  0.10,  delta_z_base =  0.80,
    cost_repro   =  0.50,
    # ── p_r(z) × f_s(z): fecundity ──────────────────────────────────────
    p_poll  = 0.30,  p_germ  = 0.001,  p_s1 = 0.45,
    # ── d(x'|x): dispersal ──────────────────────────────────────────────
    # canopy_z omitted — site-specific, overwritten by run_colonization_onesite.R
    lambda = 1,  Ut = 1
)


out_path <- file.path(PARAMS_DIR, "beta0A.rds")
saveRDS(params_beta0A, out_path)

params_pgerm <- list(
# ── s(z, e): survival ────────────────────────────────────────────────────
beta0S  = -0.24,  beta0J  =  0.41,  beta0A  =  1.73,
beta1   =  0.10,
z_S_min =  0.0,  z_S_max =  1.0,
z_J_min =  1.0,  z_J_max =  7.0,
z_A_min =  7.0,  z_A_max = 20.0,
# ── g(z'|z, e): growth / stage transitions ────────────────────────────
psi0S        = -3.30,  psi0J        = -2.70,
beta_precip  =  3e-4,  beta_rh      =  0.010,
sigma        =  0.10,  delta_z_base =  0.80,
cost_repro   =  0.50,
# ── p_r(z) × f_s(z): fecundity ──────────────────────────────────────
p_poll  = 0.30,  p_germ  = c(0.0001, 0.0005, 0.001, 0.003, 0.005),  p_s1 = 0.45,
# ── d(x'|x): dispersal ──────────────────────────────────────────────
# canopy_z omitted — site-specific, overwritten by run_colonization_onesite.R
lambda = 1,  Ut = 1
)

out_path <- file.path(PARAMS_DIR, "p_germ.rds")
saveRDS(params_pgerm, out_path)

params_costrepro <- list(
# ── s(z, e): survival ────────────────────────────────────────────────────
beta0S  = -0.24,  beta0J  =  0.41,  beta0A  =  1.73,
beta1   =  0.10,
z_S_min =  0.0,  z_S_max =  1.0,
z_J_min =  1.0,  z_J_max =  7.0,
z_A_min =  7.0,  z_A_max = 20.0,
# ── g(z'|z, e): growth / stage transitions ────────────────────────────
psi0S        = -3.30,  psi0J        = -2.70,
beta_precip  =  3e-4,  beta_rh      =  0.010,
sigma        =  0.10,  delta_z_base =  0.80,
cost_repro   =  c(0.20, 0.40, 0.60, 0.80, 1.00),
# ── p_r(z) × f_s(z): fecundity ──────────────────────────────────────
p_poll  = 0.30,  p_germ  = 0.001,  p_s1 = 0.45,
# ── d(x'|x): dispersal ──────────────────────────────────────────────
# canopy_z omitted — site-specific, overwritten by run_colonization_onesite.R
lambda = 1,  Ut = 1
)

out_path <- file.path(PARAMS_DIR, "cost_repro.rds")
saveRDS(params_costrepro, out_path)

params_betarh <- list(
# ── s(z, e): survival ────────────────────────────────────────────────────
beta0S  = -0.24,  beta0J  =  0.41,  beta0A  =  1.73,
beta1   =  0.10,
z_S_min =  0.0,  z_S_max =  1.0,
z_J_min =  1.0,  z_J_max =  7.0,
z_A_min =  7.0,  z_A_max = 20.0,
# ── g(z'|z, e): growth / stage transitions ────────────────────────────
psi0S        = -3.30,  psi0J        = -2.70,
beta_precip  =  3e-4,  beta_rh      =  c(0.002, 0.005, 0.010, 0.020, 0.040),
sigma        =  0.10,  delta_z_base =  0.80,
cost_repro   =  0.50,
# ── p_r(z) × f_s(z): fecundity ──────────────────────────────────────
p_poll  = 0.30,  p_germ  = 0.001,  p_s1 = 0.45,
# ── d(x'|x): dispersal ──────────────────────────────────────────────
# canopy_z omitted — site-specific, overwritten by run_colonization_onesite.R
lambda = 1,  Ut = 1
)

out_path <- file.path(PARAMS_DIR, "beta_rh.rds")
saveRDS(params_betarh, out_path)

params_betaprecip <- list(
# ── s(z, e): survival ────────────────────────────────────────────────────
beta0S  = -0.24,  beta0J  =  0.41,  beta0A  =  1.73,
beta1   =  0.10,
z_S_min =  0.0,  z_S_max =  1.0,
z_J_min =  1.0,  z_J_max =  7.0,
z_A_min =  7.0,  z_A_max = 20.0,
# ── g(z'|z, e): growth / stage transitions ────────────────────────────
psi0S        = -3.30,  psi0J        = -2.70,
beta_precip  =  c(1e-4, 2e-4, 3e-4, 6e-4, 1e-3),  beta_rh      =  0.010,
sigma        =  0.10,  delta_z_base =  0.80,
cost_repro   =  0.50,
# ── p_r(z) × f_s(z): fecundity ──────────────────────────────────────
p_poll  = 0.30,  p_germ  = 0.001,  p_s1 = 0.45,
# ── d(x'|x): dispersal ──────────────────────────────────────────────
# canopy_z omitted — site-specific, overwritten by run_colonization_onesite.R
lambda = 1,  Ut = 1
)

out_path <- file.path(PARAMS_DIR, "beta_precip.rds")
saveRDS(params_betaprecip, out_path)

