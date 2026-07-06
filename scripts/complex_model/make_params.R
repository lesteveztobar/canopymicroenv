source("scripts/complex_model/paths.R")

# Founders per species (see get_colonization.R::run_spinup()), used as a
# fixed baseline across every experiment below except the n_founders sweep
# itself. Update this once that sweep identifies a value that actually
# persists, then regenerate all the other params files so they pick it up.
N_FOUNDERS_DEFAULT <- 30

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
    lambda = 1,  Ut = 1,
    # ── spin-up ────────────────────────────────────────────────────────────
    n_founders = N_FOUNDERS_DEFAULT
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
    lambda = 1,  Ut = 1,
    n_founders = N_FOUNDERS_DEFAULT
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
lambda = 1,  Ut = 1,
n_founders = N_FOUNDERS_DEFAULT
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
lambda = 1,  Ut = 1,
n_founders = N_FOUNDERS_DEFAULT
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
lambda = 1,  Ut = 1,
n_founders = N_FOUNDERS_DEFAULT
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
lambda = 1,  Ut = 1,
n_founders = N_FOUNDERS_DEFAULT
)

out_path <- file.path(PARAMS_DIR, "beta_precip.rds")
saveRDS(params_betaprecip, out_path)

# ── Establishment: niche tolerance pad ──────────────────────────────────────────
# Sweeps how much to pad each species' realized climate niche (see
# get_niche()/niche_match() in get_colonization.R) beyond its raw observed
# range. 0 = no padding (often a single point given how few observations
# most species have — get_niche()'s min_width floor keeps pad meaningful
# even then); larger values progressively relax the niche gate.
params_nichepad <- list(
beta0S  = -0.24,  beta0J  =  0.41,  beta0A  =  1.73,
beta1   =  0.10,
z_S_min =  0.0,  z_S_max =  1.0,
z_J_min =  1.0,  z_J_max =  7.0,
z_A_min =  7.0,  z_A_max = 20.0,
psi0S        = -3.30,  psi0J        = -2.70,
beta_precip  =  3e-4,  beta_rh      =  0.010,
sigma        =  0.10,  delta_z_base =  0.80,
cost_repro   =  0.50,
p_poll  = 0.30,  p_germ  = 0.001,  p_s1 = 0.45,
# canopy_z omitted — site-specific, overwritten by run_colonization_onesite.R
lambda = 1,  Ut = 1,
n_founders = N_FOUNDERS_DEFAULT,
niche_pad = c(0, 0.005, 0.01, 0.05)
)

out_path <- file.path(PARAMS_DIR, "niche_pad.rds")
saveRDS(params_nichepad, out_path)

# ── Spin-up: founder count ──────────────────────────────────────────────────────
# Sweeps n_founders per species (see get_colonization.R::run_spinup()) —
# decoupled from field observation count. Run this experiment alone first:
# the fecundity math (reproduce(), get_colonization.R:226-234) means expected
# seed output per adult per year is tiny even at the most generous swept
# p_poll/p_germ/p_s1 combination, so persistence may hinge on founder count
# more than on any vital rate. Whatever value here actually persists should
# become the new N_FOUNDERS_DEFAULT above, then regenerate the other params
# files so every other experiment (including the reproduction factorial)
# uses a founder count that's actually viable.
params_nfounders <- list(
beta0S  = -0.24,  beta0J  =  0.41,  beta0A  =  1.73,
beta1   =  0.10,
z_S_min =  0.0,  z_S_max =  1.0,
z_J_min =  1.0,  z_J_max =  7.0,
z_A_min =  7.0,  z_A_max = 20.0,
psi0S        = -3.30,  psi0J        = -2.70,
beta_precip  =  3e-4,  beta_rh      =  0.010,
sigma        =  0.10,  delta_z_base =  0.80,
cost_repro   =  0.50,
p_poll  = 0.30,  p_germ  = 0.001,  p_s1 = 0.45,
# canopy_z omitted — site-specific, overwritten by run_colonization_onesite.R
lambda = 1,  Ut = 1,
n_founders = seq(5, 100, by = 5)
)

out_path <- file.path(PARAMS_DIR, "n_founders.rds")
saveRDS(params_nfounders, out_path)

# ── Reproduction factorial: n_founders x p_poll x p_germ x p_s1 ────────────────
# Full factorial across founder count and the three components of fecundity
# (pollination success, mycorrhizal germination, first-year seed-to-seedling
# survival) — 5^4 = 625 combinations, run via run_factorial_experiment().
# n_founders alone (5-100, at literature-default reproduction) came back
# 100% extinct at every level tested — the fecundity formula makes expected
# seed output per adult per year ~0.00006 at baseline, so no realistic
# founder count rescues it there. Crossing n_founders into the factorial
# instead of testing it in isolation finds the actual joint threshold
# directly, mirroring how the simple model's Fig. 5 (report.pdf sec. 3.3)
# found the p_est x lambda threshold in one factorial rather than
# sequentially. n_founders levels span from what we already know fails
# (50, 100) up to what the math says should be comfortably enough at the
# most generous reproduction corner tested (400, 800): at
# p_poll=0.7/p_germ=0.005/p_s1=0.75, expected seeds/adult/yr ~0.0011, so
# n_founders=200 gives ~5 expected seeds over 30 years, n_founders=800 ~30.
params_reprofactorial <- list(
beta0S  = -0.24,  beta0J  =  0.41,  beta0A  =  1.73,
beta1   =  0.10,
z_S_min =  0.0,  z_S_max =  1.0,
z_J_min =  1.0,  z_J_max =  7.0,
z_A_min =  7.0,  z_A_max = 20.0,
psi0S        = -3.30,  psi0J        = -2.70,
beta_precip  =  3e-4,  beta_rh      =  0.010,
sigma        =  0.10,  delta_z_base =  0.80,
cost_repro   =  0.50,
p_poll  = c(0.05, 0.15, 0.30, 0.50, 0.70),
p_germ  = c(0.0001, 0.0005, 0.001, 0.003, 0.005),
p_s1    = c(0.15, 0.30, 0.45, 0.60, 0.75),
# canopy_z omitted — site-specific, overwritten by run_colonization_onesite.R
lambda = 1,  Ut = 1,
n_founders = c(50, 100, 200, 400, 800)
)

out_path <- file.path(PARAMS_DIR, "reproduction_factorial.rds")
saveRDS(params_reprofactorial, out_path)

# ── Best case: everything pushed as favourable as possible ─────────────────────
# Every lever more generous than anything tested so far in the OAT sweeps or
# the factorial — not meant to be realistic, just to answer "can the model
# persist at all?" n_reps replicates (see run_replicated() in
# get_colonization.R) rather than a single run, since a single zero-
# establishment outcome can't distinguish "genuinely blocked" from "just an
# unlucky stochastic draw" — see project memory on the niche/founder/
# capacity-saturation investigation for why this matters here specifically.
params_bestcase <- list(
beta0S  = -0.24,  beta0J  =  0.41,  beta0A  =  3.50,   # max tested (best survival)
beta1   =  0.10,
z_S_min =  0.0,  z_S_max =  1.0,
z_J_min =  1.0,  z_J_max =  7.0,
z_A_min =  7.0,  z_A_max = 20.0,
psi0S        = -3.30,  psi0J        = -2.70,
beta_precip  =  1e-3,  beta_rh      =  0.040,          # max tested (best growth)
sigma        =  0.10,  delta_z_base =  0.80,
cost_repro   =  0.20,                                   # min tested (least reproduction cost)
p_poll  = 0.90,  p_germ  = 0.01,  p_s1 = 0.90,           # beyond anything tested so far
# canopy_z omitted — site-specific, overwritten by run_colonization_onesite.R
lambda = 1,  Ut = 1,
n_founders = 1000,                                       # beyond the 800 max tested
niche_pad  = 0.5,                                         # far beyond the 0.05 max tested
n_reps     = 5
)

out_path <- file.path(PARAMS_DIR, "best_case.rds")
saveRDS(params_bestcase, out_path)
