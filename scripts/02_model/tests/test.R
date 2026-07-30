opts <- options(stringsAsFactors = FALSE)

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_msg <- function(...) {
    msg <- paste(..., collapse = "")
    cat(sprintf("[%s] %s\n", timestamp(), msg), file = stdout())
    if (interactive()) flush.console()
    flush(stdout())
}

bench <- function(label, expr) {
    log_msg("Starting ", label, " run")
    t0 <- Sys.time()
    res <- tryCatch(
        eval(expr, envir = parent.frame()),
        error = function(e) {
            log_msg("ERROR in ", label, ": ", conditionMessage(e))
            stop(e)
        }
    )
    log_msg(
        "Finished ", label, " in ",
        round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
        " seconds"
    )
    res
}

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_args, value = TRUE)
if (length(script_file) > 0) {
    script_path <- sub("^--file=", "", script_file[1])
} else {
    stop("Cannot determine script path. Run this script with Rscript.")
}
script_dir <- dirname(normalizePath(script_path))
project_root <- normalizePath(file.path(script_dir, "..", "..", ".."))

# CORRECTED: config/paths.R now lives under scripts/02_model/config/, not
# directly under BASE_DIR -- config/ is a sibling of tests/ within 02_model,
# so it's one level up from script_dir (tests/), not three.
paths_file <- file.path(dirname(script_dir), "config", "paths.R")
if (!file.exists(paths_file)) {
    stop("config/paths.R not found at ", paths_file)
}
source(paths_file)

root_dir <- if (exists("paths") && is.list(paths) && !is.null(paths$root)) {
    normalizePath(paths$root)
} else if (exists("PROJECT_ROOT")) {
    normalizePath(PROJECT_ROOT)
} else {
    project_root
}

scripts_dir <- if (exists("paths") && is.list(paths) && !is.null(paths$scripts)) {
    file.path(root_dir, paths$scripts)
} else {
    file.path(root_dir, "scripts")
}

source_engine_files <- function() {
    engine_dir <- file.path(scripts_dir, "02_model", "engine")
    if (!dir.exists(engine_dir)) stop("Engine directory not found: ", engine_dir)
    files <- sort(list.files(engine_dir, pattern = "\\.R$", full.names = TRUE))
    for (f in files) source(f)

    extras <- c(
        file.path(scripts_dir, "02_model", "setup", "make_params.R")
    )
    for (path in extras) {
        if (!file.exists(path)) {
            stop("Required helper file not found: ", path)
        }
        source(path)
    }
}

find_microenv <- function() {
    pat <- "(la|La).*elenita.*h0\\.5.*\\.rds$"
    found <- list.files(root_dir,
        pattern = pat, recursive = TRUE,
        full.names = TRUE, ignore.case = TRUE
    )
    if (length(found) > 0) {
        return(found[[1]])
    }

    pat2 <- "elenita.*\\.rds$"
    found2 <- list.files(root_dir,
        pattern = pat2, recursive = TRUE,
        full.names = TRUE, ignore.case = TRUE
    )
    if (length(found2) > 0) {
        return(found2[[1]])
    }
    NULL
}

parse_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    opts <- list(
        compare = FALSE,
        microenv = NULL,
        timesteps = 20L,
        spinup = 5L,
        seed = 42L,
        n_founders = 100L,
        site = "La Elenita",
        visual = FALSE,
        output = file.path(project_root, "scripts", "02_model", "tests", "test_results")
    )
    for (arg in args) {
        if (arg == "--compare") {
            opts$compare <- TRUE
        } else if (grepl("^--microenv=", arg)) {
            opts$microenv <- sub("^--microenv=", "", arg)
        } else if (grepl("^--timesteps=", arg)) {
            opts$timesteps <- as.integer(sub("^--timesteps=", "", arg))
        } else if (grepl("^--spinup=", arg)) {
            opts$spinup <- as.integer(sub("^--spinup=", "", arg))
        } else if (grepl("^--seed=", arg)) {
            opts$seed <- as.integer(sub("^--seed=", "", arg))
        } else if (grepl("^--n_founders=", arg)) {
            opts$n_founders <- as.integer(sub("^--n_founders=", "", arg))
        } else if (grepl("^--output=", arg)) opts$output <- sub("^--output=", "", arg)
    }
    opts
}

summary_xy <- function(xy) {
    values <- as.numeric(xy)
    values <- values[is.finite(values)]
    if (length(values) == 0) {
        return(NULL)
    }
    c(
        min = min(values),
        max = max(values),
        range = max(values) - min(values),
        mean = mean(values),
        sd = sd(values),
        median = median(values)
    )
}

inspect_slot <- function(slot, slot_name, hours) {
    if (!is.list(slot) || is.null(slot$Tz)) {
        return(NULL)
    }
    vars <- c("Tz", "relhum", "windspeed", "Rdirdown", "Rdifdown")
    for (v in vars) {
        arr <- slot[[v]]
        if (is.null(arr) || length(dim(arr)) != 3) next
        cat(sprintf("  %s / %s: dim=%s\n", slot_name, v, paste(dim(arr), collapse = "x")))
        for (h in intersect(hours, seq_len(dim(arr)[3]))) {
            stats <- summary_xy(arr[, , h])
            cat(sprintf(
                "    hour %3d: min=%.4f max=%.4f range=%.4f sd=%.4f\n",
                h, stats["min"], stats["max"], stats["range"], stats["sd"]
            ))
        }
        overall <- summary_xy(arr)
        cat(sprintf(
            "    overall: min=%.4f max=%.4f range=%.4f sd=%.4f\n",
            overall["min"], overall["max"], overall["range"], overall["sd"]
        ))
    }
}

inspect_height <- function(height, microenv, hours = c(1, 50, 100, 150, 200, 250)) {
    cat("Inspecting height:", height, "\n")
    h <- load_height(microenv, height)
    if (is.null(h)) {
        cat("  height not found\n")
        return(NULL)
    }
    inspect_slot(h$tmax, "tmax", hours)
    inspect_slot(h$tmin, "tmin", hours)
    invisible(TRUE)
}

.make_monthly_xy <- function(arr) {
    if (length(dim(arr)) == 2) arr <- array(arr, dim = c(dim(arr), 1))
    n <- dim(arr)[3]
    if (n %% 12 == 0 && n > 24) {
        months <- split(seq_len(n), rep(1:12, each = n / 12))
        lapply(months, function(idx) apply(arr[, , idx, drop = FALSE], c(1, 2), mean, na.rm = TRUE))
    } else {
        one <- apply(arr, c(1, 2), mean, na.rm = TRUE)
        rep(list(one), 12)
    }
}

.get_clim_spatial <- function(height, microenv) {
    h <- load_height(microenv, height)
    if (is.null(h)) {
        return(NULL)
    }

    temp_tmax <- .make_monthly_xy(h$tmax$Tz)
    temp_tmin <- .make_monthly_xy(h$tmin$Tz)
    relhum_tmax <- .make_monthly_xy(h$tmax$relhum)
    relhum_tmin <- .make_monthly_xy(h$tmin$relhum)
    wind_tmax <- .make_monthly_xy(h$tmax$windspeed)
    wind_tmin <- .make_monthly_xy(h$tmin$windspeed)
    sw_tmax <- lapply(seq_len(12), function(m) .make_monthly_xy(h$tmax$Rdirdown)[[m]] + .make_monthly_xy(h$tmax$Rdifdown)[[m]])
    sw_tmin <- lapply(seq_len(12), function(m) .make_monthly_xy(h$tmin$Rdirdown)[[m]] + .make_monthly_xy(h$tmin$Rdifdown)[[m]])
    dif_tmax <- .make_monthly_xy(h$tmax$Rdifdown)
    dif_tmin <- .make_monthly_xy(h$tmin$Rdifdown)

    w <- microenv$.weather
    if (!is.null(w) && !is.null(w$obs_time)) {
        precip_by_month <- tapply(w$precip, format(w$obs_time, "%m"), mean, na.rm = TRUE)
        precip_month <- as.numeric(precip_by_month[sprintf("%02d", 1:12)])
        winddir <- mean(w$winddir, na.rm = TRUE)
    } else if (!is.null(w)) {
        precip_month <- rep(mean(w$precip, na.rm = TRUE), 12)
        winddir <- mean(w$winddir, na.rm = TRUE)
    } else {
        precip_month <- rep(NA_real_, 12)
        winddir <- NA_real_
    }

    monthly <- lapply(seq_len(12), function(m) {
        list(
            temp = (temp_tmax[[m]] + temp_tmin[[m]]) / 2,
            relhum = (relhum_tmax[[m]] + relhum_tmin[[m]]) / 2,
            windspeed = (wind_tmax[[m]] + wind_tmin[[m]]) / 2,
            swdown = (sw_tmax[[m]] + sw_tmin[[m]]) / 2,
            difrad = (dif_tmax[[m]] + dif_tmin[[m]]) / 2,
            precip = precip_month[m],
            winddir = winddir
        )
    })

    mean_matrix <- function(lst) {
        Reduce("+", lst) / length(lst)
    }

    list(
        year = list(
            temp = mean_matrix(lapply(monthly, "[[", "temp")),
            relhum = mean_matrix(lapply(monthly, "[[", "relhum")),
            windspeed = mean_matrix(lapply(monthly, "[[", "windspeed")),
            swdown = mean_matrix(lapply(monthly, "[[", "swdown")),
            difrad = mean_matrix(lapply(monthly, "[[", "difrad")),
            precip = mean(unlist(lapply(monthly, "[[", "precip")), na.rm = TRUE),
            winddir = winddir
        ),
        months = monthly
    )
}

build_clim_cache_spatial <- function(microenv, heights = NULL) {
    if (is.null(heights)) {
        heights <- microenv_heights(microenv)
    } else {
        heights <- as.numeric(heights)
    }
    if (length(heights) == 0) {
        stop("No heights supplied to build_clim_cache_spatial()")
    }

    n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1L)))
    log_msg("Building spatial climate cache for ", length(heights), " height(s): ", paste(heights, collapse = ", "))

    clim_by_height <- parallel::mclapply(
        heights,
        function(h) .get_clim_spatial(h, microenv),
        mc.cores = n_cores
    )

    clim_month_by_height <- lapply(clim_by_height, function(cl) {
        if (is.null(cl)) NULL else cl$months
    })

    list(clim_by_height = clim_by_height, clim_month_by_height = clim_month_by_height)
}

niche_axis_score <- function(value, axis) {
    if (is.null(axis) || is.null(value) || is.na(value)) {
        return(NA_real_)
    }
    approx(axis$x, axis$score, xout = value, rule = 2)$y
}

niche_axis_scores_array <- function(clim_values, niche) {
    if (is.null(niche)) {
        return(NULL)
    }
    scores <- lapply(names(niche$axes), function(v) {
        axis <- niche$axes[[v]]
        vals <- as.numeric(clim_values[[v]])
        sc <- approx(axis$x, axis$score, xout = vals, rule = 2)$y
        if (is.matrix(clim_values[[v]])) matrix(sc, nrow = nrow(clim_values[[v]]), ncol = ncol(clim_values[[v]])) else sc
    })
    names(scores) <- names(niche$axes)
    scores
}

niche_raw_score_array <- function(clim_values, niche) {
    if (is.null(niche)) {
        mat <- clim_values[[1]]
        if (is.matrix(mat)) {
            return(matrix(100, nrow = nrow(mat), ncol = ncol(mat)))
        }
        return(100)
    }
    scores <- niche_axis_scores_array(clim_values, niche)
    if (is.null(scores)) {
        return(100)
    }
    if (is.matrix(scores[[1]])) {
        stacked <- simplify2array(scores)
        if (any(stacked <= 0, na.rm = TRUE)) {
            return(matrix(0, nrow = nrow(stacked), ncol = ncol(stacked)))
        }
        exp(apply(log(stacked), c(1, 2), mean, na.rm = TRUE))
    } else {
        x <- do.call(cbind, lapply(scores, as.numeric))
        if (any(x <= 0, na.rm = TRUE)) {
            return(0)
        }
        exp(mean(log(x)))
    }
}

niche_overall_score_array <- function(clim_values, niche) {
    raw <- niche_raw_score_array(clim_values, niche)
    ceiling <- if (!is.null(niche) && !is.null(niche$ceiling)) niche$ceiling else 100
    if (is.matrix(raw)) {
        if (!is.finite(ceiling) || ceiling <= 0) {
            return(raw)
        }
        pmin(100, 100 * raw / ceiling)
    }
    if (!is.finite(ceiling) || ceiling <= 0) {
        return(raw)
    }
    min(100, 100 * raw / ceiling)
}

niche_match_array <- function(clim_values, niche) {
    score <- niche_overall_score_array(clim_values, niche)
    if (is.matrix(score)) {
        return(score / 100)
    }
    score / 100
}

.run_pass1_disperse_spatial <- function(state, abundanceA, size_A, t) {
    p <- state$params
    Disp <- array(0L, dim = c(
        state$xDim + 2L * state$maxDisp,
        state$yDim + 2L * state$maxDisp,
        state$zDim + 2L * state$maxDisp,
        state$n_species
    ))
    for (sp in seq_len(state$n_species)) {
        nonzero <- which(abundanceA[, , , t, sp] > 0, arr.ind = TRUE)
        for (i in seq_len(nrow(nonzero))) {
            x <- nonzero[i, 1]
            y <- nonzero[i, 2]
            z <- nonzero[i, 3]
            if (!state$landscape[x, y, z]) next
            N <- abundanceA[x, y, z, t, sp]
            if (is.na(N) || N == 0L) next
            psb_s <- size_A[x, y, z, t, sp]
            seeds <- reproduce(N,
                s = psb_s,
                p_poll = p$p_poll, p_germ = p$p_germ, p_s1 = p$p_s1
            )
            if (seeds == 0L) next
            clim <- state$clim_spatial_by_height[[z]]
            wind <- if (!is.null(clim$year$windspeed)) clim$year$windspeed[x, y] else mean(clim$windspeed, na.rm = TRUE)
            if (is.na(wind)) wind <- 1
            winddir <- if (!is.null(clim$year$winddir)) clim$year$winddir else mean(clim$winddir, na.rm = TRUE)
            Disp[, , , sp] <- disperse(
                x = x, y = y, z = z, seeds = seeds,
                clim = list(windspeed = wind, winddir = winddir),
                height = state$heights[z], canopy_z = p$canopy_z,
                a = p$a, lambda = p$lambda, Ut = p$Ut,
                maxDisp = state$maxDisp, Disp = Disp[, , , sp]
            )
        }
    }
    Disp
}

.run_pass2_establish_spatial <- function(state, abundanceS, abundanceJ, abundanceA,
                                         size_S, dispersalmatrix, t) {
    p <- state$params
    pad <- state$maxDisp
    tnext <- t + 1L
    if (tnext > dim(abundanceS)[4]) {
        return(list(S = abundanceS, size_S = size_S))
    }
    total_occ <- apply(abundanceS[, , , t, , drop = FALSE], 1:3, sum) +
        apply(abundanceJ[, , , t, , drop = FALSE], 1:3, sum) +
        apply(abundanceA[, , , t, , drop = FALSE], 1:3, sum)
    for (sp in seq_len(state$n_species)) {
        for (zi in seq_len(state$zDim)) {
            clim <- state$clim_spatial_by_height[[zi]]
            if (is.null(clim)) next
            seeds_slice <- dispersalmatrix[
                (pad + 1L):(pad + state$xDim),
                (pad + 1L):(pad + state$yDim),
                zi + pad, sp
            ]
            if (sum(seeds_slice) == 0L) next
            p_est <- niche_match_array(
                list(temp = clim$year$temp, relhum = clim$year$relhum, swdown = clim$year$swdown),
                state$niches_by_species[[state$species_ids[sp]]]
            )
            if (is.matrix(clim$year$relhum)) {
                p_est <- (clim$year$relhum / 100) * pmin(clim$year$swdown / p$mean_swdown_site, 1) * p_est
            } else {
                p_est <- (clim$year$relhum / 100) * pmin(clim$year$swdown / p$mean_swdown_site, 1) * p_est
            }
            can_est <- state$landscape[, , zi] &
                total_occ[, , zi] < state$carCap_voxel[, , zi] &
                seeds_slice > 0L
            if (!any(can_est)) next
            seed_vec <- seeds_slice[can_est]
            prob_vec <- p_est[can_est]
            established_vec <- rbinom(length(seed_vec), as.integer(seed_vec), pmin(pmax(prob_vec, 0), 1))
            space_avail <- pmax(0L, state$carCap_voxel[, , zi][can_est] - total_occ[, , zi][can_est])
            established_vec <- pmin(established_vec, space_avail)
            abundanceS[, , zi, tnext, sp][can_est] <-
                abundanceS[, , zi, tnext, sp][can_est] + established_vec
            size_S[, , zi, tnext, sp][can_est] <- p$s_S_min
        }
    }
    list(S = abundanceS, size_S = size_S)
}

.run_pass3_survive_grow_spatial <- function(state, abundanceS, abundanceJ, abundanceA,
                                            size_S, size_J, size_A, fruited, t) {
    p <- state$params
    xDim <- state$xDim
    yDim <- state$yDim
    for (sp in seq_len(state$n_species)) {
        for (zi in seq_len(state$zDim)) {
            if (sum(abundanceS[, , zi, t, sp]) + sum(abundanceJ[, , zi, t, sp]) + sum(abundanceA[, , zi, t, sp]) == 0L) next
            clim_year <- state$clim_spatial_by_height[[zi]]
            precip_annual <- mean(clim_year$year$precip, na.rm = TRUE) * 8760
            slS <- abundanceS[, , zi, t, sp]
            slJ <- abundanceJ[, , zi, t, sp]
            slA <- abundanceA[, , zi, t, sp]
            sS_slice <- size_S[, , zi, t, sp]
            sJ_slice <- size_J[, , zi, t, sp]
            sa_slice <- size_A[, , zi, t, sp]
            delta_s_total_S <- array(0, dim = c(xDim, yDim))
            delta_s_total_J <- array(0, dim = c(xDim, yDim))
            delta_s_total_A <- array(0, dim = c(xDim, yDim))
            for (month in seq_len(12)) {
                cm <- state$clim_spatial_month_by_height[[zi]][[month]]
                temp <- cm$temp
                relhum <- cm$relhum
                swdown_rel <- if (p$mean_swdown_site > 0) cm$swdown / p$mean_swdown_site else 1.0
                s_S_arr <- survival_logit(
                    "S", sS_slice, temp, relhum, swdown_rel,
                    p$beta0_S, p$beta0_J, p$beta0_A, p$beta1
                )
                s_J_arr <- survival_logit(
                    "J", sJ_slice, temp, relhum, swdown_rel,
                    p$beta0_S, p$beta0_J, p$beta0_A, p$beta1
                )
                s_A_arr <- survival_logit(
                    "A", sa_slice, temp, relhum, swdown_rel,
                    p$beta0_S, p$beta0_J, p$beta0_A, p$beta1
                )
                slS <- array(rbinom(length(slS), as.integer(slS), pmin(pmax(s_S_arr, 0), 1)), dim = dim(slS))
                slJ <- array(rbinom(length(slJ), as.integer(slJ), pmin(pmax(s_J_arr, 0), 1)), dim = dim(slJ))
                slA <- array(rbinom(length(slA), as.integer(slA), pmin(pmax(s_A_arr, 0), 1)), dim = dim(slA))
                p_StoJ <- transition_logit(
                    "S", precip_annual, mean(relhum, na.rm = TRUE),
                    p$psi0S, p$psi0J, p$beta_precip, p$beta_rh
                )
                p_JtoA <- transition_logit(
                    "J", precip_annual, mean(relhum, na.rm = TRUE),
                    p$psi0S, p$psi0J, p$beta_precip, p$beta_rh
                )
                p_StoJ <- pmin(1, pmax(0, p_StoJ + rnorm(1, 0, p$sigma)))
                p_JtoA <- pmin(1, pmax(0, p_JtoA + rnorm(1, 0, p$sigma)))
                n_StoJ <- array(rbinom(length(slS), as.integer(slS), p_StoJ), dim = dim(slS))
                n_JtoA <- array(rbinom(length(slJ), as.integer(slJ), p_JtoA), dim = dim(slJ))
                sa_denom <- slA + n_JtoA
                sa_mask <- sa_denom > 0
                sa_blend <- (slA * sa_slice + n_JtoA * pmax(p$s_A_min, sJ_slice)) / pmax(sa_denom, 1)
                sa_slice[sa_mask] <- sa_blend[sa_mask]
                J_stayers <- slJ - n_JtoA
                sJ_denom <- J_stayers + n_StoJ
                sJ_mask <- sJ_denom > 0
                sJ_blend <- (J_stayers * sJ_slice + n_StoJ * sS_slice) / pmax(sJ_denom, 1)
                sJ_slice[sJ_mask] <- sJ_blend[sJ_mask]
                slS <- slS - n_StoJ
                slJ <- slJ + n_StoJ - n_JtoA
                slA <- slA + n_JtoA
                ds_base <- (p$delta_s_base / 12) * (precip_annual / 2500) * (mean(relhum, na.rm = TRUE) / 85)
                noise_sd <- p$sigma * 0.5 / sqrt(12)
                delta_s_total_S <- delta_s_total_S + array(ds_base + rnorm(xDim * yDim, 0, noise_sd), dim = c(xDim, yDim))
                delta_s_total_J <- delta_s_total_J + array(ds_base + rnorm(xDim * yDim, 0, noise_sd), dim = c(xDim, yDim))
                delta_s_total_A <- delta_s_total_A + array(ds_base + rnorm(xDim * yDim, 0, noise_sd), dim = c(xDim, yDim))
            }
            abundanceS[, , zi, t, sp] <- slS
            abundanceJ[, , zi, t, sp] <- slJ
            abundanceA[, , zi, t, sp] <- slA
            delta_s_total_A[fruited[, , zi, sp]] <- delta_s_total_A[fruited[, , zi, sp]] * p$cost_repro
            sS_new <- pmin(p$s_S_max, pmax(p$s_S_min, sS_slice + delta_s_total_S))
            sJ_new <- pmin(p$s_J_max, pmax(p$s_J_min, sJ_slice + delta_s_total_J))
            sA_new <- pmin(p$s_A_max, pmax(p$s_A_min, sa_slice + delta_s_total_A))
            has_S <- abundanceS[, , zi, t, sp] > 0
            has_J <- abundanceJ[, , zi, t, sp] > 0
            has_A <- abundanceA[, , zi, t, sp] > 0
            sS_slice[has_S] <- sS_new[has_S]
            sJ_slice[has_J] <- sJ_new[has_J]
            sa_slice[has_A] <- sA_new[has_A]
            size_S[, , zi, t, sp] <- sS_slice
            size_J[, , zi, t, sp] <- sJ_slice
            size_A[, , zi, t, sp] <- sa_slice
        }
    }
    list(
        S = abundanceS, J = abundanceJ, A = abundanceA,
        size_S = size_S, size_J = size_J, size_A = size_A
    )
}

runcolonization_spatial <- function(site, niches, canopy_grid, microenv,
                                    timesteps = 50L, resolution = 10L, carCap = 1L,
                                    maxDisp = 5L, stochastic = FALSE,
                                    Visualize = FALSE, sleeptime = 0.2,
                                    visualize_dispersion = FALSE, spinup = 5L,
                                    parameters = NULL, forestparams = NULL,
                                    allsites = FALSE, train_frac = 0.70,
                                    seed = 42L, clim_cache = NULL, clim_cache_spatial = NULL) {
    if (is.null(parameters)) stop("parameters list required")
    set.seed(seed)
    site_obs_all <- if (allsites) niches else niches[niches$Area_or_Site == site, ]
    train_idx <- unlist(lapply(
        split(
            seq_len(nrow(site_obs_all)),
            site_obs_all$FinalID
        ),
        function(idx) sample(idx, max(1L, round(length(idx) * train_frac)))
    ))
    niches_train <- site_obs_all[train_idx, ]
    niches_val <- site_obs_all[-train_idx, ]
    state <- init_colonization(site, niches_train, canopy_grid, microenv,
        resolution, carCap, maxDisp,
        params = parameters,
        forestparams = forestparams, allsites = allsites,
        clim_cache = clim_cache
    )
    if (is.null(clim_cache_spatial)) clim_cache_spatial <- build_clim_cache_spatial(microenv)
    state$clim_spatial_by_height <- clim_cache_spatial$clim_by_height
    state$clim_spatial_month_by_height <- clim_cache_spatial$clim_month_by_height

    xDim <- state$xDim
    yDim <- state$yDim
    zDim <- state$zDim
    n_species <- state$n_species
    abundanceS <- array(0L, dim = c(xDim, yDim, zDim, timesteps, n_species))
    abundanceJ <- array(0L, dim = c(xDim, yDim, zDim, timesteps, n_species))
    abundanceA <- array(0L, dim = c(xDim, yDim, zDim, timesteps, n_species))
    size_S <- array(parameters$s_S_min, dim = c(xDim, yDim, zDim, timesteps, n_species))
    size_J <- array(parameters$s_J_min, dim = c(xDim, yDim, zDim, timesteps, n_species))
    size_A <- array(parameters$s_A_min, dim = c(xDim, yDim, zDim, timesteps, n_species))
    totalS <- numeric(timesteps)
    totalJ <- numeric(timesteps)
    totalA <- numeric(timesteps)

    sp <- run_spinup_spatial(state,
        n_gens = spinup, Visualize = Visualize,
        carCap = carCap, sleeptime = sleeptime,
        visualize_dispersion = visualize_dispersion
    )
    abundanceS[, , , 1, ] <- sp$S
    abundanceJ[, , , 1, ] <- sp$J
    abundanceA[, , , 1, ] <- sp$A
    size_S[, , , 1, ] <- sp$size_S
    size_J[, , , 1, ] <- sp$size_J
    size_A[, , , 1, ] <- sp$size_A
    totalS[1] <- sum(abundanceS[, , , 1])
    totalJ[1] <- sum(abundanceJ[, , , 1])
    totalA[1] <- sum(abundanceA[, , , 1])
    fruited <- array(FALSE, dim = c(xDim, yDim, zDim, n_species))
    Disp <- sp$last_disp

    for (t in seq_len(timesteps - 1L)) {
        Disp <- .run_pass1_disperse_spatial(state, abundanceA, size_A, t)
        pass2 <- .run_pass2_establish_spatial(state, abundanceS, abundanceJ, abundanceA, size_S, Disp, t)
        abundanceS <- pass2$S
        size_S <- pass2$size_S
        res <- .run_pass3_survive_grow_spatial(
            state, abundanceS, abundanceJ, abundanceA,
            size_S, size_J, size_A, fruited, t
        )
        abundanceS <- res$S
        abundanceJ <- res$J
        abundanceA <- res$A
        size_S <- res$size_S
        size_J <- res$size_J
        size_A <- res$size_A
        n_new <- abundanceS[, , , t + 1, ]
        n_carry <- abundanceS[, , , t, ]
        tot_n <- n_new + n_carry
        blended <- (n_new * size_S[, , , t + 1, ] + n_carry * size_S[, , , t, ]) / pmax(tot_n, 1)
        blended[tot_n == 0] <- parameters$s_S_min
        size_S[, , , t + 1, ] <- blended
        abundanceS[, , , t + 1, ] <- abundanceS[, , , t + 1, ] + abundanceS[, , , t, ]
        abundanceJ[, , , t + 1, ] <- abundanceJ[, , , t + 1, ] + abundanceJ[, , , t, ]
        abundanceA[, , , t + 1, ] <- abundanceA[, , , t + 1, ] + abundanceA[, , , t, ]
        size_J[, , , t + 1, ] <- size_J[, , , t, ]
        size_A[, , , t + 1, ] <- size_A[, , , t, ]
        totalS[t + 1] <- sum(abundanceS[, , , t + 1])
        totalJ[t + 1] <- sum(abundanceJ[, , , t + 1])
        totalA[t + 1] <- sum(abundanceA[, , , t + 1])
    }

    list(
        abundanceS = abundanceS, abundanceJ = abundanceJ, abundanceA = abundanceA,
        totalabundanceS = totalS, totalabundanceJ = totalJ, totalabundanceA = totalA,
        state = state, obs_train = niches_train, obs_val = niches_val
    )
}

main <- function() {
    opts <- parse_args()
    log_msg("Starting test.R")
    log_msg("Loading paths from ", paths_file)
    source_engine_files()
    log_msg("Loaded engine files and helper scripts")

    microenv_path <- opts$microenv
    if (is.null(microenv_path)) microenv_path <- find_microenv()
    if (is.null(microenv_path)) stop("Cannot find La Elenita microenv RDS. Pass --microenv=/path/to/file.rds")
    log_msg("Using microenv: ", microenv_path)
    microenv <- readRDS(microenv_path)

    log_msg("Inspecting raw spatial variance for height 0.5")
    inspect_height(0.5, microenv, hours = c(1, 50, 100, 150, 200, 250))
    log_msg("Finished raw variance inspection")

    if (!opts$compare) {
        log_msg("Skipping comparison because --compare was not supplied")
        return(invisible(NULL))
    }

    log_msg("Loading parameters from ", file.path(PARAMS_DIR, "realistic_273founders.rds"))
    params <- readRDS(file.path(PARAMS_DIR, "realistic_273founders.rds"))
    if (!is.null(opts$n_founders)) params$n_founders <- opts$n_founders

    height_subset <- NULL
    env_heights <- Sys.getenv("SPATIAL_HEIGHTS", unset = "")
    if (nzchar(env_heights)) {
        height_subset <- as.numeric(strsplit(env_heights, ",")[[1]])
        height_subset <- height_subset[is.finite(height_subset)]
    }

    log_msg("Building averaged climate cache")
    clim_cache_avg <- build_clim_cache(microenv)

    log_msg("Building spatial climate cache")
    clim_cache_spat <- build_clim_cache_spatial(microenv, heights = height_subset)
    log_msg("Finished building climate caches")

    log_msg("Running averaged colonization (timesteps=", opts$timesteps, ", spinup=", opts$spinup, ")")
    res_avg <- bench("averaged", quote({
        log_msg("Launching averaged colonization")
        runcolonization(
            opts$site, niches, canopy_grid, microenv,
            timesteps = opts$timesteps, spinup = opts$spinup,
            parameters = params, clim_cache = clim_cache_avg,
            seed = opts$seed
        )
    }))

    log_msg("Running spatial colonization (timesteps=", opts$timesteps, ", spinup=", opts$spinup, ")")
    res_spat <- bench("spatial", quote({
        log_msg("Launching spatial colonization")
        runcolonization_spatial(
            opts$site, niches, canopy_grid, microenv,
            timesteps = opts$timesteps, spinup = opts$spinup,
            parameters = params, clim_cache = clim_cache_avg,
            clim_cache_spatial = clim_cache_spat,
            seed = opts$seed
        )
    }))

    if (!dir.exists(opts$output)) dir.create(opts$output, recursive = TRUE)
    saveRDS(list(avg = res_avg, spatial = res_spat), file = file.path(opts$output, "la_elenita_spatial_vs_avg.rds"))
    log_msg("Saved comparison results to ", file.path(opts$output, "la_elenita_spatial_vs_avg.rds"))
}

main()
