# Repository entry point for the complete analysis reported in the paper.
# Run from the repository root with:
#   source("R/01_run_mobo_final.R")
# or:
#   Rscript R/01_run_mobo_final.R
#
# The full experiment is computationally intensive but supports checkpoints.

# ============================================================
# COMPLETE BO PIPELINE driven by Julia simulation outputs only
# Multi-objective BO (3 objs) + MNL surrogate + TS/SCAL/EHVI/EPI
# ============================================================

options(stringsAsFactors = FALSE)

# ===============================
# Required packages
# ===============================
req_pkgs <- c("dplyr","tidyr","lhs","nnet","JuliaCall")
for (p in req_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Please install required package: ", p)
}
library(dplyr)
library(tidyr)
library(lhs)
library(nnet)
library(JuliaCall)

# ===============================
# Step 0: Portable Julia/project setup
# ===============================
source(file.path("R", "project_setup.R"))

setup_result <- julia_call("setup_simulator", boundary_condition = "closed")
print(setup_result)

sleipner_cfg_status <- julia_call("setup_sleipner_reservoir")
print(sleipner_cfg_status)
if (!identical(sleipner_cfg_status$status, "success")) {
  stop("Julia Sleipner reservoir configuration failed: ", sleipner_cfg_status$message)
}

# ===============================
# Utilities
# ===============================

# Mt/yr -> m^3/yr (1 Mt = 1e9 kg)
mt_to_m3yr <- function(rate_mt, rho = 425.0) (as.numeric(rate_mt) * 1e9) / as.numeric(rho)

# fixed well locations, in index space
get_well_locations <- function(nx_bc, ny_bc, scenario = 1L) {
  ix_center <- nx_bc %/% 2
  iy_center <- ny_bc %/% 2
  if (scenario == 1L) {
    return(list(n_wells = 1L, locs = list(c(ix_center, iy_center))))
  }
  if (scenario == 3L) {
    # clustered triangle
    return(list(
      n_wells = 3L,
      locs = list(
        c(ix_center,     iy_center),
        c(ix_center + 2, iy_center),
        c(ix_center,     iy_center + 2)
      )
    ))
  }
  stop("scenario must be 1 or 3")
}

# Build layer split schedule (5y + 5y + 5y = 15 steps)
# Segment 1: fixed bottom layers 2/4/6 with 0.50/0.30/0.20
# Segment 2: 100% into L_mid
# Segment 3: 100% into L_top
build_layer_split <- function(t_year, L_mid, L_top, n_layers) {
  L_mid <- max(1L, min(n_layers, as.integer(L_mid)))
  L_top <- max(1L, min(n_layers, as.integer(L_top)))
  
  q_split <- rep(0, n_layers)
  if (t_year < 5) {
    bottom_layers <- c(2L, 4L, 6L)
    bottom_frac   <- c(0.50, 0.30, 0.20)
    for (k in seq_along(bottom_layers)) {
      ell <- bottom_layers[k]
      if (ell >= 1 && ell <= n_layers) q_split[ell] <- q_split[ell] + bottom_frac[k]
    }
  } else if (t_year < 10) {
    q_split[L_mid] <- 1.0
  } else {
    q_split[L_top] <- 1.0
  }
  q_split
}

# Build injection matrices for Julia run_simulation
# Dimensions per layer: (time, nx_bc, ny_bc)
build_injection_matrices <- function(scenario, rate_mt, L_mid, L_top,
                                     setup_result,
                                     years = 15L,
                                     rho_inj = 425.0,
                                     availability = 0.95) {
  nx <- setup_result$nx
  ny <- setup_result$ny
  n_layers <- setup_result$n_layers
  
  nx_bc <- nx + 2
  ny_bc <- ny + 2
  
  well_info <- get_well_locations(nx_bc, ny_bc, scenario = scenario)
  n_wells <- well_info$n_wells
  locs <- well_info$locs
  
  # decision: per-well annual injection rate
  # total injection scales with number of wells implicitly by applying per-well rate to each well
  rate_m3_yr_per_well <- mt_to_m3yr(rate_mt, rho = rho_inj) * availability
  
  n_times <- years
  injection_matrices <- vector("list", n_layers)
  for (ell in seq_len(n_layers)) {
    injection_matrices[[ell]] <- array(0, dim = c(n_times, nx_bc, ny_bc))
  }
  
  for (ti in seq_len(n_times)) {
    t_year <- ti - 1L  # 0..14
    q_split <- build_layer_split(t_year, L_mid, L_top, n_layers)
    
    # apply same schedule to each well location
    for (w in seq_len(n_wells)) {
      loc <- locs[[w]]
      for (ell in seq_len(n_layers)) {
        if (q_split[ell] > 0) {
          injection_matrices[[ell]][ti, loc[1], loc[2]] <-
            injection_matrices[[ell]][ti, loc[1], loc[2]] + q_split[ell] * rate_m3_yr_per_well
        }
      }
    }
  }
  
  # All numerical layers are retained as candidate reservoir intervals.
  # Layers 8 and 9 are NOT zeroed. Selecting L_mid = 8 or L_top = 9
  # therefore produces active injection during the corresponding phase.
  # This is a hypothetical multi-interval optimization; it should not be
  # described as the historical Sleipner well completion.
  
  list(injection_matrices = injection_matrices,
       nx_bc = nx_bc, ny_bc = ny_bc,
       n_layers = n_layers,
       n_times = n_times,
       rate_m3_yr_per_well = rate_m3_yr_per_well,
       n_wells = well_info$n_wells)
}

# Internal check: the final two numerical layers are active when selected.
stopifnot(
  build_layer_split(5, L_mid = 8L, L_top = 9L, n_layers = 9L)[8] == 1,
  build_layer_split(10, L_mid = 8L, L_top = 9L, n_layers = 9L)[9] == 1
)

# IMPORTANT:
# This cache stores only the small quantities required by the objective:
# the retained-volume time series, total injected volume, and number of wells.
.julia_cache_env <- new.env(parent = emptyenv())
.julia_sim_counter <- 0L

cache_key <- function(scenario, rate_mt, L_mid, L_top) {
  paste0(
    "s=", scenario,
    "|r=", sprintf("%.6f", rate_mt),
    "|m=", L_mid,
    "|t=", L_top
  )
}

run_julia_sim_cached <- function(scenario, rate_mt, L_mid, L_top, cfg) {
  key <- cache_key(scenario, rate_mt, L_mid, L_top)

  if (exists(key, envir = .julia_cache_env, inherits = FALSE)) {
    return(get(key, envir = .julia_cache_env, inherits = FALSE))
  }

  inj <- build_injection_matrices(
    scenario = scenario,
    rate_mt = rate_mt,
    L_mid = L_mid,
    L_top = L_top,
    setup_result = setup_result,
    years = cfg$sim_years,
    rho_inj = cfg$rho_inj,
    availability = cfg$a_bar
  )

  sim_result <- julia_call(
    "run_simulation",
    start_time = 0.0,
    end_time = as.numeric(cfg$sim_years),
    time_step = 1.0,
    injection_rate_matrices = inj$injection_matrices,
    verbose = FALSE
  )

  if (!is.list(sim_result) ||
      is.null(sim_result$status) ||
      sim_result$status != "success") {
    out <- list(
      status = "failed",
      total_co2_volumes = numeric(0),
      injected_total = NA_real_,
      n_wells = inj$n_wells
    )
  } else {
    injected_total <- sum(
      vapply(inj$injection_matrices, sum, numeric(1))
    )

    out <- list(
      status = "success",
      total_co2_volumes = as.numeric(sim_result$total_co2_volumes),
      injected_total = as.numeric(injected_total),
      n_wells = as.integer(inj$n_wells)
    )
  }

  # Store only the compact result, never the large injection matrices.
  assign(key, out, envir = .julia_cache_env)

  # Drop large temporary R objects immediately.
  rm(inj, sim_result)
  .julia_sim_counter <<- .julia_sim_counter + 1L

  # Periodic garbage collection in both R and Julia.
  if (.julia_sim_counter %% 20L == 0L) {
    invisible(gc(verbose = FALSE))
    try(julia_command("GC.gc()"), silent = TRUE)
  }

  out
}

clear_compact_simulation_cache <- function() {
  rm(list = ls(envir = .julia_cache_env, all.names = TRUE),
     envir = .julia_cache_env)
  invisible(gc(verbose = FALSE))
  try(julia_command("GC.gc()"), silent = TRUE)
  invisible(NULL)
}

# ===============================
# cfg (BO + objective parameters)
# ===============================
cfg <- list(
  # Julia / injection constants
  sim_years = 15L,
  rho_inj   = 425.0,   # kg/m^3, used for Mt/yr <-> m^3/yr conversion
  a_bar     = 0.95,

  # Objective 3: transparent representative European CCS-chain costs.
  # Currency and price basis: 2024 EUR per metric tonne of injected CO2.
  # Values are midpoints of the JRC 2024 ranges:
  #   capture: 40-90 EUR/tCO2 -> 65 EUR/tCO2
  #   transport: 2-30 EUR/tCO2 -> 16 EUR/tCO2
  #   storage: 5-35 EUR/tCO2 -> 20 EUR/tCO2
  # Source: European Commission Joint Research Centre, Clean Energy
  # Technology Observatory: CCUS in the European Union, 2024, JRC139285.
  cost_currency = "EUR",
  cost_price_year = 2024L,
  c_cap      = 65.0,   # EUR/tCO2
  c_trans    = 16.0,   # EUR/tCO2
  c_storage  = 20.0,   # EUR/tCO2

  # BO controls
  bins     = 5L,
  mc_samp  = 20L,
  ehvi_k   = 0.0,
  ucb_k    = 2.0,
  n_init   = 20L,
  n_iter   = 80L,
  n_cand   = 30L
)

# ===============================
# Objective function driven by Julia simulation (three objectives)
# f1: maximize time-integrated retained CO2 volume, m^3 yr
# f2: minimize final-time unretained CO2 volume, m^3
# f3: minimize representative CCS-chain cost, 2024 EUR
# ===============================
co2_obj_julia <- function(x, cfg) {
  scenario <- as.integer(x[1])
  rate_mt <- as.numeric(x[2])
  L_mid <- as.integer(x[3])
  L_top <- as.integer(x[4])

  if (!(scenario %in% c(1L, 3L))) scenario <- 1L
  n_layers <- as.integer(setup_result$n_layers)
  L_mid <- max(1L, min(n_layers, L_mid))
  L_top <- max(1L, min(n_layers, L_top))

  rr <- run_julia_sim_cached(
    scenario,
    rate_mt,
    L_mid,
    L_top,
    cfg
  )

  if (!is.list(rr) ||
      is.null(rr$status) ||
      rr$status != "success") {
    return(c(f1 = -1e30, f2 = 1e30, f3 = 1e30))
  }

  vol_vec <- as.numeric(rr$total_co2_volumes)
  vol_vec[!is.finite(vol_vec)] <- 0
  stored_final <- tail(vol_vec, 1)

  # f1: time-integrated retained CO2 volume (dt = 1 year)
  f1 <- sum(vol_vec)

  # f2: final-time unretained CO2 volume
  injected_total_m3 <- as.numeric(rr$injected_total)
  f2 <- max(0, injected_total_m3 - stored_final)

  # Convert injected reservoir volume to metric tonnes of CO2.
  injected_total_tonnes <- injected_total_m3 * cfg$rho_inj / 1000.0

  # f3: representative 2024-EUR CCS-chain cost.
  # Fixed and well costs are not added separately because the JRC storage
  # range is an aggregate storage cost; this avoids double counting.
  unit_chain_cost_eur_per_t <- (
    cfg$c_cap +
    cfg$c_trans +
    cfg$c_storage
  )
  f3 <- unit_chain_cost_eur_per_t * injected_total_tonnes

  c(f1 = f1, f2 = f2, f3 = f3)
}

# ============================================================
# Transform to maximization for HV/BO (3D)
# f1 max, f2/f3 min => (f1, -f2, -f3)
# ============================================================
to_max3 <- function(y) {
  y <- as.numeric(y)
  c(y[1], -y[2], -y[3])
}
to_max3_mat <- function(Y) {
  Y <- as.matrix(Y)
  cbind(Y[,1], -Y[,2], -Y[,3])
}

# ============================================================
# Nondominated set (MAX space)
# ============================================================
nondominated <- function(Y) {
  Ym <- as.matrix(Y)
  n <- nrow(Ym)
  if (n == 0) return(Ym)
  is_dom <- rep(FALSE, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i != j &&
          all(Ym[j, ] >= Ym[i, ], na.rm = TRUE) &&
          any(Ym[j, ] >  Ym[i, ], na.rm = TRUE)) {
        is_dom[i] <- TRUE
        break
      }
    }
  }
  Ym[!is_dom, , drop = FALSE]
}

# ============================================================
# Normalize to [0,1] per dim using fixed bounds
# ============================================================
normalize_to_unit <- function(Y, lo, hi) {
  Y <- as.matrix(Y)
  out <- Y
  for (j in seq_len(ncol(Y))) {
    denom <- hi[j] - lo[j]
    if (!is.finite(denom) || denom <= 0) {
      out[,j] <- 0
    } else {
      out[,j] <- (Y[,j] - lo[j]) / denom
      out[,j] <- pmin(pmax(out[,j], 0), 1)
    }
  }
  out
}

make_fixed_bounds_from_Y0 <- function(Y0_raw, pad = 0.10) {
  Y0_max <- to_max3_mat(Y0_raw)
  fmin <- apply(Y0_max, 2, min, na.rm = TRUE)
  fmax <- apply(Y0_max, 2, max, na.rm = TRUE)
  span <- fmax - fmin
  span[!is.finite(span) | span == 0] <- pmax(abs(fmin[!is.finite(span) | span == 0]), 1.0)
  lo <- fmin - pad * span
  hi <- fmax + pad * span
  list(lo = lo, hi = hi)
}

# ============================================================
# Hypervolume in normalized MAX space
# FAST EXACT 3-D implementation
# ============================================================
.hv2d_origin_exact_fast <- function(P) {
  P <- as.matrix(P)
  if (nrow(P) == 0L) return(0)

  P <- P[
    is.finite(P[, 1]) & is.finite(P[, 2]) &
      P[, 1] > 0 & P[, 2] > 0,
    ,
    drop = FALSE
  ]
  if (nrow(P) == 0L) return(0)

  y_breaks <- sort(unique(c(0, P[, 1])))
  area <- 0

  for (i in 2:length(y_breaks)) {
    y_left <- y_breaks[i - 1L]
    y_right <- y_breaks[i]
    active <- P[, 1] >= y_right

    if (any(active)) {
      area <- area +
        (y_right - y_left) *
        max(P[active, 2], na.rm = TRUE)
    }
  }

  as.numeric(area)
}

compute_hv <- function(Y, ref, n_mc = 15000) {
  P <- as.matrix(Y)
  ref <- as.numeric(ref)

  if (nrow(P) == 0L) return(0)
  if (ncol(P) != 3L) {
    stop("The exact fast hypervolume function requires exactly 3 objectives.")
  }
  if (length(ref) != 3L || any(!is.finite(ref))) return(0)

  valid <- apply(P, 1, function(z) all(is.finite(z)))
  P <- P[valid, , drop = FALSE]
  if (nrow(P) == 0L) return(0)

  P <- unique(nondominated(P))

  # Translate the reference point to the origin.
  Q <- sweep(P, 2, ref, "-")
  Q <- Q[apply(Q, 1, function(z) all(z > 0)), , drop = FALSE]
  if (nrow(Q) == 0L) return(0)

  x_breaks <- sort(unique(c(0, Q[, 1])))
  volume <- 0

  for (i in 2:length(x_breaks)) {
    x_left <- x_breaks[i - 1L]
    x_right <- x_breaks[i]
    active <- Q[, 1] >= x_right

    if (any(active)) {
      volume <- volume +
        (x_right - x_left) *
        .hv2d_origin_exact_fast(Q[active, 2:3, drop = FALSE])
    }
  }

  as.numeric(volume)
}

# Quick internal checks.
stopifnot(
  abs(
    compute_hv(
      matrix(c(1, 1, 1), nrow = 1),
      ref = c(0, 0, 0)
    ) - 1
  ) < 1e-10
)

# ============================================================
# Multinomial-logit surrogate (3 objectives)
# X = scenario, rate_mt, L_mid, L_top
# ============================================================

.make_quantile_breaks <- function(y, B) {
  y <- as.numeric(y); y <- y[is.finite(y)]
  if (length(y) < 2) return(c(-Inf, Inf))
  probs <- seq(0, 1, length.out = B + 1)
  qs <- as.numeric(stats::quantile(y, probs = probs, na.rm = TRUE, type = 7))
  for (k in 2:length(qs)) if (!is.finite(qs[k]) || qs[k] <= qs[k-1]) qs[k] <- qs[k-1] + 1e-12
  qs[1] <- -Inf; qs[length(qs)] <- Inf
  qs
}

create_surrogate_mnl <- function(X, Y, cfg = list()) {
  X <- as.data.frame(X)
  Y <- as.matrix(Y)
  if (nrow(X) != nrow(Y)) stop("create_surrogate_mnl: nrow(X)!=nrow(Y)")
  if (ncol(Y) < 3) stop("create_surrogate_mnl: expects 3 objectives")
  
  stopifnot(all(c("scenario","rate_mt","L_mid","L_top") %in% names(X)))
  
  X$scenario <- factor(X$scenario, levels = c(1,3))
  X$L_mid <- factor(X$L_mid, levels = as.character(1:setup_result$n_layers))
  X$L_top <- factor(X$L_top, levels = as.character(1:setup_result$n_layers))
  
  B <- if (!is.null(cfg$bins)) as.integer(cfg$bins) else 5L
  B <- max(2L, B)
  
  breaks_list <- vector("list", 3L)
  bins_mat    <- matrix(1L, nrow = nrow(Y), ncol = 3L)
  
  for (i in 1:3) {
    breaks_list[[i]] <- .make_quantile_breaks(Y[, i], B = B)
    br <- sort(as.numeric(breaks_list[[i]]))
    br[1] <- -Inf; br[length(br)] <- Inf
    eps <- 1e-12
    for (k in 2:length(br)) if (br[k] <= br[k-1]) br[k] <- br[k-1] + eps
    
    b <- findInterval(as.numeric(Y[, i]), vec = br, rightmost.closed = TRUE, all.inside = TRUE)
    b <- pmax(1L, pmin(length(br) - 1L, b))
    bins_mat[, i] <- as.integer(b)
    
    breaks_list[[i]] <- br
  }
  
  class_lbl <- apply(bins_mat, 1, function(b) paste0(b[1], "_", b[2], "_", b[3]))
  class_fct <- factor(class_lbl)
  
  dfY <- as.data.frame(Y); names(dfY)[1:3] <- c("f1","f2","f3")
  dfY$class <- class_lbl
  
  centroids_df <- dfY %>%
    dplyr::group_by(.data$class) %>%
    dplyr::summarise(
      mu1 = mean(.data$f1, na.rm = TRUE),
      mu2 = mean(.data$f2, na.rm = TRUE),
      mu3 = mean(.data$f3, na.rm = TRUE),
      n   = dplyr::n(),
      .groups = "drop"
    )
  
  if (length(levels(class_fct)) < 2) {
    sur <- list(X = X, Y = Y, B = B,
                breaks_list = breaks_list,
                class_levels = levels(class_fct),
                centroid_map = centroids_df,
                fit = NULL)
    class(sur) <- "sur_mnl_quantile"
    return(sur)
  }
  
  df_train <- data.frame(class = class_fct, X)
  fit <- nnet::multinom(class ~ scenario + rate_mt + L_mid + L_top, data = df_train, trace = FALSE)
  
  sur <- list(X = X, Y = Y, B = B,
              breaks_list = breaks_list,
              class_levels = levels(class_fct),
              centroid_map = centroids_df,
              fit = fit)
  class(sur) <- "sur_mnl_quantile"
  sur
}

predict_class_probs_mnl <- function(sur, x_cand) {
  if (!inherits(sur, "sur_mnl_quantile")) stop("sur wrong class")
  
  K <- length(sur$class_levels)
  if (K <= 0) stop("sur$class_levels is empty")
  
  # If fit missing/unusable: uniform
  if (is.null(sur$fit)) {
    probs <- rep(1 / K, K)
    names(probs) <- sur$class_levels
    return(probs)
  }
  
  x <- as.data.frame(x_cand)
  x$scenario <- factor(x$scenario, levels = c(1,3))
  x$L_mid    <- factor(x$L_mid, levels = as.character(1:setup_result$n_layers))
  x$L_top    <- factor(x$L_top, levels = as.character(1:setup_result$n_layers))
  
  pr <- tryCatch(
    stats::predict(sur$fit, newdata = x, type = "probs"),
    error = function(e) NULL
  )
  
  # If predict failed -> uniform
  if (is.null(pr)) {
    probs <- rep(1 / K, K)
    names(probs) <- sur$class_levels
    return(probs)
  }
  
  # Helper: normalize output to length-K vector aligned to class_levels
  normalize_out <- function(v, nm = NULL) {
    out <- rep(0, K)
    names(out) <- sur$class_levels
    
    if (!is.null(nm)) {
      common <- intersect(nm, sur$class_levels)
      out[common] <- as.numeric(v[match(common, nm)])
    } else {
      out[seq_len(min(K, length(v)))] <- as.numeric(v[seq_len(min(K, length(v)))])
    }
    
    out[!is.finite(out)] <- 0
    out <- pmax(out, 0)
    if (sum(out) <= 0) out <- rep(1 / K, K)
    out / sum(out)
  }
  
  # Case 1: matrix/data.frame
  if (!is.null(dim(pr))) {
    v <- as.numeric(pr[1, , drop = TRUE])
    nm <- colnames(pr)
    return(normalize_out(v, nm))
  }
  
  # Case 2: vector (may be named)
  if (is.atomic(pr)) {
    nm <- names(pr)
    
    # Special case: multinom sometimes returns a single prob in binary case
    if (length(pr) == 1 && K == 2) {
      p2 <- as.numeric(pr[1])
      if (!is.finite(p2)) p2 <- 0.5
      p2 <- max(min(p2, 1), 0)
      v <- c(1 - p2, p2)
      names(v) <- sur$class_levels
      return(v / sum(v))
    }
    
    return(normalize_out(pr, nm))
  }
  
  # Fallback: uniform
  probs <- rep(1 / K, K)
  names(probs) <- sur$class_levels
  probs
}

sample_Y_for_candidate_mnl <- function(sur, x_cand, cfg = list()) {
  probs <- predict_class_probs_mnl(sur, x_cand)
  cls <- sample(names(probs), size = 1, prob = probs)
  row <- sur$centroid_map[sur$centroid_map$class == cls, , drop = FALSE]
  if (nrow(row) == 0) {
    cls2 <- names(sort(probs, decreasing = TRUE))[1]
    row <- sur$centroid_map[sur$centroid_map$class == cls2, , drop = FALSE]
  }
  as.numeric(row[1, c("mu1","mu2","mu3")])
}

create_surrogate <- function(X, Y, cfg = list()) create_surrogate_mnl(X, Y, cfg)
sample_Y_for_candidate <- function(sur, x_cand, cfg = list()) sample_Y_for_candidate_mnl(sur, x_cand, cfg)

# ============================================================
# Acquisition functions (3 objectives)
# ============================================================

acq_scal_ucb <- function(sur, Xcand_feat, cfg, bounds,
                         w = c(1,1,1), n_samp = NULL, kappa = NULL) {
  if (is.null(n_samp)) n_samp <- min(200L, cfg$mc_samp)
  if (is.null(kappa))  kappa  <- cfg$ucb_k
  w <- as.numeric(w)
  
  score_one <- function(j) {
    util <- numeric(n_samp)
    for (s in seq_len(n_samp)) {
      yraw <- sample_Y_for_candidate(sur, Xcand_feat[j, , drop = FALSE], cfg)
      ymax <- to_max3(yraw)
      yn   <- as.numeric(normalize_to_unit(matrix(ymax, 1), bounds$lo, bounds$hi))
      util[s] <- sum(w * yn)
    }
    m <- mean(util, na.rm = TRUE)
    sdv <- stats::sd(util, na.rm = TRUE)
    if (!is.finite(m)) m <- -Inf
    if (!is.finite(sdv)) sdv <- 0
    m + kappa * sdv
  }
  
  vapply(seq_len(nrow(Xcand_feat)), score_one, numeric(1))
}

acq_ts_hv <- function(sur, Xcand_feat, front0_norm, ref_norm, cfg, bounds) {
  hv0 <- compute_hv(front0_norm, ref = ref_norm, n_mc = 1000)
  score_one <- function(j) {
    yraw <- sample_Y_for_candidate(sur, Xcand_feat[j, , drop = FALSE], cfg)
    pt   <- to_max3(yraw)
    pt_n <- normalize_to_unit(matrix(pt, 1), bounds$lo, bounds$hi)
    hv1 <- compute_hv(nondominated(rbind(front0_norm, pt_n)), ref = ref_norm, n_mc = 1000)
    hv1 - hv0
  }
  vapply(seq_len(nrow(Xcand_feat)), score_one, numeric(1))
}

# Standard Expected Preference Improvement (EPI):
# expected positive improvement in the scalarized normalized utility.
# No additional UCB-style standard-deviation term is added.
acq_epi <- function(sur, Xcand_feat, history_Y_raw, cfg, bounds,
                    w = c(1,1,1), n_samp = NULL) {
  if (is.null(n_samp)) n_samp <- min(200L, cfg$mc_samp)
  w <- as.numeric(w)

  if (is.null(history_Y_raw) || nrow(as.matrix(history_Y_raw)) == 0) {
    u_best <- -Inf
  } else {
    Ymax <- to_max3_mat(history_Y_raw)
    Yn <- normalize_to_unit(Ymax, bounds$lo, bounds$hi)
    u_best <- max(as.numeric(Yn %*% w), na.rm = TRUE)
    if (!is.finite(u_best)) u_best <- -Inf
  }

  score_one <- function(j) {
    imp <- numeric(n_samp)

    for (s in seq_len(n_samp)) {
      yraw <- sample_Y_for_candidate(
        sur,
        Xcand_feat[j, , drop = FALSE],
        cfg
      )
      pt <- to_max3(yraw)
      pt_n <- as.numeric(
        normalize_to_unit(
          matrix(pt, 1),
          bounds$lo,
          bounds$hi
        )
      )
      u <- sum(w * pt_n)
      imp[s] <- max(0, u - u_best)
    }

    value <- mean(imp, na.rm = TRUE)
    if (!is.finite(value)) value <- -Inf
    value
  }

  vapply(seq_len(nrow(Xcand_feat)), score_one, numeric(1))
}

# Standard Monte Carlo Expected Hypervolume Improvement (EHVI):
# the acquisition value is the sample mean of positive hypervolume
# improvements. No additional UCB-style standard-deviation term is added.
acq_ehvi_mc <- function(sur, Xcand_feat, front0_norm, ref_norm, cfg, bounds,
                        n_samp = NULL, hv_mc = 1000) {
  if (is.null(n_samp)) n_samp <- min(200L, cfg$mc_samp)

  hv0 <- compute_hv(
    front0_norm,
    ref = ref_norm,
    n_mc = hv_mc
  )

  score_one <- function(j) {
    imps <- numeric(n_samp)

    for (s in seq_len(n_samp)) {
      yraw <- sample_Y_for_candidate(
        sur,
        Xcand_feat[j, , drop = FALSE],
        cfg
      )
      pt <- to_max3(yraw)
      pt_n <- normalize_to_unit(
        matrix(pt, 1),
        bounds$lo,
        bounds$hi
      )
      hv1 <- compute_hv(
        nondominated(rbind(front0_norm, pt_n)),
        ref = ref_norm,
        n_mc = hv_mc
      )
      imps[s] <- max(0, hv1 - hv0)
    }

    value <- mean(imps, na.rm = TRUE)
    if (!is.finite(value)) value <- -Inf
    value
  }

  vapply(seq_len(nrow(Xcand_feat)), score_one, numeric(1))
}


# ============================================================
# - Julia setup and simulator code
# - 10 paired seeds and 80 sequential iterations are unchanged.
# - Candidate set is 30 for every method and every seed.
# - Acquisition samples are 20 for every method and every seed.
# - EHVI and EPI use their standard sample-mean definitions (no +0.1 SD term).
# - Hypervolume is exact in 3-D instead of repeatedly estimated with
#   thousands of Monte Carlo points.
# Multi-seed BO + common reference front + HV/IGD+ statistics
# ============================================================

# ------------------------------------------------------------
# 1) Reproducible experiment settings
# ------------------------------------------------------------
cfg$n_init  <- 20L
cfg$n_iter  <- 80L
cfg$n_cand  <- 30L  
cfg$bins    <- 5L
cfg$mc_samp <- 20L   
cfg$ucb_k   <- 2.0
cfg$ehvi_k  <- 0.0

# Primary repeated-run experiment.
# Ten paired seeds meet the request for repeated stochastic runs.
# Change to 1:20 if computationally affordable.
SEEDS <- 1:10
METHODS <- c("TS", "SCAL", "EHVI", "EPI")

# For the final paper:
# "pooled" = reference front from the nondominated union of all repeated runs.
# "grid"   = add a dense grid over scenario, rate, L_mid, and L_top.
# Grid mode is more expensive but directly implements Reviewer 3's suggestion.
REFERENCE_MODE <- "grid_only"  # primary IGD+ reference is the grid-only nondominated front
REFERENCE_RATE_GRID_N <- 21L   # 0.80, 0.81, ..., 1.00 when set to 21

# Set TRUE for the final revision after the primary multi-seed run succeeds.
RUN_HYPERPARAMETER_SENSITIVITY <- TRUE
SENSITIVITY_SEEDS <- 1:5
KAPPA_VALUES <- c(0.5, 1.0, 2.0, 4.0)
EHVI_MC_VALUES <- c(10L, 20L, 40L)

# Output folder is created inside the current Julia project directory.
OUTPUT_DIR <- RESULT_ROOT
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "runs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "reference"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "sensitivity"), recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(OUTPUT_DIR, "physical_cache"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "surrogate_validation"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "simulator_validation"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "cost_sensitivity"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "text_for_manuscript"), recursive = TRUE, showWarnings = FALSE)

# Add a persistent disk layer around the existing compact simulator cache.
# This allows the independent reference-grid evaluation to resume without
# rerunning completed physical simulations.
.run_julia_sim_cached_compact_memory <- run_julia_sim_cached
PHYSICAL_CACHE_DIR <- file.path(OUTPUT_DIR, "physical_cache")

run_julia_sim_cached <- function(scenario, rate_mt, L_mid, L_top, cfg) {
  key <- cache_key(scenario, rate_mt, L_mid, L_top)
  disk_file <- file.path(PHYSICAL_CACHE_DIR, paste0(key, ".rds"))

  if (exists(key, envir = .julia_cache_env, inherits = FALSE)) {
    return(get(key, envir = .julia_cache_env, inherits = FALSE))
  }

  if (file.exists(disk_file)) {
    out <- readRDS(disk_file)
    assign(key, out, envir = .julia_cache_env)
    return(out)
  }

  out <- .run_julia_sim_cached_compact_memory(
    scenario,
    rate_mt,
    L_mid,
    L_top,
    cfg
  )
  saveRDS(out, disk_file)

  # Keep the in-memory cache bounded. Older results remain available on disk.
  cache_names <- ls(envir = .julia_cache_env, all.names = TRUE)
  if (length(cache_names) > 250L) {
    remove_names <- head(cache_names, length(cache_names) - 150L)
    rm(list = remove_names, envir = .julia_cache_env)
    invisible(gc(verbose = FALSE))
    try(julia_command("GC.gc()"), silent = TRUE)
  }

  out
}

# Required only for the new statistical figures.
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Please install required package: ggplot2")
}
library(ggplot2)

# Save all implementation settings explicitly for reproducibility.
implementation_settings <- data.frame(
  setting = c(
    "simulation_years",
    "initial_design_size",
    "sequential_iterations",
    "total_evaluations_per_run",
    "candidate_set_size",
    "MNL_bins_B",
    "acquisition_MC_samples",
    "SCAL_UCB_kappa",
    "EPI_weights",
    "SCAL_UCB_weights",
    "normalized_HV_reference_point",
    "number_of_primary_seeds",
    "primary_seeds",
    "primary_IGDplus_reference",
    "supplementary_IGDplus_reference",
    "reference_rate_grid_size",
    "injectable_numerical_layers",
    "layers_8_and_9_status",
    "cost_currency",
    "cost_price_year",
    "cost_basis",
    "EHVI_definition",
    "EHVI_extra_exploration_coefficient",
    "EPI_definition",
    "EPI_extra_exploration_coefficient"
  ),
  value = c(
    as.character(cfg$sim_years),
    as.character(cfg$n_init),
    as.character(cfg$n_iter),
    as.character(cfg$n_init + cfg$n_iter),
    as.character(cfg$n_cand),
    as.character(cfg$bins),
    as.character(cfg$mc_samp),
    as.character(cfg$ucb_k),
    "1,1,1 (equal normalized-objective weights)",
    "1,1,1 (equal normalized-objective weights)",
    "-0.10,-0.10,-0.10",
    as.character(length(SEEDS)),
    paste(SEEDS, collapse = ","),
    "nondominated front of the 3402-design dense grid only",
    "nondominated union of grid and all BO evaluations",
    as.character(REFERENCE_RATE_GRID_N),
    paste(seq_len(as.integer(setup_result$n_layers)), collapse = ","),
    "active candidate reservoir intervals; not forcibly zeroed",
    cfg$cost_currency,
    as.character(cfg$cost_price_year),
    "representative midpoint of JRC 2024 EUR/tCO2 ranges",
    "Monte Carlo mean of nonnegative hypervolume improvement",
    "0",
    "Monte Carlo mean of nonnegative scalarized preference improvement",
    "0"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  implementation_settings,
  file.path(OUTPUT_DIR, "tables", "implementation_and_hyperparameter_settings.csv"),
  row.names = FALSE
)

# Sleipner interpretation used in the manuscript and response letter.
sleipner_layer_interpretation <- paste(
  "Published Sleipner studies describe injection near the base of the Utsira",
  "Formation followed by buoyant ascent and ponding beneath intraformational",
  "mudstone barriers, producing nine seismically imaged plume accumulations.",
  "The nine numerical layers in this optimization are therefore treated as",
  "candidate reservoir subintervals in a hypothetical multi-interval control",
  "study; direct injection into layers 8 and 9 is not presented as a",
  "reconstruction of the historical Sleipner well completion, and the",
  "overlying caprock remains represented by the closed/sealing boundary."
)
writeLines(
  sleipner_layer_interpretation,
  file.path(OUTPUT_DIR, "text_for_manuscript", "sleipner_layer_interpretation.txt")
)

# ------------------------------------------------------------
# 2) Shared initial design for a paired comparison
# ------------------------------------------------------------
make_initial_history <- function(cfg, seed, debug = TRUE) {
  set.seed(as.integer(seed))
  n_layers <- as.integer(setup_result$n_layers)

  X0 <- data.frame(
    scenario = sample(rep(c(1, 3), length.out = cfg$n_init)),
    rate_mt  = runif(cfg$n_init, 0.8, 1.0),
    L_mid    = sample(seq_len(n_layers), cfg$n_init, replace = TRUE),
    L_top    = sample(seq_len(n_layers), cfg$n_init, replace = TRUE),
    stringsAsFactors = FALSE
  )

  Y0 <- do.call(
    rbind,
    lapply(seq_len(nrow(X0)), function(i) {
      co2_obj_julia(as.numeric(X0[i, ]), cfg)
    })
  )
  Y0 <- as.matrix(Y0)
  colnames(Y0) <- c("f1", "f2", "f3")

  if (debug) {
    message(
      sprintf(
        "[seed %d] shared initial design evaluated; objective SD = %s",
        seed,
        paste(signif(apply(Y0, 2, stats::sd), 5), collapse = " | ")
      )
    )
  }

  list(X = X0, Y = Y0)
}

# ------------------------------------------------------------
# 3) Minimally extended BO loop
# This is the original run_bo logic with only three additions:
#   (i) seed is an argument rather than being fixed at 42;
#   (ii) a shared initial_history can be supplied;
#   (iii) method/iteration-specific seeds make interrupted reruns reproducible.
# The Julia calls and all model calculations remain unchanged.
# ------------------------------------------------------------
run_bo <- function(
    cfg,
    method = c("TS", "SCAL", "EHVI", "EPI"),
    seed = 42L,
    initial_history = NULL,
    debug = FALSE,
    checkpoint_file = NULL,
    checkpoint_every = 5L
) {
  method <- match.arg(method)
  seed <- as.integer(seed)
  n_layers <- as.integer(setup_result$n_layers)

  run_signature <- paste(
    cfg$n_init,
    cfg$n_iter,
    cfg$n_cand,
    cfg$mc_samp,
    cfg$bins,
    cfg$ucb_k,
    cfg$ehvi_k,
    sep = "|"
  )

  resume_state <- NULL
  if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) {
    candidate_state <- tryCatch(
      readRDS(checkpoint_file),
      error = function(e) NULL
    )

    if (is.list(candidate_state) &&
        identical(candidate_state$run_signature, run_signature) &&
        identical(as.integer(candidate_state$seed), seed) &&
        identical(as.character(candidate_state$method), method)) {
      resume_state <- candidate_state
      message(
        "Resuming checkpoint: seed=", seed,
        ", method=", method,
        ", completed iteration=", resume_state$completed_iteration
      )
    } else {
      warning(
        "Ignoring incompatible checkpoint: ",
        checkpoint_file
      )
      unlink(checkpoint_file)
    }
  }

  if (!is.null(resume_state)) {
    history <- resume_state$history
    bounds <- resume_state$bounds
    ref_fixed_norm <- resume_state$ref
    hv_trace <- resume_state$hv_trace
    n0 <- as.integer(resume_state$n0)
    start_iteration <- as.integer(resume_state$completed_iteration) + 1L
  } else {
    if (is.null(initial_history)) {
      initial_history <- make_initial_history(cfg, seed, debug = debug)
    }

    history <- list(
      X = as.data.frame(initial_history$X),
      Y = as.matrix(initial_history$Y)
    )
    colnames(history$Y) <- c("f1", "f2", "f3")

    bounds <- make_fixed_bounds_from_Y0(history$Y, pad = 0.10)
    ref_fixed_norm <- rep(-0.10, 3)

    n0 <- nrow(history$X)
    hv_trace <- rep(NA_real_, n0 + cfg$n_iter)

    Y0_norm <- normalize_to_unit(
      to_max3_mat(history$Y),
      bounds$lo,
      bounds$hi
    )
    front0 <- nondominated(Y0_norm)
    hv_trace[n0] <- compute_hv(
      front0,
      ref = ref_fixed_norm,
      n_mc = 15000
    )

    start_iteration <- 1L
  }

  feat_cols <- c("scenario", "rate_mt", "L_mid", "L_top")
  method_offset <- match(method, METHODS)

  if (start_iteration <= cfg$n_iter) {
    for (it in seq.int(start_iteration, cfg$n_iter)) {
      set.seed(seed * 100000L + method_offset * 1000L + it)

      X_feat_hist <- as.matrix(
        history$X[, feat_cols, drop = FALSE]
      )
      sur <- create_surrogate(
        X_feat_hist,
        history$Y,
        cfg
      )

      Xcand <- data.frame(
        scenario = sample(
          c(1, 3),
          cfg$n_cand,
          replace = TRUE
        ),
        rate_mt = runif(
          cfg$n_cand,
          0.8,
          1.0
        ),
        L_mid = sample(
          seq_len(n_layers),
          cfg$n_cand,
          replace = TRUE
        ),
        L_top = sample(
          seq_len(n_layers),
          cfg$n_cand,
          replace = TRUE
        ),
        stringsAsFactors = FALSE
      )
      Xcand_feat <- as.matrix(
        Xcand[, feat_cols, drop = FALSE]
      )

      Yhist_norm <- normalize_to_unit(
        to_max3_mat(history$Y),
        bounds$lo,
        bounds$hi
      )
      front0_norm <- nondominated(Yhist_norm)

      scores <- switch(
        method,
        TS = acq_ts_hv(
          sur,
          Xcand_feat,
          front0_norm,
          ref_fixed_norm,
          cfg,
          bounds
        ),
        SCAL = acq_scal_ucb(
          sur,
          Xcand_feat,
          cfg,
          bounds,
          w = c(1, 1, 1),
          n_samp = cfg$mc_samp,
          kappa = cfg$ucb_k
        ),
        EPI = acq_epi(
          sur,
          Xcand_feat,
          history$Y,
          cfg,
          bounds,
          w = c(1, 1, 1),
          n_samp = cfg$mc_samp
        ),
        EHVI = acq_ehvi_mc(
          sur,
          Xcand_feat,
          front0_norm,
          ref_fixed_norm,
          cfg,
          bounds,
          n_samp = cfg$mc_samp,
          hv_mc = 1000
        )
      )

      scores[!is.finite(scores)] <- -Inf
      best <- which.max(scores)

      x_next <- Xcand[best, , drop = FALSE]
      y_next <- co2_obj_julia(
        as.numeric(x_next[1, ]),
        cfg
      )

      history$X <- dplyr::bind_rows(
        history$X,
        x_next
      )
      history$Y <- rbind(
        history$Y,
        y_next
      )
      colnames(history$Y) <- c("f1", "f2", "f3")

      Ynorm_now <- normalize_to_unit(
        to_max3_mat(history$Y),
        bounds$lo,
        bounds$hi
      )
      front_now <- nondominated(Ynorm_now)
      hv_trace[n0 + it] <- compute_hv(
        front_now,
        ref = ref_fixed_norm,
        n_mc = 15000
      )

      if (debug) {
        message(sprintf(
          "[seed %02d | %s] iter %3d/%3d HV=%.6f",
          seed,
          method,
          it,
          cfg$n_iter,
          hv_trace[n0 + it]
        ))
      } else if (
        it == 1L ||
        it %% 10L == 0L ||
        it == cfg$n_iter
      ) {
        message(sprintf(
          "[seed %02d | %s] completed %3d/%3d iterations",
          seed,
          method,
          it,
          cfg$n_iter
        ))
      }

      if (!is.null(checkpoint_file) &&
          (it %% checkpoint_every == 0L ||
           it == cfg$n_iter)) {
        saveRDS(
          list(
            run_signature = run_signature,
            seed = seed,
            method = method,
            completed_iteration = it,
            n0 = n0,
            history = history,
            bounds = bounds,
            ref = ref_fixed_norm,
            hv_trace = hv_trace
          ),
          checkpoint_file
        )
      }

      # Remove per-iteration temporary objects and collect memory.
      rm(
        X_feat_hist,
        sur,
        Xcand,
        Xcand_feat,
        Yhist_norm,
        front0_norm,
        scores,
        x_next,
        y_next,
        Ynorm_now,
        front_now
      )
      if (it %% 10L == 0L) {
        invisible(gc(verbose = FALSE))
        try(julia_command("GC.gc()"), silent = TRUE)
      }
    }
  }

  result <- list(
    seed = seed,
    method = method,
    history = history,
    hv_trace = hv_trace,
    bounds = bounds,
    ref = ref_fixed_norm,
    cfg = cfg
  )

  if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) {
    unlink(checkpoint_file)
  }

  invisible(gc(verbose = FALSE))
  try(julia_command("GC.gc()"), silent = TRUE)

  result
}

# ------------------------------------------------------------
# 4) Run all methods for every seed
#
# The same X0 and Y0 are reused by TS, SCAL-UCB, EHVI, and EPI for
# each seed.
# ------------------------------------------------------------
run_multiseed_experiment <- function(cfg, seeds, methods, output_dir, debug = TRUE) {
  out <- setNames(vector("list", length(seeds)), sprintf("seed_%03d", seeds))

  for (s in seeds) {
    seed_start_time <- Sys.time()
    seed_name <- sprintf("seed_%03d", s)
    message("\n============================================================")
    message("STARTING ", seed_name)
    message("============================================================")

    initial_file <- file.path(
      output_dir, "runs",
      sprintf("initial_history_seed_%03d.rds", s)
    )

    if (file.exists(initial_file)) {
      initial_history <- readRDS(initial_file)
      message("Loaded existing shared initial design: ", initial_file)
    } else {
      initial_history <- make_initial_history(cfg, s, debug = debug)
      saveRDS(initial_history, initial_file)
    }

    out[[seed_name]] <- setNames(vector("list", length(methods)), methods)

    for (m in methods) {
      run_file <- file.path(
        output_dir, "runs",
        sprintf("result_seed_%03d_%s.rds", s, m)
      )

      if (file.exists(run_file)) {
        message("Loading completed run: seed=", s, ", method=", m)
        res <- readRDS(run_file)
      } else {
        message("Running: seed=", s, ", method=", m)
        checkpoint_file <- file.path(
          output_dir,
          "runs",
          sprintf("checkpoint_seed_%03d_%s.rds", s, m)
        )

        res <- run_bo(
          cfg = cfg,
          method = m,
          seed = s,
          initial_history = initial_history,
          debug = debug,
          checkpoint_file = checkpoint_file,
          checkpoint_every = 5L
        )
        saveRDS(res, run_file)
      }

      out[[seed_name]][[m]] <- res

      # Do not retain large temporary simulator allocations between methods.
      invisible(gc(verbose = FALSE))
      try(julia_command("GC.gc()"), silent = TRUE)
    }
    message(
      "Completed ", seed_name,
      " in ", round(as.numeric(difftime(Sys.time(), seed_start_time, units = "mins")), 2),
      " minutes."
    )
  }

  saveRDS(out, file.path(output_dir, "multiseed_results_all.rds"))
  out
}

multi_results <- run_multiseed_experiment(
  cfg = cfg,
  seeds = SEEDS,
  methods = METHODS,
  output_dir = OUTPUT_DIR,
  debug = FALSE
)

# Compatibility aliases: these allow the user's previous single-run plotting
primary_seed_name <- sprintf("seed_%03d", SEEDS[1])
res_TS   <- multi_results[[primary_seed_name]][["TS"]]
res_SCAL <- multi_results[[primary_seed_name]][["SCAL"]]
res_EHVI <- multi_results[[primary_seed_name]][["EHVI"]]
res_EPI  <- multi_results[[primary_seed_name]][["EPI"]]
res_real <- list(TS = res_TS, SCAL = res_SCAL, EHVI = res_EHVI, EPI = res_EPI)
cfg_real <- cfg
results <- res_real

# ------------------------------------------------------------
# 5) Collect every evaluated design from all repeated runs
# ------------------------------------------------------------
extract_run_history <- function(res) {
  X <- as.data.frame(res$history$X)
  Y <- as.data.frame(res$history$Y)
  names(Y)[1:3] <- c("f1", "f2", "f3")

  data.frame(
    seed = as.integer(res$seed),
    method = as.character(res$method),
    evaluation = seq_len(nrow(Y)),
    phase = ifelse(
      seq_len(nrow(Y)) <= cfg$n_init,
      "initial",
      "sequential"
    ),
    scenario = as.integer(X$scenario),
    rate_mt = as.numeric(X$rate_mt),
    L_mid = as.integer(X$L_mid),
    L_top = as.integer(X$L_top),
    f1_retained_integral_m3yr = as.numeric(Y$f1),
    f2_unretained_final_m3 = as.numeric(Y$f2),
    f3_total_cost_eur = as.numeric(Y$f3),
    stringsAsFactors = FALSE
  )
}

all_histories <- do.call(
  rbind,
  lapply(multi_results, function(seed_res) {
    do.call(rbind, lapply(seed_res, extract_run_history))
  })
)
rownames(all_histories) <- NULL

write.csv(
  all_histories,
  file.path(OUTPUT_DIR, "tables", "all_evaluated_designs_all_seeds.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 6) Optional dense grid for a common reference Pareto front
# ------------------------------------------------------------
build_reference_grid <- function(cfg, rate_grid_n = 21L) {
  expand.grid(
    scenario = c(1L, 3L),
    rate_mt = seq(0.8, 1.0, length.out = rate_grid_n),
    L_mid = seq_len(as.integer(setup_result$n_layers)),
    L_top = seq_len(as.integer(setup_result$n_layers)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}

evaluate_reference_grid <- function(cfg, output_dir, rate_grid_n = 21L) {
  grid_file <- file.path(
    output_dir, "reference",
    sprintf("reference_grid_nrate_%d.rds", rate_grid_n)
  )

  if (file.exists(grid_file)) {
    message("Loading existing reference grid: ", grid_file)
    return(readRDS(grid_file))
  }

  grid <- build_reference_grid(cfg, rate_grid_n)
  grid$f1 <- NA_real_
  grid$f2 <- NA_real_
  grid$f3 <- NA_real_

  checkpoint_file <- file.path(
    output_dir, "reference",
    sprintf("reference_grid_checkpoint_nrate_%d.rds", rate_grid_n)
  )

  if (file.exists(checkpoint_file)) {
    old <- readRDS(checkpoint_file)
    if (nrow(old) == nrow(grid)) {
      grid <- old
      message("Resuming the reference grid from its checkpoint.")
    }
  }

  todo <- which(!is.finite(grid$f1) | !is.finite(grid$f2) | !is.finite(grid$f3))

  for (ii in seq_along(todo)) {
    i <- todo[ii]
    y <- co2_obj_julia(
      as.numeric(grid[i, c("scenario", "rate_mt", "L_mid", "L_top")]),
      cfg
    )
    grid[i, c("f1", "f2", "f3")] <- y

    if (ii %% 25L == 0L || ii == length(todo)) {
      message(sprintf(
        "Reference grid: completed %d/%d pending designs.",
        ii, length(todo)
      ))
      saveRDS(grid, checkpoint_file)
    }
  }

  saveRDS(grid, grid_file)
  write.csv(
    grid,
    file.path(
      output_dir, "reference",
      sprintf("reference_grid_nrate_%d.csv", rate_grid_n)
    ),
    row.names = FALSE
  )
  grid
}

reference_grid <- NULL
if (REFERENCE_MODE %in% c("grid", "grid_only")) {
  reference_grid <- evaluate_reference_grid(
    cfg,
    OUTPUT_DIR,
    REFERENCE_RATE_GRID_N
  )
}

# ------------------------------------------------------------
# 7) Deterministic exact hypervolume for THREE objectives
#
# The acquisition functions above retain the user's working Monte Carlo HV.
# For reporting and statistical tests, this exact 3-D implementation removes
# Monte Carlo noise and guarantees that cumulative HV cannot decrease.
# ------------------------------------------------------------
hv2d_origin_exact <- function(P) {
  P <- as.matrix(P)
  if (nrow(P) == 0L) return(0)
  P <- P[is.finite(P[, 1]) & is.finite(P[, 2]), , drop = FALSE]
  P <- P[P[, 1] > 0 & P[, 2] > 0, , drop = FALSE]
  if (nrow(P) == 0L) return(0)

  y_breaks <- sort(unique(c(0, P[, 1])))
  area <- 0

  for (i in 2:length(y_breaks)) {
    y_left <- y_breaks[i - 1L]
    y_right <- y_breaks[i]
    active <- P[, 1] >= y_right
    if (any(active)) {
      max_z <- max(P[active, 2])
      area <- area + (y_right - y_left) * max_z
    }
  }
  area
}

compute_hv_exact_3d <- function(Y, ref = c(-0.10, -0.10, -0.10)) {
  P <- as.matrix(Y)
  if (nrow(P) == 0L) return(0)
  if (ncol(P) != 3L) stop("compute_hv_exact_3d requires exactly 3 objectives.")

  P <- P[apply(P, 1, function(z) all(is.finite(z))), , drop = FALSE]
  if (nrow(P) == 0L) return(0)

  P <- nondominated(P)
  P <- unique(P)
  Q <- sweep(P, 2, as.numeric(ref), "-")
  Q <- Q[apply(Q, 1, function(z) all(z > 0)), , drop = FALSE]
  if (nrow(Q) == 0L) return(0)

  x_breaks <- sort(unique(c(0, Q[, 1])))
  volume <- 0

  for (i in 2:length(x_breaks)) {
    x_left <- x_breaks[i - 1L]
    x_right <- x_breaks[i]
    active <- Q[, 1] >= x_right
    if (any(active)) {
      area_yz <- hv2d_origin_exact(Q[active, 2:3, drop = FALSE])
      volume <- volume + (x_right - x_left) * area_yz
    }
  }

  as.numeric(volume)
}

# Small internal checks for the exact HV routine.
stopifnot(abs(compute_hv_exact_3d(matrix(c(1, 1, 1), nrow = 1), c(0, 0, 0)) - 1) < 1e-10)
.test_hv <- compute_hv_exact_3d(
  rbind(c(1, 0.5, 1), c(0.5, 1, 1)),
  c(0, 0, 0)
)
stopifnot(abs(.test_hv - 0.75) < 1e-10)
rm(.test_hv)

# ------------------------------------------------------------
# 8) Common normalization, grid-only reference front, and union reference
# ------------------------------------------------------------
run_Y_raw <- as.matrix(
  all_histories[, c(
    "f1_retained_integral_m3yr",
    "f2_unretained_final_m3",
    "f3_total_cost_eur"
  )]
)

if (is.null(reference_grid)) {
  stop("REFERENCE_MODE requires the completed dense reference grid.")
}

grid_Y_raw <- as.matrix(reference_grid[, c("f1", "f2", "f3")])

# A common normalization range is taken from the union so that an BO point
# that falls between grid rate levels is not silently clipped merely because
# the grid is discrete. The PRIMARY IGD+ reference front itself is grid-only.
normalization_pool_raw <- rbind(grid_Y_raw, run_Y_raw)
normalization_pool_max <- to_max3_mat(normalization_pool_raw)

global_lo <- apply(normalization_pool_max, 2, min, na.rm = TRUE)
global_hi <- apply(normalization_pool_max, 2, max, na.rm = TRUE)
global_span <- global_hi - global_lo
global_span[!is.finite(global_span) | global_span <= 0] <- 1
global_hi <- global_lo + global_span

global_bounds <- list(lo = global_lo, hi = global_hi)
evaluation_ref_point <- rep(-0.10, 3)

# Primary independent grid-only reference front.
grid_norm <- normalize_to_unit(
  to_max3_mat(grid_Y_raw),
  global_bounds$lo,
  global_bounds$hi
)
reference_front_grid_norm <- unique(nondominated(grid_norm))
reference_front_norm <- reference_front_grid_norm

# Supplementary union front for a robustness check.
union_norm <- normalize_to_unit(
  to_max3_mat(normalization_pool_raw),
  global_bounds$lo,
  global_bounds$hi
)
reference_front_union_norm <- unique(nondominated(union_norm))

colnames(reference_front_grid_norm) <- c(
  "normalized_storage",
  "normalized_low_unretained_volume",
  "normalized_low_cost"
)
colnames(reference_front_union_norm) <- colnames(reference_front_grid_norm)

write.csv(
  reference_front_grid_norm,
  file.path(OUTPUT_DIR, "reference", "grid_only_reference_pareto_front_normalized.csv"),
  row.names = FALSE
)
write.csv(
  reference_front_union_norm,
  file.path(OUTPUT_DIR, "reference", "grid_plus_BO_union_reference_front_normalized.csv"),
  row.names = FALSE
)

reference_definition <- data.frame(
  item = c(
    "primary_reference_mode",
    "primary_grid_reference_points",
    "supplementary_union_reference_points",
    "dense_grid_designs",
    "normalization_lower_bound_max_space",
    "normalization_upper_bound_max_space",
    "hypervolume_reference_point_normalized",
    "primary_IGD_plus_reference",
    "supplementary_IGD_plus_reference"
  ),
  value = c(
    "grid_only",
    nrow(reference_front_grid_norm),
    nrow(reference_front_union_norm),
    nrow(reference_grid),
    paste(signif(global_bounds$lo, 8), collapse = ","),
    paste(signif(global_bounds$hi, 8), collapse = ","),
    paste(evaluation_ref_point, collapse = ","),
    "nondominated front from the 3402-design dense grid only",
    "nondominated union of the grid and all repeated BO evaluations"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  reference_definition,
  file.path(OUTPUT_DIR, "tables", "reference_front_definition.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9) IGD+ for maximization objectives
#
# For a reference point r and approximation point a, the positive
# directional distance in maximization space is max(r - a, 0).
# ------------------------------------------------------------
compute_igd_plus_max <- function(approx_front, reference_front) {
  A <- as.matrix(approx_front)
  R <- as.matrix(reference_front)

  if (nrow(A) == 0L || nrow(R) == 0L) return(NA_real_)

  distances <- vapply(seq_len(nrow(R)), function(i) {
    r <- R[i, ]
    diff <- matrix(r, nrow = nrow(A), ncol = ncol(A), byrow = TRUE) - A
    diff <- pmax(diff, 0)
    min(sqrt(rowSums(diff^2)))
  }, numeric(1))

  mean(distances)
}

# ------------------------------------------------------------
# 10) Pareto diversity measures
#
# spacing_sd is the SD of nearest-neighbour distances.
# objective_extent is the sum of the ranges covered in the 3 normalized
# objectives; it helps distinguish an evenly spaced but very narrow cluster
# from broad Pareto-front coverage.
# ------------------------------------------------------------
pareto_diversity_metrics <- function(front) {
  F <- unique(as.matrix(front))
  if (nrow(F) == 0L) {
    return(c(
      pareto_size = 0,
      mean_nearest_neighbor = NA_real_,
      spacing_sd = NA_real_,
      objective_extent = 0
    ))
  }
  if (nrow(F) == 1L) {
    return(c(
      pareto_size = 1,
      mean_nearest_neighbor = NA_real_,
      spacing_sd = NA_real_,
      objective_extent = 0
    ))
  }

  D <- as.matrix(stats::dist(F))
  diag(D) <- Inf
  nearest <- apply(D, 1, min)

  c(
    pareto_size = nrow(F),
    mean_nearest_neighbor = mean(nearest),
    spacing_sd = stats::sd(nearest),
    objective_extent = sum(apply(F, 2, function(z) diff(range(z))))
  )
}

# ------------------------------------------------------------
# 11) Recompute HV and IGD+ at every evaluation for every run
#     using the same bounds and reference front.
# ------------------------------------------------------------
calculate_run_metrics <- function(res, global_bounds, reference_front, ref_point) {
  Yraw <- as.matrix(res$history$Y)
  Ymax <- to_max3_mat(Yraw)
  Ynorm <- normalize_to_unit(Ymax, global_bounds$lo, global_bounds$hi)

  n <- nrow(Ynorm)
  out <- data.frame(
    seed = rep(as.integer(res$seed), n),
    method = rep(as.character(res$method), n),
    evaluation = seq_len(n),
    hypervolume = NA_real_,
    IGD_plus = NA_real_,
    pareto_size = NA_integer_,
    stringsAsFactors = FALSE
  )

  for (k in seq_len(n)) {
    front_k <- nondominated(Ynorm[seq_len(k), , drop = FALSE])
    front_k <- unique(front_k)

    out$hypervolume[k] <- compute_hv_exact_3d(front_k, ref = ref_point)
    out$IGD_plus[k] <- compute_igd_plus_max(front_k, reference_front)
    out$pareto_size[k] <- nrow(front_k)
  }

  final_front <- unique(nondominated(Ynorm))
  diversity <- pareto_diversity_metrics(final_front)

  final <- data.frame(
    seed = as.integer(res$seed),
    method = as.character(res$method),
    final_hypervolume = tail(out$hypervolume, 1),
    final_IGD_plus = tail(out$IGD_plus, 1),
    final_pareto_size = unname(diversity["pareto_size"]),
    mean_nearest_neighbor = unname(diversity["mean_nearest_neighbor"]),
    spacing_sd = unname(diversity["spacing_sd"]),
    objective_extent = unname(diversity["objective_extent"]),
    best_f1_retained_integral_m3yr = max(Yraw[, 1], na.rm = TRUE),
    best_f2_unretained_final_m3 = min(Yraw[, 2], na.rm = TRUE),
    best_f3_total_cost_eur = min(Yraw[, 3], na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  list(trace = out, final = final)
}

metric_objects <- lapply(multi_results, function(seed_res) {
  lapply(seed_res, calculate_run_metrics,
         global_bounds = global_bounds,
         reference_front = reference_front_norm,
         ref_point = evaluation_ref_point)
})

metric_trace <- do.call(
  rbind,
  lapply(metric_objects, function(seed_obj) {
    do.call(rbind, lapply(seed_obj, `[[`, "trace"))
  })
)
rownames(metric_trace) <- NULL

final_metrics <- do.call(
  rbind,
  lapply(metric_objects, function(seed_obj) {
    do.call(rbind, lapply(seed_obj, `[[`, "final"))
  })
)
rownames(final_metrics) <- NULL

# Supplementary robustness check: final IGD+ against the grid+BO union front.
calculate_final_union_igd <- function(res) {
  Yraw <- as.matrix(res$history$Y)
  Ynorm <- normalize_to_unit(
    to_max3_mat(Yraw),
    global_bounds$lo,
    global_bounds$hi
  )
  final_front <- unique(nondominated(Ynorm))
  data.frame(
    seed = as.integer(res$seed),
    method = as.character(res$method),
    final_IGD_plus_grid_only = compute_igd_plus_max(
      final_front,
      reference_front_grid_norm
    ),
    final_IGD_plus_grid_plus_BO_union = compute_igd_plus_max(
      final_front,
      reference_front_union_norm
    ),
    stringsAsFactors = FALSE
  )
}

igd_reference_robustness <- do.call(
  rbind,
  lapply(multi_results, function(seed_res) {
    do.call(rbind, lapply(seed_res, calculate_final_union_igd))
  })
)
igd_reference_robustness$difference_union_minus_grid <- (
  igd_reference_robustness$final_IGD_plus_grid_plus_BO_union -
  igd_reference_robustness$final_IGD_plus_grid_only
)

igd_reference_robustness_summary <- do.call(
  rbind,
  lapply(split(igd_reference_robustness, igd_reference_robustness$method), function(d) {
    data.frame(
      method = d$method[1],
      grid_only_mean = mean(d$final_IGD_plus_grid_only),
      grid_only_sd = sd(d$final_IGD_plus_grid_only),
      union_mean = mean(d$final_IGD_plus_grid_plus_BO_union),
      union_sd = sd(d$final_IGD_plus_grid_plus_BO_union),
      mean_difference_union_minus_grid = mean(d$difference_union_minus_grid),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  igd_reference_robustness,
  file.path(OUTPUT_DIR, "tables", "IGDplus_grid_only_vs_union_by_seed_method.csv"),
  row.names = FALSE
)
write.csv(
  igd_reference_robustness_summary,
  file.path(OUTPUT_DIR, "tables", "IGDplus_grid_only_vs_union_mean_sd.csv"),
  row.names = FALSE
)

write.csv(
  metric_trace,
  file.path(OUTPUT_DIR, "tables", "HV_IGDplus_convergence_all_seeds.csv"),
  row.names = FALSE
)
write.csv(
  final_metrics,
  file.path(OUTPUT_DIR, "tables", "final_metrics_by_seed_and_method.csv"),
  row.names = FALSE
)

# Confirm that exact cumulative HV is non-decreasing, apart from numerical
# roundoff. This directly addresses the earlier negative-increment issue.
hv_monotonicity_check <- aggregate(
  hypervolume ~ seed + method,
  metric_trace,
  function(x) min(diff(x))
)
names(hv_monotonicity_check)[3] <- "minimum_increment"
hv_monotonicity_check$passes <- (
  hv_monotonicity_check$minimum_increment >= -1e-12
)

write.csv(
  hv_monotonicity_check,
  file.path(OUTPUT_DIR, "tables", "hypervolume_monotonicity_check.csv"),
  row.names = FALSE
)

if (!all(hv_monotonicity_check$passes)) {
  warning("At least one exact-HV trace decreased; inspect hypervolume_monotonicity_check.csv")
}

# ------------------------------------------------------------
# 12) Mean, SD, SE, and 95% confidence intervals
# ------------------------------------------------------------
summarize_mean_sd_ci <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0L) {
    return(c(n = 0, mean = NA, sd = NA, se = NA, ci_low = NA, ci_high = NA))
  }
  m <- mean(x)
  s <- if (n > 1L) stats::sd(x) else NA_real_
  se <- if (n > 1L) s / sqrt(n) else NA_real_
  crit <- if (n > 1L) stats::qt(0.975, df = n - 1L) else NA_real_
  c(
    n = n,
    mean = m,
    sd = s,
    se = se,
    ci_low = if (n > 1L) m - crit * se else NA_real_,
    ci_high = if (n > 1L) m + crit * se else NA_real_
  )
}

summary_metric_names <- c(
  "final_hypervolume",
  "final_IGD_plus",
  "final_pareto_size",
  "mean_nearest_neighbor",
  "spacing_sd",
  "objective_extent",
  "best_f1_retained_integral_m3yr",
  "best_f2_unretained_final_m3",
  "best_f3_total_cost_eur"
)

final_summary <- do.call(
  rbind,
  lapply(METHODS, function(m) {
    dm <- final_metrics[final_metrics$method == m, , drop = FALSE]
    do.call(
      rbind,
      lapply(summary_metric_names, function(metric) {
        sm <- summarize_mean_sd_ci(dm[[metric]])
        data.frame(
          method = m,
          metric = metric,
          n = unname(sm["n"]),
          mean = unname(sm["mean"]),
          sd = unname(sm["sd"]),
          se = unname(sm["se"]),
          ci_low_95 = unname(sm["ci_low"]),
          ci_high_95 = unname(sm["ci_high"]),
          stringsAsFactors = FALSE
        )
      })
    )
  })
)

write.csv(
  final_summary,
  file.path(OUTPUT_DIR, "tables", "final_performance_mean_sd_95CI.csv"),
  row.names = FALSE
)

# Compact manuscript table: mean ± SD for the principal metrics.
manuscript_summary <- reshape(
  final_summary[
    final_summary$metric %in% c(
      "final_hypervolume",
      "final_IGD_plus",
      "final_pareto_size",
      "spacing_sd",
      "objective_extent"
    ),
    c("method", "metric", "mean", "sd", "ci_low_95", "ci_high_95")
  ],
  idvar = "method",
  timevar = "metric",
  direction = "wide"
)

write.csv(
  manuscript_summary,
  file.path(OUTPUT_DIR, "tables", "manuscript_multiseed_summary_table.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 13) Paired statistical tests
#
# Seeds are blocks. Friedman tests compare all four methods, followed by
# paired Wilcoxon signed-rank tests with Holm correction.
# ------------------------------------------------------------
make_seed_method_matrix <- function(data, metric, methods = METHODS) {
  wide <- reshape(
    data[, c("seed", "method", metric)],
    idvar = "seed",
    timevar = "method",
    direction = "wide"
  )
  wanted <- paste0(metric, ".", methods)
  missing <- setdiff(wanted, names(wide))
  if (length(missing) > 0L) {
    stop("Missing seed-method columns for metric ", metric, ": ",
         paste(missing, collapse = ", "))
  }
  M <- as.matrix(wide[, wanted, drop = FALSE])
  colnames(M) <- methods
  rownames(M) <- wide$seed
  M <- M[stats::complete.cases(M), , drop = FALSE]
  M
}

run_friedman_test <- function(data, metric) {
  M <- make_seed_method_matrix(data, metric)

  if (nrow(M) < 2L) {
    return(data.frame(
      metric = metric,
      n_complete_seeds = nrow(M),
      statistic = NA_real_,
      df = NA_real_,
      p_value = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  test <- stats::friedman.test(M)
  data.frame(
    metric = metric,
    n_complete_seeds = nrow(M),
    statistic = unname(test$statistic),
    df = unname(test$parameter),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}

run_paired_wilcoxon <- function(data, metric, methods = METHODS) {
  pairs <- combn(methods, 2, simplify = FALSE)

  ans <- do.call(rbind, lapply(pairs, function(pair) {
    d1 <- data[data$method == pair[1], c("seed", metric)]
    d2 <- data[data$method == pair[2], c("seed", metric)]
    names(d1)[2] <- "x"
    names(d2)[2] <- "y"
    merged <- merge(d1, d2, by = "seed")
    merged <- merged[is.finite(merged$x) & is.finite(merged$y), ]

    wt <- tryCatch(
      stats::wilcox.test(
        merged$x,
        merged$y,
        paired = TRUE,
        exact = FALSE,
        conf.int = FALSE
      ),
      error = function(e) NULL
    )

    data.frame(
      metric = metric,
      method_1 = pair[1],
      method_2 = pair[2],
      n_pairs = nrow(merged),
      mean_difference_method1_minus_method2 = mean(merged$x - merged$y),
      median_difference_method1_minus_method2 = stats::median(merged$x - merged$y),
      W = if (is.null(wt)) NA_real_ else unname(wt$statistic),
      p_value_raw = if (is.null(wt)) NA_real_ else wt$p.value,
      stringsAsFactors = FALSE
    )
  }))

  ans$p_value_Holm <- stats::p.adjust(ans$p_value_raw, method = "holm")
  ans
}

test_metrics <- c(
  "final_hypervolume",
  "final_IGD_plus",
  "final_pareto_size",
  "spacing_sd",
  "objective_extent"
)

friedman_results <- do.call(
  rbind,
  lapply(test_metrics, function(metric) {
    run_friedman_test(final_metrics, metric)
  })
)

wilcoxon_results <- do.call(
  rbind,
  lapply(test_metrics, function(metric) {
    run_paired_wilcoxon(final_metrics, metric)
  })
)

write.csv(
  friedman_results,
  file.path(OUTPUT_DIR, "tables", "friedman_tests_all_methods.csv"),
  row.names = FALSE
)
write.csv(
  wilcoxon_results,
  file.path(OUTPUT_DIR, "tables", "paired_wilcoxon_Holm_tests.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 14) Convergence summaries across seeds
# ------------------------------------------------------------
summarize_trace_metric <- function(data, metric) {
  groups <- split(data, interaction(data$method, data$evaluation, drop = TRUE))

  out <- do.call(rbind, lapply(groups, function(d) {
    sm <- summarize_mean_sd_ci(d[[metric]])
    data.frame(
      method = d$method[1],
      evaluation = d$evaluation[1],
      metric = metric,
      n = unname(sm["n"]),
      mean = unname(sm["mean"]),
      sd = unname(sm["sd"]),
      se = unname(sm["se"]),
      ci_low_95 = unname(sm["ci_low"]),
      ci_high_95 = unname(sm["ci_high"]),
      stringsAsFactors = FALSE
    )
  }))
  out[order(out$method, out$evaluation), ]
}

hv_convergence_summary <- summarize_trace_metric(metric_trace, "hypervolume")
igd_convergence_summary <- summarize_trace_metric(metric_trace, "IGD_plus")

write.csv(
  hv_convergence_summary,
  file.path(OUTPUT_DIR, "tables", "HV_convergence_mean_sd_95CI.csv"),
  row.names = FALSE
)
write.csv(
  igd_convergence_summary,
  file.path(OUTPUT_DIR, "tables", "IGDplus_convergence_mean_sd_95CI.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 15) Consistent tables for all designs and Pareto designs
# ------------------------------------------------------------
all_design_summary <- do.call(
  rbind,
  lapply(METHODS, function(m) {
    d <- all_histories[all_histories$method == m, ]
    data.frame(
      method = m,
      population = "all evaluated designs across all repeated seeds",
      n_evaluations = nrow(d),
      max_f1_retained_integral_m3yr =
        max(d$f1_retained_integral_m3yr, na.rm = TRUE),
      min_f2_unretained_final_m3 =
        min(d$f2_unretained_final_m3, na.rm = TRUE),
      min_f3_total_cost_eur =
        min(d$f3_total_cost_eur, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

pooled_pareto_by_method <- do.call(
  rbind,
  lapply(METHODS, function(m) {
    d <- all_histories[all_histories$method == m, ]
    Yraw <- as.matrix(d[, c(
      "f1_retained_integral_m3yr",
      "f2_unretained_final_m3",
      "f3_total_cost_eur"
    )])
    Ynorm <- normalize_to_unit(
      to_max3_mat(Yraw),
      global_bounds$lo,
      global_bounds$hi
    )
    keep <- apply(Yraw, 1, function(z) all(is.finite(z)))
    idx_valid <- which(keep)
    nd_local <- nondominated(Ynorm[keep, , drop = FALSE])

    # Match unique nondominated normalized rows back to the valid data.
    keys_all <- apply(
      round(Ynorm[keep, , drop = FALSE], 12),
      1,
      paste,
      collapse = "|"
    )
    keys_nd <- unique(apply(round(nd_local, 12), 1, paste, collapse = "|"))
    selected <- idx_valid[keys_all %in% keys_nd]
    p <- d[selected, , drop = FALSE]

    data.frame(
      method = m,
      population = "pooled method-specific Pareto designs across all repeated seeds",
      n_pareto_rows = nrow(p),
      max_f1_retained_integral_m3yr =
        max(p$f1_retained_integral_m3yr, na.rm = TRUE),
      min_f2_unretained_final_m3 =
        min(p$f2_unretained_final_m3, na.rm = TRUE),
      min_f3_total_cost_eur =
        min(p$f3_total_cost_eur, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  all_design_summary,
  file.path(OUTPUT_DIR, "tables", "summary_all_evaluated_designs.csv"),
  row.names = FALSE
)
write.csv(
  pooled_pareto_by_method,
  file.path(OUTPUT_DIR, "tables", "summary_pooled_Pareto_designs.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 17) statistical plots
# ------------------------------------------------------------
method_cols <- c(
  TS = "#00A6D6",
  SCAL = "#F28E2B",
  EHVI = "#8E44AD",
  EPI = "#E31A1C"
)

theme_revision <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, color = "grey50"),
      axis.text = ggplot2::element_text(color = "black")
    )
}

plot_trace_with_ci <- function(summary_data, y_label, file_stem) {
  p <- ggplot2::ggplot(
    summary_data[summary_data$evaluation >= cfg$n_init, ],
    ggplot2::aes(
      x = evaluation,
      y = mean,
      color = method,
      fill = method
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_low_95, ymax = ci_high_95),
      alpha = 0.16,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_manual(values = method_cols) +
    ggplot2::scale_fill_manual(values = method_cols) +
    ggplot2::labs(
      x = "Evaluation",
      y = y_label
    ) +
    theme_revision()

  ggplot2::ggsave(
    file.path(OUTPUT_DIR, "figures", paste0(file_stem, ".pdf")),
    p, width = 7.2, height = 4.8
  )
  ggplot2::ggsave(
    file.path(OUTPUT_DIR, "figures", paste0(file_stem, ".png")),
    p, width = 7.2, height = 4.8, dpi = 350
  )
  print(p)
  invisible(p)
}

p_hv_multiseed <- plot_trace_with_ci(
  hv_convergence_summary,
  "Hypervolume (mean and 95% CI)",
  "multiseed_HV_convergence_95CI"
)

p_igd_multiseed <- plot_trace_with_ci(
  igd_convergence_summary,
  "IGD+ (mean and 95% CI; lower is better)",
  "multiseed_IGDplus_convergence_95CI"
)

p_hv_box <- ggplot2::ggplot(
  final_metrics,
  ggplot2::aes(x = method, y = final_hypervolume, fill = method)
) +
  ggplot2::geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.65) +
  ggplot2::geom_jitter(width = 0.08, size = 2) +
  ggplot2::scale_fill_manual(values = method_cols) +
  ggplot2::labs(x = "Method", y = "Final hypervolume") +
  theme_revision() +
  ggplot2::theme(legend.position = "none")

p_igd_box <- ggplot2::ggplot(
  final_metrics,
  ggplot2::aes(x = method, y = final_IGD_plus, fill = method)
) +
  ggplot2::geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.65) +
  ggplot2::geom_jitter(width = 0.08, size = 2) +
  ggplot2::scale_fill_manual(values = method_cols) +
  ggplot2::labs(x = "Method", y = "Final IGD+ (lower is better)") +
  theme_revision() +
  ggplot2::theme(legend.position = "none")

ggplot2::ggsave(
  file.path(OUTPUT_DIR, "figures", "final_HV_boxplot_all_seeds.pdf"),
  p_hv_box, width = 6.5, height = 4.8
)
ggplot2::ggsave(
  file.path(OUTPUT_DIR, "figures", "final_HV_boxplot_all_seeds.png"),
  p_hv_box, width = 6.5, height = 4.8, dpi = 350
)
ggplot2::ggsave(
  file.path(OUTPUT_DIR, "figures", "final_IGDplus_boxplot_all_seeds.pdf"),
  p_igd_box, width = 6.5, height = 4.8
)
ggplot2::ggsave(
  file.path(OUTPUT_DIR, "figures", "final_IGDplus_boxplot_all_seeds.png"),
  p_igd_box, width = 6.5, height = 4.8, dpi = 350
)

print(p_hv_box)
print(p_igd_box)

# ------------------------------------------------------------
# 18) SCAL-UCB kappa and EHVI MC sensitivity
# ------------------------------------------------------------
evaluate_one_result_final <- function(res) {
  obj <- calculate_run_metrics(
    res,
    global_bounds,
    reference_front_norm,
    evaluation_ref_point
  )
  obj$final
}

run_sensitivity_family <- function(
    family_name,
    method,
    parameter_name,
    parameter_values,
    seeds,
    cfg_base
) {
  rows <- list()
  counter <- 0L

  for (value in parameter_values) {
    for (s in seeds) {
      cfg_s <- cfg_base
      if (identical(parameter_name, "ucb_k")) cfg_s$ucb_k <- as.numeric(value)
      if (identical(parameter_name, "mc_samp")) cfg_s$mc_samp <- as.integer(value)

      initial_history <- make_initial_history(cfg_s, s, debug = FALSE)
      run_file <- file.path(
        OUTPUT_DIR, "sensitivity",
        sprintf(
          "%s_%s_%s_seed_%03d.rds",
          family_name,
          parameter_name,
          gsub("\\.", "p", as.character(value)),
          s
        )
      )

      if (file.exists(run_file)) {
        res <- readRDS(run_file)
      } else {
        res <- run_bo(
          cfg_s,
          method = method,
          seed = s,
          initial_history = initial_history,
          debug = FALSE
        )
        saveRDS(res, run_file)
      }

      fm <- evaluate_one_result_final(res)
      counter <- counter + 1L
      fm$sensitivity_family <- family_name
      fm$parameter <- parameter_name
      fm$parameter_value <- as.numeric(value)
      rows[[counter]] <- fm
      message(
        "Sensitivity complete: ", family_name,
        ", ", parameter_name, "=", value,
        ", seed=", s
      )
    }
  }

  do.call(rbind, rows)
}

if (isTRUE(RUN_HYPERPARAMETER_SENSITIVITY)) {
  kappa_sensitivity <- run_sensitivity_family(
    family_name = "SCAL_UCB_kappa",
    method = "SCAL",
    parameter_name = "ucb_k",
    parameter_values = KAPPA_VALUES,
    seeds = SENSITIVITY_SEEDS,
    cfg_base = cfg
  )

  ehvi_mc_sensitivity <- run_sensitivity_family(
    family_name = "EHVI_MC_samples",
    method = "EHVI",
    parameter_name = "mc_samp",
    parameter_values = EHVI_MC_VALUES,
    seeds = SENSITIVITY_SEEDS,
    cfg_base = cfg
  )

  hyperparameter_sensitivity <- rbind(
    kappa_sensitivity,
    ehvi_mc_sensitivity
  )

  write.csv(
    hyperparameter_sensitivity,
    file.path(
      OUTPUT_DIR,
      "tables",
      "hyperparameter_sensitivity_results_by_seed.csv"
    ),
    row.names = FALSE
  )

  sensitivity_summary <- do.call(
    rbind,
    lapply(split(
      hyperparameter_sensitivity,
      interaction(
        hyperparameter_sensitivity$sensitivity_family,
        hyperparameter_sensitivity$parameter_value,
        drop = TRUE
      )
    ), function(d) {
      hv <- summarize_mean_sd_ci(d$final_hypervolume)
      igd <- summarize_mean_sd_ci(d$final_IGD_plus)

      data.frame(
        sensitivity_family = d$sensitivity_family[1],
        parameter = d$parameter[1],
        parameter_value = d$parameter_value[1],
        n_seeds = nrow(d),
        HV_mean = unname(hv["mean"]),
        HV_sd = unname(hv["sd"]),
        HV_ci_low_95 = unname(hv["ci_low"]),
        HV_ci_high_95 = unname(hv["ci_high"]),
        IGDplus_mean = unname(igd["mean"]),
        IGDplus_sd = unname(igd["sd"]),
        IGDplus_ci_low_95 = unname(igd["ci_low"]),
        IGDplus_ci_high_95 = unname(igd["ci_high"]),
        stringsAsFactors = FALSE
      )
    })
  )

  write.csv(
    sensitivity_summary,
    file.path(
      OUTPUT_DIR,
      "tables",
      "hyperparameter_sensitivity_mean_sd_95CI.csv"
    ),
    row.names = FALSE
  )
}


# ============================================================
# 19) HELD-OUT SURROGATE VALIDATION
#     Fair common-training comparison: MNL vs mixed-kernel GP vs RF
# ============================================================

RUN_SURROGATE_VALIDATION <- TRUE
SURROGATE_SPLIT_SEEDS <- 2001:2005
SURROGATE_TEST_FRACTION <- 0.20
SURROGATE_COMMON_TRAIN_N <- 180L
SURROGATE_RF_TREES <- 500L

if (isTRUE(RUN_SURROGATE_VALIDATION)) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop(
      "The held-out surrogate comparison requires package 'ranger'. ",
      "Install it with install.packages('ranger') and rerun."
    )
  }

  make_unique_surrogate_dataset <- function(history_df) {
    d <- history_df
    d$rate_key <- round(d$rate_mt, 6)
    aggregate(
      cbind(
        f1_retained_integral_m3yr,
        f2_unretained_final_m3,
        f3_total_cost_eur
      ) ~ scenario + rate_key + L_mid + L_top,
      data = d,
      FUN = mean
    )
  }

  surrogate_data <- make_unique_surrogate_dataset(all_histories)
  names(surrogate_data)[names(surrogate_data) == "rate_key"] <- "rate_mt"

  prepare_surrogate_predictors <- function(d) {
    d <- as.data.frame(d)
    d$scenario <- factor(d$scenario, levels = c(1, 3))
    d$L_mid <- factor(
      d$L_mid,
      levels = as.character(seq_len(setup_result$n_layers))
    )
    d$L_top <- factor(
      d$L_top,
      levels = as.character(seq_len(setup_result$n_layers))
    )
    d
  }

  regression_metrics <- function(observed, predicted, train_sd) {
    observed <- as.numeric(observed)
    predicted <- as.numeric(predicted)
    keep <- is.finite(observed) & is.finite(predicted)
    observed <- observed[keep]
    predicted <- predicted[keep]

    if (length(observed) < 2L) {
      return(c(RMSE = NA, NRMSE_trainSD = NA, MAE = NA, R2 = NA))
    }

    rmse <- sqrt(mean((observed - predicted)^2))
    mae <- mean(abs(observed - predicted))
    denom <- sum((observed - mean(observed))^2)
    r2 <- if (denom > 0) {
      1 - sum((observed - predicted)^2) / denom
    } else {
      NA_real_
    }

    c(
      RMSE = rmse,
      NRMSE_trainSD = if (is.finite(train_sd) && train_sd > 0) rmse / train_sd else NA_real_,
      MAE = mae,
      R2 = r2
    )
  }

  fit_predict_mnl_validation <- function(train, test, B = 5L) {
    objective_names <- c(
      "f1_retained_integral_m3yr",
      "f2_unretained_final_m3",
      "f3_total_cost_eur"
    )
    Y_train <- as.matrix(train[, objective_names, drop = FALSE])
    breaks_list <- vector("list", length(objective_names))
    bin_matrix <- matrix(NA_integer_, nrow(Y_train), length(objective_names))

    for (j in seq_along(objective_names)) {
      breaks_list[[j]] <- .make_quantile_breaks(Y_train[, j], B)
      bin_matrix[, j] <- findInterval(
        Y_train[, j], breaks_list[[j]],
        rightmost.closed = TRUE, all.inside = TRUE
      )
    }

    class_label <- apply(bin_matrix, 1, paste, collapse = "_")
    class_factor <- factor(class_label)
    train_model <- prepare_surrogate_predictors(train)
    test_model <- prepare_surrogate_predictors(test)
    centroid_df <- aggregate(
      Y_train,
      by = list(class = class_label),
      FUN = mean
    )
    class_levels <- levels(class_factor)

    if (length(class_levels) < 2L) {
      pred <- matrix(colMeans(Y_train), nrow(test), ncol(Y_train), byrow = TRUE)
      colnames(pred) <- objective_names
      return(list(prediction = pred, accuracy = NA_real_, log_loss = NA_real_))
    }

    fit_df <- data.frame(
      class = class_factor,
      train_model[, c("scenario", "rate_mt", "L_mid", "L_top")]
    )
    fit <- tryCatch(
      nnet::multinom(
        class ~ scenario + rate_mt + L_mid + L_top,
        data = fit_df,
        trace = FALSE,
        decay = 1e-4,
        maxit = 1000,
        MaxNWts = 100000
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      pred <- matrix(colMeans(Y_train), nrow(test), ncol(Y_train), byrow = TRUE)
      colnames(pred) <- objective_names
      return(list(prediction = pred, accuracy = NA_real_, log_loss = NA_real_))
    }

    raw_prob <- tryCatch(
      predict(fit, newdata = test_model, type = "probs"),
      error = function(e) NULL
    )
    if (is.null(raw_prob)) {
      pred <- matrix(colMeans(Y_train), nrow(test), ncol(Y_train), byrow = TRUE)
      colnames(pred) <- objective_names
      return(list(prediction = pred, accuracy = NA_real_, log_loss = NA_real_))
    }

    if (is.null(dim(raw_prob))) {
      if (length(class_levels) == 2L) {
        p2 <- as.numeric(raw_prob)
        raw_prob <- cbind(1 - p2, p2)
        colnames(raw_prob) <- class_levels
      } else {
        raw_prob <- matrix(raw_prob, nrow = 1)
      }
    }

    aligned_prob <- matrix(
      0, nrow(test), length(class_levels),
      dimnames = list(NULL, class_levels)
    )
    common <- intersect(colnames(raw_prob), class_levels)
    aligned_prob[, common] <- raw_prob[, common, drop = FALSE]
    aligned_prob[!is.finite(aligned_prob)] <- 0
    zero_rows <- rowSums(aligned_prob) <= 0
    if (any(zero_rows)) aligned_prob[zero_rows, ] <- 1 / ncol(aligned_prob)
    aligned_prob <- aligned_prob / rowSums(aligned_prob)

    centroid_matrix <- matrix(
      NA_real_, length(class_levels), length(objective_names),
      dimnames = list(class_levels, objective_names)
    )
    global_centroid <- colMeans(Y_train)
    for (cl in class_levels) {
      row <- centroid_df[centroid_df$class == cl, , drop = FALSE]
      centroid_matrix[cl, ] <- if (nrow(row) == 1L) {
        as.numeric(row[1, objective_names, drop = TRUE])
      } else {
        global_centroid
      }
    }
    pred <- aligned_prob %*% centroid_matrix

    test_Y <- as.matrix(test[, objective_names, drop = FALSE])
    test_bins <- matrix(NA_integer_, nrow(test_Y), length(objective_names))
    for (j in seq_along(objective_names)) {
      test_bins[, j] <- findInterval(
        test_Y[, j], breaks_list[[j]],
        rightmost.closed = TRUE, all.inside = TRUE
      )
    }
    true_class <- apply(test_bins, 1, paste, collapse = "_")
    predicted_class <- colnames(aligned_prob)[max.col(aligned_prob)]
    accuracy <- mean(predicted_class == true_class)
    true_probability <- vapply(seq_len(nrow(test)), function(i) {
      if (true_class[i] %in% colnames(aligned_prob)) {
        aligned_prob[i, true_class[i]]
      } else {
        1e-15
      }
    }, numeric(1))

    list(
      prediction = pred,
      accuracy = accuracy,
      log_loss = -mean(log(pmax(true_probability, 1e-15)))
    )
  }

  mixed_kernel_matrix <- function(A, B, ell_rate, rho_cat) {
    A <- as.data.frame(A)
    B <- as.data.frame(B)
    rate_diff <- outer(as.numeric(A$rate_mt), as.numeric(B$rate_mt), "-")
    K <- exp(-0.5 * (rate_diff / ell_rate)^2)
    K <- K * ifelse(outer(A$scenario, B$scenario, "=="), 1, rho_cat)
    K <- K * ifelse(outer(A$L_mid, B$L_mid, "=="), 1, rho_cat)
    K <- K * ifelse(outer(A$L_top, B$L_top, "=="), 1, rho_cat)
    K
  }

  tune_mixed_gp_kernel <- function(train_x, train_y_matrix) {
    grid <- expand.grid(
      ell_rate = c(0.04, 0.08, 0.16),
      rho_cat = c(0.30, 0.70),
      stringsAsFactors = FALSE
    )
    noise <- 1e-4
    scores <- rep(Inf, nrow(grid))
    n <- nrow(train_x)

    for (g in seq_len(nrow(grid))) {
      K <- mixed_kernel_matrix(train_x, train_x, grid$ell_rate[g], grid$rho_cat[g])
      diag(K) <- diag(K) + noise
      chol_K <- tryCatch(chol(K), error = function(e) NULL)
      if (is.null(chol_K)) next
      log_det <- 2 * sum(log(diag(chol_K)))
      score <- 0
      for (j in seq_len(ncol(train_y_matrix))) {
        y <- train_y_matrix[, j]
        alpha <- backsolve(chol_K, forwardsolve(t(chol_K), y))
        score <- score + 0.5 * sum(y * alpha) + 0.5 * log_det + 0.5 * n * log(2 * pi)
      }
      scores[g] <- score
    }
    grid[which.min(scores), , drop = FALSE]
  }

  fit_predict_mixed_gp <- function(train, test) {
    objective_names <- c(
      "f1_retained_integral_m3yr",
      "f2_unretained_final_m3",
      "f3_total_cost_eur"
    )
    x_train <- train[, c("scenario", "rate_mt", "L_mid", "L_top")]
    x_test <- test[, c("scenario", "rate_mt", "L_mid", "L_top")]
    y_raw <- as.matrix(train[, objective_names, drop = FALSE])
    y_mean <- colMeans(y_raw)
    y_sd <- apply(y_raw, 2, sd)
    y_sd[!is.finite(y_sd) | y_sd <= 0] <- 1
    y <- sweep(sweep(y_raw, 2, y_mean, "-"), 2, y_sd, "/")

    best <- tune_mixed_gp_kernel(x_train, y)
    K <- mixed_kernel_matrix(x_train, x_train, best$ell_rate[1], best$rho_cat[1])
    diag(K) <- diag(K) + 1e-4
    chol_K <- chol(K)
    K_star <- mixed_kernel_matrix(x_train, x_test, best$ell_rate[1], best$rho_cat[1])
    pred_std <- matrix(NA_real_, nrow(test), length(objective_names))
    for (j in seq_along(objective_names)) {
      alpha <- backsolve(chol_K, forwardsolve(t(chol_K), y[, j]))
      pred_std[, j] <- as.numeric(t(K_star) %*% alpha)
    }
    pred <- sweep(sweep(pred_std, 2, y_sd, "*"), 2, y_mean, "+")
    colnames(pred) <- objective_names
    list(prediction = pred, ell_rate = best$ell_rate[1], rho_cat = best$rho_cat[1])
  }

  fit_predict_random_forest <- function(train, test, trees, seed) {
    objective_names <- c(
      "f1_retained_integral_m3yr",
      "f2_unretained_final_m3",
      "f3_total_cost_eur"
    )
    train_model <- prepare_surrogate_predictors(train)
    test_model <- prepare_surrogate_predictors(test)
    pred <- matrix(NA_real_, nrow(test), length(objective_names), dimnames = list(NULL, objective_names))
    for (j in seq_along(objective_names)) {
      fit <- ranger::ranger(
        as.formula(paste(objective_names[j], "~ scenario + rate_mt + L_mid + L_top")),
        data = train_model,
        num.trees = trees,
        seed = seed + j,
        respect.unordered.factors = "order"
      )
      pred[, j] <- predict(fit, data = test_model)$predictions
    }
    pred
  }

  objective_names <- c(
    "f1_retained_integral_m3yr",
    "f2_unretained_final_m3",
    "f3_total_cost_eur"
  )
  validation_rows <- list()
  class_rows <- list()
  gp_rows <- list()
  selected_training_rows <- list()
  counter <- class_counter <- gp_counter <- selected_counter <- 0L

  for (split_seed in SURROGATE_SPLIT_SEEDS) {
    set.seed(split_seed)
    test_index <- unlist(
      lapply(split(seq_len(nrow(surrogate_data)), surrogate_data$scenario), function(index) {
        sample(index, max(1L, floor(length(index) * SURROGATE_TEST_FRACTION)))
      })
    )
    test <- surrogate_data[test_index, , drop = FALSE]
    train_pool <- surrogate_data[-test_index, , drop = FALSE]

    # The exact same training rows are supplied to MNL, GP, and RF.
    common_n <- min(SURROGATE_COMMON_TRAIN_N, nrow(train_pool))
    common_index <- sample(seq_len(nrow(train_pool)), common_n)
    train <- train_pool[common_index, , drop = FALSE]

    selected_counter <- selected_counter + 1L
    selected_training_rows[[selected_counter]] <- data.frame(
      split_seed = split_seed,
      train_row_index_in_pool = common_index,
      stringsAsFactors = FALSE
    )

    mnl_result <- fit_predict_mnl_validation(train, test, B = cfg$bins)
    gp_result <- fit_predict_mixed_gp(train, test)
    rf_prediction <- fit_predict_random_forest(train, test, SURROGATE_RF_TREES, split_seed)

    predictions <- list(
      MNL = mnl_result$prediction,
      MixedKernelGP = gp_result$prediction,
      RandomForest = rf_prediction
    )

    for (model_name in names(predictions)) {
      for (objective_name in objective_names) {
        mm <- regression_metrics(
          test[[objective_name]],
          predictions[[model_name]][, objective_name],
          sd(train[[objective_name]])
        )
        counter <- counter + 1L
        validation_rows[[counter]] <- data.frame(
          split_seed = split_seed,
          model = model_name,
          objective = objective_name,
          train_n = nrow(train),
          test_n = nrow(test),
          RMSE = unname(mm["RMSE"]),
          NRMSE_trainSD = unname(mm["NRMSE_trainSD"]),
          MAE = unname(mm["MAE"]),
          R2 = unname(mm["R2"]),
          stringsAsFactors = FALSE
        )
      }
    }

    class_counter <- class_counter + 1L
    class_rows[[class_counter]] <- data.frame(
      split_seed = split_seed,
      MNL_classification_accuracy = mnl_result$accuracy,
      MNL_multiclass_log_loss = mnl_result$log_loss,
      stringsAsFactors = FALSE
    )
    gp_counter <- gp_counter + 1L
    gp_rows[[gp_counter]] <- data.frame(
      split_seed = split_seed,
      GP_ell_rate = gp_result$ell_rate,
      GP_rho_categorical = gp_result$rho_cat,
      common_training_rows = nrow(train),
      stringsAsFactors = FALSE
    )
  }

  surrogate_validation_metrics <- do.call(rbind, validation_rows)
  surrogate_classification_metrics <- do.call(rbind, class_rows)
  gp_selected_parameters <- do.call(rbind, gp_rows)
  common_training_indices <- do.call(rbind, selected_training_rows)

  surrogate_validation_summary <- do.call(
    rbind,
    lapply(split(
      surrogate_validation_metrics,
      interaction(surrogate_validation_metrics$model, surrogate_validation_metrics$objective, drop = TRUE)
    ), function(d) {
      data.frame(
        model = d$model[1],
        objective = d$objective[1],
        train_n = unique(d$train_n)[1],
        test_n_mean = mean(d$test_n),
        RMSE_mean = mean(d$RMSE),
        RMSE_sd = sd(d$RMSE),
        NRMSE_trainSD_mean = mean(d$NRMSE_trainSD),
        NRMSE_trainSD_sd = sd(d$NRMSE_trainSD),
        MAE_mean = mean(d$MAE),
        MAE_sd = sd(d$MAE),
        R2_mean = mean(d$R2),
        R2_sd = sd(d$R2),
        stringsAsFactors = FALSE
      )
    })
  )

  # Paired tests across the five split seeds for each objective.
  surrogate_friedman <- do.call(rbind, lapply(objective_names, function(obj) {
    d <- surrogate_validation_metrics[surrogate_validation_metrics$objective == obj, ]
    wide <- reshape(
      d[, c("split_seed", "model", "NRMSE_trainSD")],
      idvar = "split_seed", timevar = "model", direction = "wide"
    )
    M <- as.matrix(wide[, paste0("NRMSE_trainSD.", c("MNL", "MixedKernelGP", "RandomForest"))])
    ft <- friedman.test(M)
    data.frame(
      objective = obj,
      statistic = unname(ft$statistic),
      df = unname(ft$parameter),
      p_value = ft$p.value,
      stringsAsFactors = FALSE
    )
  }))

  surrogate_pairwise <- do.call(rbind, lapply(objective_names, function(obj) {
    d <- surrogate_validation_metrics[surrogate_validation_metrics$objective == obj, ]
    pairs <- combn(c("MNL", "MixedKernelGP", "RandomForest"), 2, simplify = FALSE)
    out <- do.call(rbind, lapply(pairs, function(pair) {
      a <- d[d$model == pair[1], c("split_seed", "NRMSE_trainSD")]
      b <- d[d$model == pair[2], c("split_seed", "NRMSE_trainSD")]
      names(a)[2] <- "x"; names(b)[2] <- "y"
      z <- merge(a, b, by = "split_seed")
      wt <- wilcox.test(z$x, z$y, paired = TRUE, exact = FALSE)
      data.frame(
        objective = obj,
        model_1 = pair[1],
        model_2 = pair[2],
        mean_difference_model1_minus_model2 = mean(z$x - z$y),
        p_value_raw = wt$p.value,
        stringsAsFactors = FALSE
      )
    }))
    out$p_value_Holm <- p.adjust(out$p_value_raw, method = "holm")
    out
  }))

  write.csv(surrogate_data, file.path(OUTPUT_DIR, "surrogate_validation", "unique_design_dataset_used_for_validation.csv"), row.names = FALSE)
  write.csv(common_training_indices, file.path(OUTPUT_DIR, "surrogate_validation", "common_training_indices_by_split.csv"), row.names = FALSE)
  write.csv(surrogate_validation_metrics, file.path(OUTPUT_DIR, "surrogate_validation", "fair_heldout_surrogate_metrics_by_split.csv"), row.names = FALSE)
  write.csv(surrogate_validation_summary, file.path(OUTPUT_DIR, "surrogate_validation", "fair_heldout_surrogate_metrics_mean_sd.csv"), row.names = FALSE)
  write.csv(surrogate_friedman, file.path(OUTPUT_DIR, "surrogate_validation", "fair_surrogate_Friedman_tests.csv"), row.names = FALSE)
  write.csv(surrogate_pairwise, file.path(OUTPUT_DIR, "surrogate_validation", "fair_surrogate_pairwise_Wilcoxon_Holm.csv"), row.names = FALSE)
  write.csv(surrogate_classification_metrics, file.path(OUTPUT_DIR, "surrogate_validation", "MNL_classification_accuracy_logloss.csv"), row.names = FALSE)
  write.csv(gp_selected_parameters, file.path(OUTPUT_DIR, "surrogate_validation", "mixed_kernel_GP_selected_parameters.csv"), row.names = FALSE)

  p_surrogate <- ggplot2::ggplot(
    surrogate_validation_metrics,
    ggplot2::aes(x = model, y = NRMSE_trainSD, fill = model)
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.08, size = 1.8) +
    ggplot2::facet_wrap(~ objective, scales = "free_y") +
    ggplot2::labs(
      x = "Surrogate (same training rows)",
      y = "Held-out normalized RMSE"
    ) +
    theme_revision() +
    ggplot2::theme(legend.position = "none")

  ggplot2::ggsave(
    file.path(OUTPUT_DIR, "figures", "fair_heldout_surrogate_comparison_NRMSE.pdf"),
    p_surrogate, width = 9, height = 4.8
  )
  ggplot2::ggsave(
    file.path(OUTPUT_DIR, "figures", "fair_heldout_surrogate_comparison_NRMSE.png"),
    p_surrogate, width = 9, height = 4.8, dpi = 350
  )

  mnl_interpretation <- paste(
    "MNL was not selected because it produced the lowest held-out prediction",
    "error. The retrospective baseline comparison is reported transparently.",
    "MNL was retained as an exploratory probabilistic surrogate because it",
    "represents discrete joint objective states, accommodates the mixed",
    "categorical-continuous decision variables, and directly supplies class",
    "probabilities for acquisition sampling. Any lower predictive accuracy",
    "relative to the GP or random-forest baselines is acknowledged as a",
    "limitation rather than interpreted as empirical superiority."
  )
  writeLines(
    mnl_interpretation,
    file.path(OUTPUT_DIR, "text_for_manuscript", "MNL_justification_and_limitation.txt")
  )
}


# ============================================================
# 20) SIMULATOR REPEATABILITY, MASS-BALANCE DIAGNOSTICS,
#     AND OPTIONAL EXTERNAL BENCHMARK COMPARISON
# ============================================================

RUN_SIMULATOR_DIAGNOSTICS <- TRUE
SIMULATOR_REPEATABILITY_REPEATS <- 3L

run_julia_uncached_compact <- function(
    scenario,
    rate_mt,
    L_mid,
    L_top,
    cfg
) {
  inj <- build_injection_matrices(
    scenario = scenario,
    rate_mt = rate_mt,
    L_mid = L_mid,
    L_top = L_top,
    setup_result = setup_result,
    years = cfg$sim_years,
    rho_inj = cfg$rho_inj,
    availability = cfg$a_bar
  )

  sim_result <- julia_call(
    "run_simulation",
    start_time = 0.0,
    end_time = as.numeric(cfg$sim_years),
    time_step = 1.0,
    injection_rate_matrices = inj$injection_matrices,
    verbose = FALSE
  )

  if (!is.list(sim_result) ||
      is.null(sim_result$status) ||
      sim_result$status != "success") {
    stop("Uncached simulator diagnostic run failed.")
  }

  volume_series <- as.numeric(sim_result$total_co2_volumes)
  injected_total <- sum(
    vapply(inj$injection_matrices, sum, numeric(1))
  )
  final_retained <- tail(volume_series, 1)

  out <- data.frame(
    f1_retained_integral_m3yr = sum(volume_series),
    final_retained_m3 = final_retained,
    injected_total_m3 = injected_total,
    unretained_final_m3 = max(0, injected_total - final_retained),
    raw_mass_balance_difference_m3 = injected_total - final_retained,
    stringsAsFactors = FALSE
  )

  rm(inj, sim_result, volume_series)
  invisible(gc(verbose = FALSE))
  try(julia_command("GC.gc()"), silent = TRUE)
  out
}

if (isTRUE(RUN_SIMULATOR_DIAGNOSTICS)) {
  diagnostic_designs <- data.frame(
    case_id = c("low_one_well", "mid_one_well", "high_three_well"),
    scenario = c(1L, 1L, 3L),
    rate_mt = c(0.80, 0.90, 1.00),
    L_mid = c(4L, 5L, 6L),
    L_top = c(5L, 7L, 4L),
    stringsAsFactors = FALSE
  )

  repeatability_rows <- list()
  counter <- 0L

  for (i in seq_len(nrow(diagnostic_designs))) {
    for (repeat_id in seq_len(SIMULATOR_REPEATABILITY_REPEATS)) {
      result_i <- run_julia_uncached_compact(
        diagnostic_designs$scenario[i],
        diagnostic_designs$rate_mt[i],
        diagnostic_designs$L_mid[i],
        diagnostic_designs$L_top[i],
        cfg
      )
      counter <- counter + 1L
      repeatability_rows[[counter]] <- cbind(
        diagnostic_designs[i, , drop = FALSE],
        repeat_id = repeat_id,
        result_i
      )
    }
  }

  simulator_repeatability <- do.call(rbind, repeatability_rows)

  repeatability_summary <- do.call(
    rbind,
    lapply(split(simulator_repeatability, simulator_repeatability$case_id),
           function(d) {
             data.frame(
               case_id = d$case_id[1],
               f1_range = diff(range(d$f1_retained_integral_m3yr)),
               final_retained_range = diff(range(d$final_retained_m3)),
               unretained_range = diff(range(d$unretained_final_m3)),
               maximum_absolute_repeat_difference = max(
                 c(
                   diff(range(d$f1_retained_integral_m3yr)),
                   diff(range(d$final_retained_m3)),
                   diff(range(d$unretained_final_m3))
                 )
               ),
               stringsAsFactors = FALSE
             )
           })
  )

  write.csv(
    simulator_repeatability,
    file.path(
      OUTPUT_DIR,
      "simulator_validation",
      "simulator_repeatability_raw.csv"
    ),
    row.names = FALSE
  )
  write.csv(
    repeatability_summary,
    file.path(
      OUTPUT_DIR,
      "simulator_validation",
      "simulator_repeatability_summary.csv"
    ),
    row.names = FALSE
  )

  # External benchmark interface.
  benchmark_template_path <- file.path(
    OUTPUT_DIR,
    "simulator_validation",
    "sleipner_external_benchmark_template.csv"
  )

  if (!file.exists(benchmark_template_path)) {
    benchmark_template <- data.frame(
      case_id = character(0),
      scenario = integer(0),
      rate_mt = numeric(0),
      L_mid = integer(0),
      L_top = integer(0),
      reference_f1_retained_integral_m3yr = numeric(0),
      reference_final_retained_m3 = numeric(0),
      reference_unretained_final_m3 = numeric(0),
      reference_source = character(0),
      stringsAsFactors = FALSE
    )
    write.csv(
      benchmark_template,
      benchmark_template_path,
      row.names = FALSE
    )
  }

  benchmark_input_path <- file.path(
    PROJECT_ROOT,
    "data",
    "sleipner_external_benchmark.csv"
  )

  if (file.exists(benchmark_input_path)) {
    benchmark_data <- read.csv(
      benchmark_input_path,
      stringsAsFactors = FALSE
    )

    required_benchmark_columns <- c(
      "case_id",
      "scenario",
      "rate_mt",
      "L_mid",
      "L_top",
      "reference_f1_retained_integral_m3yr",
      "reference_final_retained_m3",
      "reference_unretained_final_m3",
      "reference_source"
    )

    missing_columns <- setdiff(
      required_benchmark_columns,
      names(benchmark_data)
    )
    if (length(missing_columns) > 0L) {
      stop(
        "External benchmark CSV is missing columns: ",
        paste(missing_columns, collapse = ", ")
      )
    }

    model_rows <- lapply(seq_len(nrow(benchmark_data)), function(i) {
      rr <- run_julia_sim_cached(
        benchmark_data$scenario[i],
        benchmark_data$rate_mt[i],
        benchmark_data$L_mid[i],
        benchmark_data$L_top[i],
        cfg
      )
      volume_series <- rr$total_co2_volumes
      final_retained <- tail(volume_series, 1)

      data.frame(
        case_id = benchmark_data$case_id[i],
        model_f1_retained_integral_m3yr = sum(volume_series),
        model_final_retained_m3 = final_retained,
        model_unretained_final_m3 = max(
          0,
          rr$injected_total - final_retained
        ),
        stringsAsFactors = FALSE
      )
    })

    benchmark_comparison <- merge(
      benchmark_data,
      do.call(rbind, model_rows),
      by = "case_id"
    )

    benchmark_long <- rbind(
      data.frame(
        metric = "f1_retained_integral_m3yr",
        observed = benchmark_comparison$reference_f1_retained_integral_m3yr,
        predicted = benchmark_comparison$model_f1_retained_integral_m3yr
      ),
      data.frame(
        metric = "final_retained_m3",
        observed = benchmark_comparison$reference_final_retained_m3,
        predicted = benchmark_comparison$model_final_retained_m3
      ),
      data.frame(
        metric = "unretained_final_m3",
        observed = benchmark_comparison$reference_unretained_final_m3,
        predicted = benchmark_comparison$model_unretained_final_m3
      )
    )

    benchmark_metrics <- do.call(
      rbind,
      lapply(split(benchmark_long, benchmark_long$metric), function(d) {
        mm <- regression_metrics(
          d$observed,
          d$predicted,
          sd(d$observed)
        )
        data.frame(
          metric = d$metric[1],
          RMSE = unname(mm["RMSE"]),
          NRMSE_referenceSD = unname(mm["NRMSE_trainSD"]),
          MAE = unname(mm["MAE"]),
          R2 = unname(mm["R2"]),
          stringsAsFactors = FALSE
        )
      })
    )

    write.csv(
      benchmark_comparison,
      file.path(
        OUTPUT_DIR,
        "simulator_validation",
        "external_benchmark_case_comparison.csv"
      ),
      row.names = FALSE
    )
    write.csv(
      benchmark_metrics,
      file.path(
        OUTPUT_DIR,
        "simulator_validation",
        "external_benchmark_error_metrics.csv"
      ),
      row.names = FALSE
    )

    p_benchmark <- ggplot2::ggplot(
      benchmark_long,
      ggplot2::aes(x = observed, y = predicted)
    ) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
      ggplot2::facet_wrap(~ metric, scales = "free") +
      ggplot2::labs(
        x = "External benchmark/reference value",
        y = "Simplified simulator value"
      ) +
      theme_revision()

    ggplot2::ggsave(
      file.path(
        OUTPUT_DIR,
        "figures",
        "external_simulator_benchmark_comparison.pdf"
      ),
      p_benchmark,
      width = 9,
      height = 4.5
    )
  } else {
    limitation_text <- paste(
      "The numerical model is a simplified, deterministic,",
      "Sleipner-inspired test environment used to evaluate the",
      "optimization workflow. It was configured using representative",
      "reservoir geometry and properties but was not history matched to",
      "Sleipner field observations and was not validated as a predictive",
      "field-scale model. Consequently, the reported objective values",
      "should be interpreted as internally consistent numerical benchmark",
      "outputs rather than field forecasts. Geological uncertainty,",
      "multiphase-pressure calibration, geomechanical effects, and",
      "ensemble-based uncertainty propagation remain outside the scope of",
      "the present study and constitute priorities for future work."
    )

    writeLines(
      limitation_text,
      file.path(
        OUTPUT_DIR,
        "text_for_manuscript",
        "simulator_validation_limitation_paragraph.txt"
      )
    )

    writeLines(
      c(
        "No external benchmark CSV was found.",
        paste0(
          "To perform quantitative simulator validation, place a completed ",
          "'sleipner_external_benchmark.csv' in the Julia project directory."
        ),
        paste0(
          "A column template is available at: ",
          benchmark_template_path
        ),
        "",
        "Until reference data are supplied, use the generated limitation paragraph and do not claim field validation."
      ),
      file.path(
        OUTPUT_DIR,
        "simulator_validation",
        "external_benchmark_status.txt"
      )
    )
  }
}


# ============================================================
# 21) COST METADATA AND ENHANCED POST-HOC COST SENSITIVITY
#     2024 EUR per metric tonne of CO2
# ============================================================

RUN_ENHANCED_COST_SENSITIVITY <- TRUE
COST_CURRENCY <- "EUR"
COST_PRICE_YEAR <- 2024L
COST_SOURCE <- paste(
  "European Commission Joint Research Centre (2024),",
  "Clean Energy Technology Observatory: Carbon Capture, Utilisation and",
  "Storage in the European Union, JRC139285. Reported ranges:",
  "capture 40-90 EUR/tCO2, transport 2-30 EUR/tCO2, storage 5-35 EUR/tCO2."
)

cost_coefficient_metadata <- data.frame(
  parameter = c(
    "CO2 capture unit cost",
    "CO2 transport unit cost",
    "CO2 storage unit cost",
    "Total representative chain cost"
  ),
  R_symbol = c("c_cap", "c_trans", "c_storage", "c_total"),
  numerical_value = c(
    cfg$c_cap,
    cfg$c_trans,
    cfg$c_storage,
    cfg$c_cap + cfg$c_trans + cfg$c_storage
  ),
  unit = rep("2024 EUR per metric tonne CO2", 4),
  currency = rep(COST_CURRENCY, 4),
  price_year = rep(COST_PRICE_YEAR, 4),
  source_or_derivation = c(
    "Midpoint of JRC capture range 40-90 EUR/tCO2",
    "Midpoint of JRC transport range 2-30 EUR/tCO2",
    "Midpoint of JRC storage range 5-35 EUR/tCO2",
    "Sum of the three representative midpoint values"
  ),
  bibliographic_source = rep(COST_SOURCE, 4),
  stringsAsFactors = FALSE
)
write.csv(
  cost_coefficient_metadata,
  file.path(OUTPUT_DIR, "cost_sensitivity", "cost_coefficients_2024_EUR_per_tCO2.csv"),
  row.names = FALSE
)

if (isTRUE(RUN_ENHANCED_COST_SENSITIVITY)) {
  enhanced_cost_scenarios <- data.frame(
    scenario_name = c(
      "base",
      "all_chain_costs_minus25pct",
      "all_chain_costs_plus25pct",
      "capture_low_JRC_40",
      "capture_high_JRC_90",
      "transport_low_JRC_2",
      "transport_high_JRC_30",
      "storage_low_JRC_5",
      "storage_high_JRC_35"
    ),
    capture_eur_per_t = c(65, 48.75, 81.25, 40, 90, 65, 65, 65, 65),
    transport_eur_per_t = c(16, 12, 20, 16, 16, 2, 30, 16, 16),
    storage_eur_per_t = c(20, 15, 25, 20, 20, 20, 20, 5, 35),
    stringsAsFactors = FALSE
  )
  enhanced_cost_scenarios$total_chain_eur_per_t <- rowSums(
    enhanced_cost_scenarios[, c("capture_eur_per_t", "transport_eur_per_t", "storage_eur_per_t")]
  )

  calculate_injected_total_m3 <- function(d, cfg) {
    n_wells <- ifelse(d$scenario == 3L, 3, 1)
    rate_per_well <- mt_to_m3yr(d$rate_mt, cfg$rho_inj) * cfg$a_bar
    # All three five-year phases remain active, including selections 8 and 9.
    n_wells * rate_per_well * cfg$sim_years
  }

  history_cost_base <- all_histories
  history_cost_base$injected_total_m3 <- calculate_injected_total_m3(history_cost_base, cfg)
  history_cost_base$injected_total_tonnes <- (
    history_cost_base$injected_total_m3 * cfg$rho_inj / 1000.0
  )

  enhanced_cost_data <- do.call(rbind, lapply(seq_len(nrow(enhanced_cost_scenarios)), function(i) {
    sc <- enhanced_cost_scenarios[i, ]
    d <- history_cost_base
    d$cost_scenario <- sc$scenario_name
    d$f3_sensitivity_cost_eur <- (
      sc$total_chain_eur_per_t * d$injected_total_tonnes
    )
    d
  }))

  cost_metric_rows <- list()
  cost_front_rows <- list()
  cost_counter <- front_counter <- 0L

  for (cost_scenario_name in enhanced_cost_scenarios$scenario_name) {
    ds <- enhanced_cost_data[enhanced_cost_data$cost_scenario == cost_scenario_name, , drop = FALSE]
    Y_all_raw <- as.matrix(ds[, c(
      "f1_retained_integral_m3yr",
      "f2_unretained_final_m3",
      "f3_sensitivity_cost_eur"
    )])
    Y_all_max <- cbind(Y_all_raw[, 1], -Y_all_raw[, 2], -Y_all_raw[, 3])
    lo_s <- apply(Y_all_max, 2, min)
    hi_s <- apply(Y_all_max, 2, max)
    span_s <- hi_s - lo_s
    span_s[span_s <= 0 | !is.finite(span_s)] <- 1
    hi_s <- lo_s + span_s
    Y_all_norm <- normalize_to_unit(Y_all_max, lo_s, hi_s)
    reference_s <- unique(nondominated(Y_all_norm))

    for (seed_value in SEEDS) {
      for (method_value in METHODS) {
        dsm <- ds[ds$seed == seed_value & ds$method == method_value, , drop = FALSE]
        Y_raw <- as.matrix(dsm[, c(
          "f1_retained_integral_m3yr",
          "f2_unretained_final_m3",
          "f3_sensitivity_cost_eur"
        )])
        Y_max <- cbind(Y_raw[, 1], -Y_raw[, 2], -Y_raw[, 3])
        Y_norm <- normalize_to_unit(Y_max, lo_s, hi_s)
        front <- unique(nondominated(Y_norm))
        cost_counter <- cost_counter + 1L
        cost_metric_rows[[cost_counter]] <- data.frame(
          cost_scenario = cost_scenario_name,
          seed = seed_value,
          method = method_value,
          final_hypervolume = compute_hv_exact_3d(front, evaluation_ref_point),
          final_IGD_plus = compute_igd_plus_max(front, reference_s),
          pareto_size = nrow(front),
          minimum_cost_eur = min(Y_raw[, 3]),
          maximum_storage_m3yr = max(Y_raw[, 1]),
          minimum_unretained_m3 = min(Y_raw[, 2]),
          stringsAsFactors = FALSE
        )
      }
    }

    for (method_value in METHODS) {
      dm <- ds[ds$method == method_value, , drop = FALSE]
      Y_raw <- as.matrix(dm[, c(
        "f1_retained_integral_m3yr",
        "f2_unretained_final_m3",
        "f3_sensitivity_cost_eur"
      )])
      Y_max <- cbind(Y_raw[, 1], -Y_raw[, 2], -Y_raw[, 3])
      Y_norm <- normalize_to_unit(Y_max, lo_s, hi_s)
      nd <- nondominated(Y_norm)
      design_keys <- paste(dm$scenario, round(dm$rate_mt, 6), dm$L_mid, dm$L_top, sep = "|")
      norm_keys <- apply(round(Y_norm, 12), 1, paste, collapse = "|")
      nd_keys <- unique(apply(round(nd, 12), 1, paste, collapse = "|"))
      selected <- unique(design_keys[norm_keys %in% nd_keys])
      front_counter <- front_counter + 1L
      cost_front_rows[[front_counter]] <- data.frame(
        cost_scenario = cost_scenario_name,
        method = method_value,
        design_key = selected,
        stringsAsFactors = FALSE
      )
    }
  }

  enhanced_cost_metrics <- do.call(rbind, cost_metric_rows)
  enhanced_cost_front_keys <- do.call(rbind, cost_front_rows)
  enhanced_cost_summary <- do.call(rbind, lapply(split(
    enhanced_cost_metrics,
    interaction(enhanced_cost_metrics$cost_scenario, enhanced_cost_metrics$method, drop = TRUE)
  ), function(d) {
    data.frame(
      cost_scenario = d$cost_scenario[1],
      method = d$method[1],
      n_seeds = nrow(d),
      HV_mean = mean(d$final_hypervolume),
      HV_sd = sd(d$final_hypervolume),
      IGDplus_mean = mean(d$final_IGD_plus),
      IGDplus_sd = sd(d$final_IGD_plus),
      mean_pareto_size = mean(d$pareto_size),
      mean_minimum_cost_eur = mean(d$minimum_cost_eur),
      stringsAsFactors = FALSE
    )
  }))

  base_keys_df <- enhanced_cost_front_keys[enhanced_cost_front_keys$cost_scenario == "base", ]
  stability_rows <- list(); stability_counter <- 0L
  for (scenario_name in enhanced_cost_scenarios$scenario_name) {
    for (method_value in METHODS) {
      base_keys <- unique(base_keys_df$design_key[base_keys_df$method == method_value])
      scenario_keys <- unique(enhanced_cost_front_keys$design_key[
        enhanced_cost_front_keys$cost_scenario == scenario_name &
        enhanced_cost_front_keys$method == method_value
      ])
      union_n <- length(union(base_keys, scenario_keys))
      stability_counter <- stability_counter + 1L
      stability_rows[[stability_counter]] <- data.frame(
        cost_scenario = scenario_name,
        method = method_value,
        base_front_size = length(base_keys),
        scenario_front_size = length(scenario_keys),
        Pareto_design_Jaccard_vs_base = if (union_n > 0) length(intersect(base_keys, scenario_keys)) / union_n else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  cost_front_stability <- do.call(rbind, stability_rows)

  cost_method_ranking <- do.call(rbind, lapply(split(enhanced_cost_summary, enhanced_cost_summary$cost_scenario), function(d) {
    d$HV_rank <- rank(-d$HV_mean, ties.method = "average")
    d$IGDplus_rank <- rank(d$IGDplus_mean, ties.method = "average")
    d
  }))

  write.csv(enhanced_cost_scenarios, file.path(OUTPUT_DIR, "cost_sensitivity", "JRC_range_cost_scenarios_2024_EUR.csv"), row.names = FALSE)
  write.csv(enhanced_cost_metrics, file.path(OUTPUT_DIR, "cost_sensitivity", "cost_sensitivity_metrics_by_seed_method.csv"), row.names = FALSE)
  write.csv(enhanced_cost_summary, file.path(OUTPUT_DIR, "cost_sensitivity", "cost_sensitivity_metrics_mean_sd.csv"), row.names = FALSE)
  write.csv(cost_method_ranking, file.path(OUTPUT_DIR, "cost_sensitivity", "method_ranking_under_cost_scenarios.csv"), row.names = FALSE)
  write.csv(cost_front_stability, file.path(OUTPUT_DIR, "cost_sensitivity", "Pareto_design_stability_under_cost_scenarios.csv"), row.names = FALSE)

  p_cost_hv <- ggplot2::ggplot(
    enhanced_cost_summary,
    ggplot2::aes(x = cost_scenario, y = HV_mean, group = method, color = method)
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(values = method_cols) +
    ggplot2::labs(x = "2024-EUR/tCO2 cost scenario", y = "Mean final hypervolume") +
    theme_revision() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(
    file.path(OUTPUT_DIR, "figures", "cost_sensitivity_method_HV.pdf"),
    p_cost_hv, width = 10, height = 5.5
  )

  cost_limitation_text <- paste(
    "The economic objective is a transparent representative cost index in",
    "2024 EUR per metric tonne of injected CO2. The capture, transport, and",
    "storage coefficients are the midpoints of the ranges reported in the",
    "European Commission JRC 2024 CCUS status report. The sensitivity analysis",
    "recalculates objective values and Pareto metrics for all evaluated designs",
    "across the full published ranges. It is a post-hoc robustness analysis",
    "conditional on the sampled designs, not a complete BO re-optimization",
    "under every economic scenario."
  )
  writeLines(cost_limitation_text, file.path(OUTPUT_DIR, "text_for_manuscript", "cost_model_and_sensitivity_scope.txt"))
}

reviewer_coverage_manifest <- data.frame(
  reviewer_concern = c(
    "Repeated stochastic BO runs",
    "Statistical comparison of acquisitions",
    "Independent reference Pareto front",
    "Standard EHVI implementation",
    "Standard EPI implementation",
    "SCAL-UCB kappa sensitivity",
    "EHVI predictive-sample sensitivity",
    "Held-out MNL validation",
    "MNL versus mixed-kernel GP",
    "MNL versus random forest",
    "Simulator numerical repeatability",
    "External simulator benchmark",
    "Cost coefficient metadata",
    "Cost sensitivity of Pareto metrics",
    "Geological uncertainty"
  ),
  code_or_output = c(
    "10 paired seeds; shared 20-point initial designs",
    "Friedman and paired Wilcoxon-Holm tables",
    "grid-only primary front plus grid+BO union supplementary front",
    "Sample mean of nonnegative HV improvement; extra coefficient = 0",
    "Sample mean of nonnegative preference improvement; extra coefficient = 0",
    "kappa = 0.5, 1, 2, 4",
    "predictive samples = 10, 20, 40",
    "five grouped splits with the same 180 training designs for every surrogate",
    "custom continuous-categorical product kernel GP",
    "ranger regression baseline",
    "three repeated uncached runs for three fixed designs",
    "optional CSV interface; requires external benchmark data",
    "completed 2024-EUR/tCO2 table based on JRC 2024 published ranges",
    "HV, IGD+, Pareto size, ranking, and front stability under scenarios",
    "manuscript limitation only; no geological ensembles are available"
  ),
  status_after_running = c(
    "addressed",
    "addressed",
    "addressed",
    "addressed",
    "addressed",
    "addressed",
    "addressed",
    "addressed",
    "addressed",
    "addressed",
    "addressed",
    "key limitation acknowledged; optional external data interface retained",
    "addressed with EUR, 2024 price basis, units, and JRC source",
    "addressed as post-hoc robustness analysis",
    "acknowledge as limitation and future ensemble extension"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  reviewer_coverage_manifest,
  file.path(
    OUTPUT_DIR,
    "tables",
    "reviewer_concern_analysis_manifest.csv"
  ),
  row.names = FALSE
)

writeLines(
  c(
    "SIMULATOR VALIDATION LIMITATION:",
    "The deterministic repeatability checks establish numerical consistency but do not constitute physical validation. The Sleipner-inspired simulator has not been history matched against field observations and its outputs are interpreted as numerical benchmark responses rather than field forecasts.",
    "",
    "GEOLOGICAL UNCERTAINTY LIMITATION:",
    "The present analysis uses one deterministic geological realization. Uncertainty in facies architecture, permeability, porosity, capillary entry pressure, and seal continuity is therefore not propagated through the optimization. Future work should evaluate each candidate over an ensemble of geological realizations and optimize expected performance together with downside-risk measures.",
    "",
    "LAYER INTERPRETATION:",
    sleipner_layer_interpretation
  ),
  file.path(OUTPUT_DIR, "text_for_manuscript", "consolidated_limitations_and_layer_interpretation.txt")
)

# ------------------------------------------------------------
# 19) Save the complete R session objects and session information
# ------------------------------------------------------------
saveRDS(
  list(
    cfg = cfg,
    seeds = SEEDS,
    methods = METHODS,
    multi_results = multi_results,
    global_bounds = global_bounds,
    reference_front_norm = reference_front_norm,
    metric_trace = metric_trace,
    final_metrics = final_metrics,
    final_summary = final_summary,
    friedman_results = friedman_results,
    wilcoxon_results = wilcoxon_results,
    surrogate_validation_summary = if (exists("surrogate_validation_summary")) surrogate_validation_summary else NULL,
    surrogate_friedman = if (exists("surrogate_friedman")) surrogate_friedman else NULL,
    surrogate_pairwise = if (exists("surrogate_pairwise")) surrogate_pairwise else NULL,
    igd_reference_robustness = igd_reference_robustness,
    simulator_repeatability = if (exists("simulator_repeatability")) simulator_repeatability else NULL,
    enhanced_cost_summary = if (exists("enhanced_cost_summary")) enhanced_cost_summary else NULL,
    cost_front_stability = if (exists("cost_front_stability")) cost_front_stability else NULL,
    reviewer_coverage_manifest = reviewer_coverage_manifest
  ),
  file.path(OUTPUT_DIR, "complete_reviewer_revision_analysis.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "sessionInfo.txt")
)

message("\n============================================================")
message("MULTI-SEED REVIEWER ANALYSIS COMPLETED")
message("Results folder: ", normalizePath(OUTPUT_DIR, winslash = "/"))
message("Primary outputs:")
message("  tables/final_performance_mean_sd_95CI.csv")
message("  tables/friedman_tests_all_methods.csv")
message("  tables/paired_wilcoxon_Holm_tests.csv")
message("  tables/HV_convergence_mean_sd_95CI.csv")
message("  tables/IGDplus_convergence_mean_sd_95CI.csv")
message("  figures/multiseed_HV_convergence_95CI.pdf")
message("  figures/multiseed_IGDplus_convergence_95CI.pdf")
message("  surrogate_validation/fair_heldout_surrogate_metrics_mean_sd.csv")
message("  simulator_validation/simulator_repeatability_summary.csv")
message("  cost_sensitivity/cost_sensitivity_metrics_mean_sd.csv")
message("  tables/hyperparameter_sensitivity_mean_sd_95CI.csv")
message("  reference/grid_only_reference_pareto_front_normalized.csv")
message("============================================================")

# ============================================================
# Generate the final pooled Pareto trade-off figure
# ============================================================

options(stringsAsFactors = FALSE)

required_packages <- c("ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the missing package(s) before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 1) Paths
# ------------------------------------------------------------
RESULT_DIR <- RESULT_ROOT

INPUT_CSV <- file.path(
  RESULT_DIR,
  "tables",
  "all_evaluated_designs_all_seeds.csv"
)

if (!file.exists(INPUT_CSV)) {
  stop(
    "Could not find the final archived result table at:\n",
    INPUT_CSV,
    "\nRun this script from the directory containing the final output folder, ",
    "or edit RESULT_DIR near the top of this file."
  )
}

FIGURE_DIR <- file.path(RESULT_ROOT, "figures")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

OUTPUT_PNG <- file.path(
  FIGURE_DIR,
  "Figure_final_Pareto_tradeoffs.png"
)
OUTPUT_PDF <- file.path(
  FIGURE_DIR,
  "Figure_final_Pareto_tradeoffs.pdf"
)
OUTPUT_CSV <- file.path(
  RESULT_ROOT,
  "tables",
  "pooled_method_specific_pareto_designs_for_figure.csv"
)

# ------------------------------------------------------------
# 2) Read and validate the final data
# ------------------------------------------------------------
all_designs <- read.csv(INPUT_CSV, check.names = FALSE)

required_columns <- c(
  "seed",
  "method",
  "scenario",
  "rate_mt",
  "L_mid",
  "L_top",
  "f1_retained_integral_m3yr",
  "f2_unretained_final_m3",
  "f3_total_cost_eur"
)
missing_columns <- setdiff(required_columns, names(all_designs))
if (length(missing_columns) > 0L) {
  stop(
    "The archived table is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (nrow(all_designs) != 4000L) {
  warning(
    "The table contains ", nrow(all_designs),
    " rows rather than the expected 4,000 primary evaluations."
  )
}

all_designs$method <- factor(
  all_designs$method,
  levels = c("EHVI", "EPI", "SCAL", "TS")
)

# ------------------------------------------------------------
# 3) Exact pooled method-specific non-dominated set
# ------------------------------------------------------------
# Objectives are transformed to a common maximization convention:
#   maximize f1, maximize -f2, maximize -f3.
non_dominated_mask <- function(Y) {
  Y <- as.matrix(Y)
  n <- nrow(Y)
  keep <- rep(TRUE, n)
  
  for (i in seq_len(n)) {
    dominates_i <- apply(
      Y,
      1,
      function(candidate) {
        all(candidate >= Y[i, ]) && any(candidate > Y[i, ])
      }
    )
    dominates_i[i] <- FALSE
    if (any(dominates_i)) keep[i] <- FALSE
  }
  
  keep
}

method_levels <- levels(all_designs$method)
pareto_parts <- vector("list", length(method_levels))

for (j in seq_along(method_levels)) {
  method_name <- method_levels[j]
  d <- all_designs[all_designs$method == method_name, , drop = FALSE]
  
  # of each objective-identical design before extracting the pooled front.
  unique_key <- paste(
    d$scenario,
    sprintf("%.10f", d$rate_mt),
    d$L_mid,
    d$L_top,
    sprintf("%.10f", d$f1_retained_integral_m3yr),
    sprintf("%.10f", d$f2_unretained_final_m3),
    sprintf("%.10f", d$f3_total_cost_eur),
    sep = "|"
  )
  d <- d[!duplicated(unique_key), , drop = FALSE]
  
  Ymax <- cbind(
    d$f1_retained_integral_m3yr,
    -d$f2_unretained_final_m3,
    -d$f3_total_cost_eur
  )
  d$pooled_pareto <- non_dominated_mask(Ymax)
  pareto_parts[[j]] <- d
}

plot_data <- do.call(rbind, pareto_parts)
plot_data$method_label <- factor(
  ifelse(plot_data$method == "SCAL", "SCAL-UCB", as.character(plot_data$method)),
  levels = c("EHVI", "EPI", "SCAL-UCB", "TS")
)

write.csv(plot_data, OUTPUT_CSV, row.names = FALSE)

# ------------------------------------------------------------
# 4) Publication styling
# ------------------------------------------------------------
method_colors <- c(
  EHVI = "#8E44AD",
  EPI = "#E31A1C",
  `SCAL-UCB` = "#F28E2B",
  TS = "#00A6D6"
)

base_theme <- ggplot2::theme_minimal(base_size = 11.5) +
  ggplot2::theme(
    legend.position = "top",
    legend.title = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(
      size = 11.5,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = ggplot2::element_text(size = 10.5)
  )

make_panel <- function(
    x_column,
    y_column,
    x_divisor,
    y_divisor,
    x_label,
    y_label,
    panel_title,
    show_legend = FALSE
) {
  d <- plot_data
  d$x_scaled <- d[[x_column]] / x_divisor
  d$y_scaled <- d[[y_column]] / y_divisor
  
  ggplot2::ggplot(d, ggplot2::aes(x = x_scaled, y = y_scaled)) +
    ggplot2::geom_point(
      ggplot2::aes(color = method_label),
      size = 0.9,
      alpha = 0.09,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = d[d$pooled_pareto, , drop = FALSE],
      ggplot2::aes(fill = method_label),
      color = "black",
      size = 2.4,
      alpha = 0.95,
      shape = 21,
      stroke = 0.25,
      show.legend = show_legend
    ) +
    ggplot2::scale_color_manual(values = method_colors, drop = FALSE) +
    ggplot2::scale_fill_manual(values = method_colors, drop = FALSE) +
    ggplot2::labs(
      title = panel_title,
      x = x_label,
      y = y_label
    ) +
    base_theme +
    ggplot2::theme(
      legend.position = if (show_legend) "top" else "none"
    )
}

p1 <- make_panel(
  "f1_retained_integral_m3yr",
  "f2_unretained_final_m3",
  1e8,
  1e7,
  expression("Time-integrated retained volume (" * 10^8 * " m"^3 * " yr)"),
  expression("Final-time unretained volume (" * 10^7 * " m"^3 * ")"),
  "(a) Retained volume versus unretained volume",
  show_legend = TRUE
)

p2 <- make_panel(
  "f1_retained_integral_m3yr",
  "f3_total_cost_eur",
  1e8,
  1e9,
  expression("Time-integrated retained volume (" * 10^8 * " m"^3 * " yr)"),
  expression("Full-chain cost (" * 10^9 * " EUR)"),
  "(b) Retained volume versus cost"
)

p3 <- make_panel(
  "f2_unretained_final_m3",
  "f3_total_cost_eur",
  1e7,
  1e9,
  expression("Final-time unretained volume (" * 10^7 * " m"^3 * ")"),
  expression("Full-chain cost (" * 10^9 * " EUR)"),
  "(c) Unretained volume versus cost"
)

# ------------------------------------------------------------
# 5) Combine the three panels using the base grid package
# ------------------------------------------------------------
draw_three_panels <- function() {
  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 1,
    ncol = 3,
    widths = grid::unit(c(1, 1, 1), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  
  print(
    p1,
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1),
    newpage = FALSE
  )
  print(
    p2,
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2),
    newpage = FALSE
  )
  print(
    p3,
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 3),
    newpage = FALSE
  )
}

png(
  filename = OUTPUT_PNG,
  width = 3900,
  height = 1250,
  res = 300,
  bg = "white"
)
draw_three_panels()
dev.off()

pdf(
  file = OUTPUT_PDF,
  width = 13,
  height = 4.2,
  useDingbats = FALSE
)
draw_three_panels()
dev.off()

cat("Created:\n")
cat("  ", OUTPUT_PNG, "\n", sep = "")
cat("  ", OUTPUT_PDF, "\n", sep = "")
cat("  ", OUTPUT_CSV, "\n", sep = "")
cat("Pooled Pareto counts by method:\n")
print(table(plot_data$method_label[plot_data$pooled_pareto]))

