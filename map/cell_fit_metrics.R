#!/usr/bin/env Rscript
# Per-cell fit metrics, replicating the core math of puck_cb_promrank (gamma) and
# puck_c1_intden (radial density) from the aBO135 tags_puck_plots functions — but
# returning numbers instead of plots, and merging them back into cell_coords.csv.
#
# Columns added (in place) to cell_coords.csv:
#   n_peaks           : number of finite, positive prominence peaks for the cell.
#   -- single gamma (original) --
#   gamma_resid_mad   : MAD of background residuals around the gamma 1:1 line
#                       (rank > 2, prom in [prom_lo, bg_prom_max]); log10 space.
#   p1_z, p2_z        : (log10 obs - log10 gamma-expected) / gamma_resid_mad for
#                       peaks 1 and 2.
#   -- (1) overall best family, lowest BIC across ALL candidates --
#   best_family       : winner among {gamma, lognormal, powerlaw, gamma_mix_k2/3,
#                       lnorm_mix_k2/3}. powerlaw = Pareto type I (straight line in
#                       log-log rank). Mixtures fit only when bg n >= 30.
#                       DEGENERACY GUARD: a candidate whose order-stat residual MAD is
#                       >3x the median candidate MAD (or non-finite) is dropped before
#                       choosing — stops a heavy concentrated mixture component from
#                       winning BIC with a pathological expected axis.
#   best_resid_mad    : MAD of bg residuals around the chosen family's 1:1 line.
#   best_p1_z, best_p2_z : peak z-scores under the chosen family.
#   -- (2) gamma-family mixture, k auto-selected by BIC (single gamma + k=2/3) --
#   mix_k             : gamma components chosen (1 = single gamma).
#   mix_resid_mad     : MAD of bg residuals around the chosen-mixture order-stat line.
#   mix_p1_z, mix_p2_z : peak z-scores under the chosen gamma mixture.
#   mix_vs_gamma_dbic : BIC(single gamma) - BIC(chosen mixture); >0 => mixture preferred
#                       (0 when k=1 is selected).
#   -- (3) parametric-free empirical background --
#   emp_bg_mad        : MAD of log10(bg prominences) themselves (rank>2, prom in
#                       [prom_lo, bg_prom_max]); no distributional assumption.
#   prom1_bg_z, prom2_bg_z : (log10 prom_k - median(log10 bg)) / emp_bg_mad.
#   -- intden radial density --
#   intden_peak_ratio : strongest density (1st-deriv) local max within 50um
#                       divided by the median of the outer local maxima.
#   intden_peak_z     : (peak - median(rest)) / MAD(rest) of the outer maxima.
#
# Usage: Rscript cell_fit_metrics.R <cell_coords.csv>   # merged in place
#        Rscript cell_fit_metrics.R <cell_coords.csv> <out.csv>   # write elsewhere
#
# NOTE: fit constants below must stay in sync with the plotting functions in
# scripts/utils/rmd/tags_puck_plots.Rmd (smooth_spar / bg_prom_max / split_r /
# bin_width). This script is the single source of truth for the *pipeline* metrics.

suppressMessages(library(MASS))
have_mixtools <- requireNamespace("mixtools", quietly = TRUE)

args    <- commandArgs(trailingOnly = TRUE)
map_csv <- args[1]
out_csv <- if (length(args) >= 2) args[2] else map_csv   # default: merge in place

bg_prom_max <- 0.1
prom_lo     <- 1e-5
bin_width   <- 3
smooth_spar <- 0.55
split_r     <- 30

prof <- read.table(map_csv, header = TRUE, sep = ",",
                   stringsAsFactors = FALSE, check.names = FALSE)
cb_col <- names(prof)[1]   # pandas writes the index (cb) as the first column

prom_cols <- grep("^prom[0-9]+$", names(prof), value = TRUE)
prom_cols <- prom_cols[order(as.integer(sub("^prom", "", prom_cols)))]
cum_cols  <- grep("^cumumi_r[0-9]+$", names(prof), value = TRUE)
cum_rad   <- as.integer(sub("^cumumi_r", "", cum_cols))
o <- order(cum_rad); cum_cols <- cum_cols[o]; cum_rad <- cum_rad[o]

cgrad <- function(x, yv) {
  n <- length(yv); g <- numeric(n)
  g[2:(n - 1)] <- (yv[3:n] - yv[1:(n - 2)]) / (x[3:n] - x[1:(n - 2)])
  g[1] <- (yv[2] - yv[1]) / (x[2] - x[1])
  g[n] <- (yv[n] - yv[n - 1]) / (x[n] - x[n - 1])
  g
}

# ---- background prominence metrics: single gamma, best-of{gamma,lognormal},
#      2-gamma mixture, and a parametric-free empirical reference ----
# typed one-row template (string best_family kept separate from numerics)
na_gamma <- function(n = NA_integer_) data.frame(
  n_peaks = n, gamma_resid_mad = NA_real_, p1_z = NA_real_, p2_z = NA_real_,
  best_family = NA_character_, best_resid_mad = NA_real_,
  best_p1_z = NA_real_, best_p2_z = NA_real_,
  mix_k = NA_integer_,
  mix_resid_mad = NA_real_, mix_p1_z = NA_real_, mix_p2_z = NA_real_,
  mix_vs_gamma_dbic = NA_real_,
  emp_bg_mad = NA_real_, prom1_bg_z = NA_real_, prom2_bg_z = NA_real_,
  stringsAsFactors = FALSE)

# Pareto type I (power law) MLE with xmin fixed at min(v). Returns alpha, xmin,
# AIC (1 free param), and a quantile function for order-statistic projection.
powerlaw_fit <- function(v) {
  xmin <- min(v); n <- length(v)
  s <- sum(log(v / xmin))
  if (!is.finite(s) || s <= 0) return(NULL)
  alpha <- n / s                              # shape of x^-(alpha+1) density
  loglik <- n * log(alpha) + n * alpha * log(xmin) - (alpha + 1) * sum(log(v))
  list(alpha = alpha, xmin = xmin, aic = -2 * loglik + 2,
       q = function(pp) xmin * (1 - pp)^(-1 / alpha))
}

# numeric order-statistic quantiles of a k-component mixture CDF, evaluated at pp.
.qmix <- function(pp, Fx_fun, lo, hi) {
  xs <- 10^seq(log10(lo), log10(hi), length.out = 6000)
  approx(Fx_fun(xs), xs, xout = pp, rule = 2)$y
}
qmix_gamma <- function(pp, lam, sh, sc, lo, hi)
  .qmix(pp, function(xs) Reduce(`+`, lapply(seq_along(lam), function(j)
        lam[j] * pgamma(xs, shape = sh[j], scale = sc[j]))), lo, hi)
qmix_lnorm <- function(pp, lam, mu, sg, lo, hi)
  .qmix(pp, function(xs) Reduce(`+`, lapply(seq_along(lam), function(j)
        lam[j] * plnorm(xs, meanlog = mu[j], sdlog = sg[j]))), lo, hi)

# Fit every candidate family to the background and return a named list of
# candidates, each = list(family, k, bic, exp=<order-stat expected vector at pp>).
# Families: single gamma / lognormal / powerlaw, gamma mixtures k=2..3, and
# lognormal mixtures k=2..3 (normalmixEM on natural-log data; raw-data loglik via
# Jacobian -sum(log bg) so BIC is comparable across all families).
build_candidates <- function(bgv, p, pp) {
  nb <- length(bgv); lo <- min(bgv) / 10; hi <- max(p) * 2
  cands <- list()
  addq <- function(name, family, k, npar, loglik, expv) {
    if (!is.finite(loglik) || all(!is.finite(expv))) return()
    cands[[name]] <<- list(family = family, k = k,
                           bic = -2 * loglik + log(nb) * npar, exp = expv)
  }
  g <- tryCatch(suppressWarnings(fitdistr(bgv, "gamma")), error = function(e) NULL)
  if (!is.null(g)) addq("gamma", "gamma", 1L, 2, g$loglik,
                        qgamma(pp, g$estimate[["shape"]], g$estimate[["rate"]]))
  l <- tryCatch(suppressWarnings(fitdistr(bgv, "lognormal")), error = function(e) NULL)
  if (!is.null(l)) addq("lognormal", "lognormal", 1L, 2, l$loglik,
                        qlnorm(pp, l$estimate[["meanlog"]], l$estimate[["sdlog"]]))
  pl <- tryCatch(powerlaw_fit(bgv), error = function(e) NULL)
  if (!is.null(pl)) {
    ll <- nb * log(pl$alpha) + nb * pl$alpha * log(pl$xmin) -
          (pl$alpha + 1) * sum(log(bgv))
    addq("powerlaw", "powerlaw", 1L, 1, ll, pl$q(pp))
  }
  if (have_mixtools && nb >= 30) {
    for (k in 2:3) {
      mm <- tryCatch(suppressWarnings({
              utils::capture.output(
                .mm <- mixtools::gammamixEM(bgv, k = k, maxit = 200, verb = FALSE))
              .mm}), error = function(e) NULL)
      if (!is.null(mm) && is.finite(mm$loglik))
        addq(sprintf("gamma_mix_k%d", k), "gamma_mix", k, 3 * k - 1, mm$loglik,
             qmix_gamma(pp, mm$lambda, mm$gamma.pars[1, ], mm$gamma.pars[2, ], lo, hi))
    }
    ly <- log(bgv); jac <- sum(ly)
    for (k in 2:3) {
      mm <- tryCatch(suppressWarnings({
              utils::capture.output(
                .mm <- mixtools::normalmixEM(ly, k = k, maxit = 400, verb = FALSE))
              .mm}), error = function(e) NULL)
      if (!is.null(mm) && is.finite(mm$loglik))
        addq(sprintf("lnorm_mix_k%d", k), "lnorm_mix", k, 3 * k - 1,
             mm$loglik - jac,                      # raw-data loglik (Jacobian)
             qmix_lnorm(pp, mm$lambda, mm$mu, mm$sigma, lo, hi))
    }
  }
  cands
}

# residual-MAD z-scores from a set of order-statistic expected prominences
resid_z <- function(p, exp_all, rk) {
  resid  <- log10(p) - log10(pmax(exp_all, 1e-300))
  bg_m   <- rk > 2 & p <= bg_prom_max & p >= prom_lo
  mad_bg <- if (sum(bg_m) >= 5) mad(resid[bg_m]) else NA_real_
  list(mad = mad_bg, z1 = resid[1] / mad_bg,
       z2 = if (length(resid) > 1) resid[2] / mad_bg else NA_real_)
}

gamma_metric <- function(prow) {
  p <- as.numeric(prow[prom_cols])
  p <- p[is.finite(p) & p > 0]
  out <- na_gamma(length(p))
  if (length(p) < 12) return(out)
  rk  <- seq_along(p)
  N   <- length(p)
  bgv <- p[rk > 2 & p <= bg_prom_max]
  if (length(bgv) < 10) return(out)
  pp  <- pmin(pmax(1 - (rk - 0.5) / N, 1e-7), 1 - 1e-9)

  # --- fit all candidate families once ---
  cands <- tryCatch(build_candidates(bgv, p, pp), error = function(e) list())
  if (!length(cands)) return(out)
  # attach residual-MAD z-scores to each candidate
  for (nm in names(cands)) {
    rz <- resid_z(p, cands[[nm]]$exp, rk)
    cands[[nm]]$mad <- rz$mad; cands[[nm]]$z1 <- rz$z1; cands[[nm]]$z2 <- rz$z2
  }
  bics <- sapply(cands, `[[`, "bic")
  mads <- sapply(cands, `[[`, "mad")

  # --- single gamma (original columns) ---
  gamma_bic <- NA_real_
  if (!is.null(cands$gamma)) {
    g <- cands$gamma
    out$gamma_resid_mad <- g$mad; out$p1_z <- g$z1; out$p2_z <- g$z2
    gamma_bic <- g$bic
  }

  # --- gamma-family mixture (single gamma + gamma_mix k2/3), BIC-best ---
  gfam <- cands[sapply(cands, function(c) c$family %in% c("gamma", "gamma_mix"))]
  if (length(gfam)) {
    ms <- gfam[[which.min(sapply(gfam, `[[`, "bic"))]]
    out$mix_k <- ms$k
    out$mix_resid_mad <- ms$mad; out$mix_p1_z <- ms$z1; out$mix_p2_z <- ms$z2
    if (is.finite(gamma_bic)) out$mix_vs_gamma_dbic <- gamma_bic - ms$bic
  }

  # --- overall best family by BIC, guarding against degenerate fits ---
  # (a heavy concentrated mixture component can win BIC yet give a pathological
  #  order-stat axis: huge resid MAD. drop any candidate whose MAD is >3x the
  #  median candidate MAD, or non-finite, before picking the lowest BIC.)
  ok <- is.finite(mads)
  if (any(ok)) {
    med_mad <- median(mads[ok])
    keep <- ok & mads <= 3 * med_mad
    if (!any(keep)) keep <- ok
    bsub <- bics[keep]
    bn <- names(bsub)[which.min(bsub)]
    b  <- cands[[bn]]
    out$best_family <- bn
    out$best_resid_mad <- b$mad; out$best_p1_z <- b$z1; out$best_p2_z <- b$z2
  }

  # --- (3) parametric-free empirical background ---
  bg_m  <- rk > 2 & p <= bg_prom_max & p >= prom_lo
  if (sum(bg_m) >= 5) {
    lbg  <- log10(p[bg_m])
    cen  <- median(lbg); sc <- mad(lbg)
    out$emp_bg_mad <- sc
    if (is.finite(sc) && sc > 0) {
      out$prom1_bg_z <- (log10(p[1]) - cen) / sc
      out$prom2_bg_z <- if (length(p) > 1) (log10(p[2]) - cen) / sc else NA_real_
    }
  }
  out
}

# ---- intden radial-density local-max ratio + z ----
NA_INTDEN <- c(intden_peak_ratio = NA, intden_peak_z = NA,
               intden_top_z = NA, intden_top_r = NA)
intden_metric <- function(prow) {
  y <- as.numeric(prow[cum_cols])
  if (all(!is.finite(y))) return(NA_INTDEN)
  d <- data.frame(radius = cum_rad, cumumi = y)
  d <- d[is.finite(d$radius) & is.finite(d$cumumi), , drop = FALSE]
  if (nrow(d) < 4) return(NA_INTDEN)
  bc <- (floor(d$radius / bin_width) + 0.5) * bin_width
  d  <- aggregate(cumumi ~ bc, data = d, FUN = mean); names(d) <- c("radius", "cumumi")
  d  <- d[order(d$radius), , drop = FALSE]
  if (nrow(d) < 4) return(NA_INTDEN)
  rg <- seq(0, max(d$radius), length.out = 400)
  # inner spline (reflected, weighted)
  base <- rbind(data.frame(radius = 0, cumumi = 0),
                d[d$radius > 0 & d$radius <= split_r, , drop = FALSE])
  refl <- data.frame(radius = -base$radius, cumumi = base$cumumi)
  refl <- refl[refl$radius < 0, , drop = FALSE]
  d_in <- rbind(refl, base); d_in <- d_in[order(d_in$radius), , drop = FALSE]
  if (length(unique(d_in$radius)) < 4) return(NA_INTDEN)
  w_in <- ifelse(d_in$radius == 0, 1e4, ifelse(abs(d_in$radius) < 25, 1e5, 1))
  ss_in <- smooth.spline(d_in$radius, d_in$cumumi, w = w_in, spar = smooth_spar)
  lev_in <- predict(ss_in, rg)$y
  outsub <- d[d$radius >= split_r, , drop = FALSE]
  if (length(unique(outsub$radius)) >= 4) {
    ss_out <- smooth.spline(outsub$radius, outsub$cumumi, spar = smooth_spar)
    off <- predict(ss_in, split_r)$y - predict(ss_out, split_r)$y
    lev_out <- predict(ss_out, rg)$y + off
    lev <- ifelse(rg <= split_r, lev_in, lev_out)
  } else lev <- lev_in
  lev <- cummax(lev)
  d1 <- cgrad(rg, lev)
  n1 <- length(d1); is_max <- rep(FALSE, n1)
  ii <- 2:(n1 - 1)
  is_max[ii] <- d1[ii] > d1[ii - 1] & d1[ii] > d1[ii + 1]
  lm_r <- rg[is_max]; lm_d <- d1[is_max]
  out <- NA_INTDEN
  # global highest local maximum vs all the other maxima (any radius)
  if (length(lm_d) >= 3) {
    gi  <- which.max(lm_d)
    others <- lm_d[-gi]
    mg <- mad(others)
    out["intden_top_r"] <- lm_r[gi]
    out["intden_top_z"] <- if (is.finite(mg) && mg > 0)
                             (lm_d[gi] - median(others)) / mg else NA_real_
  }
  # within-50um peak vs outer maxima (original metric)
  in50 <- lm_d[lm_r <= 50]; rest <- lm_d[lm_r > 50]
  if (length(in50) >= 1 && length(rest) >= 2) {
    peak <- max(in50); med <- median(rest); m <- mad(rest)
    out["intden_peak_ratio"] <- peak / med
    out["intden_peak_z"] <- if (m > 0) (peak - med) / m else NA_real_
  }
  out
}

n <- nrow(prof)
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "0"))
if (is.na(ncores) || ncores < 1) ncores <- max(1L, parallel::detectCores() - 1L)

one_row <- function(i) {
  prow <- prof[i, ]
  gm <- tryCatch(gamma_metric(prow), error = function(e) na_gamma(NA_integer_))
  im <- tryCatch(intden_metric(prow), error = function(e) NA_INTDEN)
  cbind(gm, data.frame(
    intden_peak_ratio = unname(im["intden_peak_ratio"]),
    intden_peak_z     = unname(im["intden_peak_z"]),
    intden_top_z      = unname(im["intden_top_z"]),
    intden_top_r      = unname(im["intden_top_r"])
  ))
}

cat("computing metrics for", n, "cells on", ncores, "cores\n")
res <- parallel::mclapply(seq_len(n), one_row, mc.cores = ncores)
metrics <- do.call(rbind, res)

# Merge in place: drop any pre-existing metric columns, then cbind fresh ones.
prof <- prof[, !names(prof) %in% names(metrics), drop = FALSE]
prof <- cbind(prof, metrics)

# Preserve the leading cb column header (pandas index style) on write.
write.table(prof, out_csv, sep = ",", quote = FALSE, row.names = FALSE)
ok <- sum(!is.na(metrics$gamma_resid_mad))
cat("merged metrics into", out_csv, " rows:", nrow(prof),
    " gamma non-NA:", ok,
    " intden non-NA:", sum(!is.na(metrics$intden_peak_ratio)), "\n")
