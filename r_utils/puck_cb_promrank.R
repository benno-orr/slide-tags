# log10(prominence) vs peak rank, with the pipeline's leftmost-line fit(s)
# overlaid. Slope/intercept are not stored, so each line is re-fit from the prom
# values within its stored rank window ({linrank,logrank}_lr_sec_lo/_sec_hi),
# reproducing map_cells.py. peakK_*_lr_fc = vertical gap from the line up to peak
# K (dashed = linear-rank, dotted = log-rank drop-lines). When show_regions, two
# extra panels reproduce the region-counting metric ({linrank,logrank}_n_regions):
# the smoothed-log-prom derivative + its spline, with detected transition minima.
# Columns expected (cell_coords.csv, map_cells.py): peakK_prom, peakK_height,
# {linrank,logrank}_lr_sec_lo/_sec_hi, peak{1,2}_bestfit_lr_fc, {linrank,logrank}_n_regions.
# Depends: .get_cell_row(), .prom_region_fit().
puck_cb_promrank <- function(cell_id, profiles_df, suffix = "-1",
                           fit_prom_max = 0.1,
                           max_rank = NULL, show_fit = TRUE,
                           fit_bg = FALSE, bg_prom_max = 0.1,
                           fit_gamma = FALSE,
                           gamma_model = c("mixture", "single", "powerlaw"),
                           spline_spar = 0.45,
                           show_regions = TRUE, region_min_prom = 0.05,
                           region_top_exclude = 3,
                           show_linregions = TRUE, linreg_win = 11,
                           linreg_r2 = 0.97, linreg_min_len = 8) {
  gamma_model <- match.arg(gamma_model)
  row <- .get_cell_row(profiles_df, cell_id, suffix)
  if (is.null(row)) return(ggplot() + ggtitle(paste0(cell_id, " — not in table")))

  prom_cols <- grep("^peak[0-9]+_prom$", names(row), value = TRUE)
  prom_cols <- prom_cols[order(as.integer(gsub("[^0-9]", "", prom_cols)))]
  p <- as.numeric(row[1, prom_cols])
  keep <- is.finite(p) & p > 0
  p <- p[keep]
  if (length(p) == 0) return(ggplot() + ggtitle(paste0(cell_id, " — no peaks")))

  d <- data.frame(rank = seq_along(p), log_prom = log10(p))

  # gamma background model: fit a gamma to the background prominences
  # (rank > 2, prom <= bg_prom_max), project the expected prominence-by-rank via
  # the order-statistic plotting position qgamma(1 - (rank-0.5)/N), and build a
  # linearised view (observed vs gamma-expected prominence) where the background
  # falls on the 1:1 line and real peaks sit above it.
  gproj <- glin <- gpar <- exp_all <- NULL
  if (fit_gamma) {
    bgv <- p[seq_along(p) > 2 & p <= bg_prom_max]
    N  <- length(p)
    pp <- pmin(pmax(1 - (d$rank - 0.5) / N, 1e-7), 1 - 1e-9)
    if (gamma_model == "single" && length(bgv) >= 10) {
      gf <- tryCatch(suppressWarnings(MASS::fitdistr(bgv, "gamma")),
                     error = function(e) NULL)
      if (!is.null(gf)) {
        sh <- gf$estimate[["shape"]]; rt <- gf$estimate[["rate"]]
        exp_all <- qgamma(pp, sh, rt)
        gpar <- sprintf("gamma: k=%.2f rate=%.0f", sh, rt)
      }
    } else if (gamma_model == "mixture" && length(bgv) >= 30) {
      # gamma mixture, #components auto-selected by BIC over k = 1..3
      n_bg <- length(bgv)
      fit_k <- function(k) {
        if (k == 1) {
          gf <- tryCatch(suppressWarnings(MASS::fitdistr(bgv, "gamma")),
                         error = function(e) NULL)
          if (is.null(gf)) return(NULL)
          list(k = 1, lam = 1, shapes = gf$estimate[["shape"]],
               scales = 1 / gf$estimate[["rate"]],
               bic = -2 * gf$loglik + log(n_bg) * 2)
        } else {
          mm <- tryCatch(suppressWarnings(
                  mixtools::gammamixEM(bgv, k = k, maxit = 400, verb = FALSE)),
                  error = function(e) NULL)
          if (is.null(mm) || is.null(mm$loglik) || !is.finite(mm$loglik)) return(NULL)
          list(k = k, lam = mm$lambda, shapes = mm$gamma.pars[1, ],
               scales = mm$gamma.pars[2, ],
               bic = -2 * mm$loglik + log(n_bg) * (3 * k - 1))
        }
      }
      cands <- Filter(Negate(is.null), lapply(1:3, fit_k))
      if (length(cands)) {
        ms <- cands[[which.min(sapply(cands, `[[`, "bic"))]]
        xs <- 10^seq(log10(min(bgv) / 10), log10(max(p) * 2), length.out = 4000)
        Fx <- Reduce(`+`, lapply(seq_along(ms$lam), function(j)
                ms$lam[j] * pgamma(xs, shape = ms$shapes[j], scale = ms$scales[j])))
        exp_all <- approx(Fx, xs, xout = pp, rule = 2)$y
        gpar <- sprintf("gamma mix (BIC k=%d): w=%s k=%s", ms$k,
                        paste(sprintf("%.2f", ms$lam), collapse = "/"),
                        paste(sprintf("%.1f", ms$shapes), collapse = "/"))
      }
    } else if (gamma_model == "powerlaw" && length(bgv) >= 10) {
      # Pareto type I (xmin = min(bg); MLE alpha) — a straight line in log-log
      # rank-prominence. order-stat quantile = xmin * (1 - pp)^(-1/alpha).
      xmin <- min(bgv); s <- sum(log(bgv / xmin))
      if (is.finite(s) && s > 0) {
        alpha <- length(bgv) / s
        exp_all <- xmin * (1 - pp)^(-1 / alpha)
        gpar <- sprintf("power law: alpha=%.2f xmin=%.1e", alpha, xmin)
      }
    }
    if (!is.null(exp_all)) {
      gproj <- data.frame(rank = d$rank, log_prom = log10(pmax(exp_all, 1e-300)))
      glin  <- data.frame(rank = d$rank, obs = p, exp = exp_all)
    }
  }

  # peak-vs-gamma residual z (source of the gamma_z metric): residual of each
  # point from the gamma line in log space, scaled by the MAD of the background
  # residuals (rank > 2, 1e-5 <= prom <= bg_prom_max). computed here so the
  # main prom-rank plot can annotate peak1_gamma_z.
  resid <- bg_m <- NULL; mad_bg <- NA_real_
  if (!is.null(glin)) {
    resid  <- log10(glin$obs) - log10(glin$exp)
    bg_m   <- glin$rank > 2 & glin$obs <= bg_prom_max & glin$obs >= 1e-5
    mad_bg <- if (sum(bg_m) >= 5) mad(resid[bg_m]) else NA_real_
    glin$z <- resid / mad_bg
  }

  num <- function(nm) suppressWarnings(as.numeric(row[[nm]][1]))
  kk  <- seq_len(min(2, length(p)))
  lr_lfc_lin <- lr_lfc_log <- rep(NA_real_, 2)

  # leftmost-line ("lr") fits, re-fit from each line's STORED rank window so the
  # overlay reproduces map_cells.py. The linear-rank fit (log10(prom) ~ rank,
  # orange) is drawn ONLY on the linear-rank panel; the log-rank fit
  # (log10(prom) ~ log10(rank), blue) ONLY on the log-rank panel. The points used
  # in each regression (rank within [sec_lo, sec_hi], prom <= fit_prom_max) are
  # ringed, and the panel reports the fit R². No green smoothing curve.
  lin_fit <- NULL
  lin_lo <- num("linrank_lr_sec_lo"); lin_hi <- num("linrank_lr_sec_hi")
  if (length(lin_lo) == 1 && is.finite(lin_lo) && is.finite(lin_hi)) {
    w <- d$rank >= lin_lo & d$rank <= lin_hi & p <= fit_prom_max
    if (sum(w) >= 2) {
      m  <- stats::lm(log_prom ~ rank, data = d[w, , drop = FALSE])
      ca <- stats::coef(m)
      la <- data.frame(rank = seq(1, lin_hi, length.out = 200))
      la$log_prom <- ca[[1]] + ca[[2]] * la$rank
      lr_lfc_lin[kk] <- log10(p[kk]) - (ca[[1]] + ca[[2]] * kk)
      lin_fit <- list(
        line = la,
        used = d[w, , drop = FALSE],
        drop = data.frame(rank = kk, y = log10(p[kk]), yend = ca[[1]] + ca[[2]] * kk),
        r2 = summary(m)$r.squared, npts = sum(w))
    }
  }

  log_fit <- NULL
  log_lo <- num("logrank_lr_sec_lo"); log_hi <- num("logrank_lr_sec_hi")
  if (length(log_lo) == 1 && is.finite(log_lo) && is.finite(log_hi)) {
    w <- d$rank >= log_lo & d$rank <= log_hi & p <= fit_prom_max
    if (sum(w) >= 2) {
      m  <- stats::lm(log_prom ~ log10(rank), data = d[w, , drop = FALSE])
      cb <- stats::coef(m)
      lb <- data.frame(rank = seq(1, log_hi, length.out = 200))
      lb$log_prom <- cb[[1]] + cb[[2]] * log10(lb$rank)
      lr_lfc_log[kk] <- log10(p[kk]) - (cb[[1]] + cb[[2]] * log10(kk))
      log_fit <- list(
        line = lb,
        used = d[w, , drop = FALSE],
        drop = data.frame(rank = kk, y = log10(p[kk]), yend = cb[[1]] + cb[[2]] * log10(kk)),
        r2 = summary(m)$r.squared, npts = sum(w))
    }
  }

  # peaks 1 & 2 highlight (shared)
  pk <- data.frame(rank = c(1, 2), log_prom = log10(p[1:2]))[seq_len(min(2, length(p))), ]

  # y range from the data points; restrict to ranks <= max_rank when set
  dv <- if (is.null(max_rank)) d else d[d$rank <= max_rank, , drop = FALSE]
  yr <- range(dv$log_prom, na.rm = TRUE)
  xr <- c(1, if (is.null(max_rank)) max(d$rank) else max_rank)
  bf1 <- num("peak1_bestfit_lr_fc"); bf2 <- num("peak2_bestfit_lr_fc")
  n_reg_lin <- num("linrank_n_regions"); n_reg_log <- num("logrank_n_regions")

  # roughly-linear regions in each rank space (sliding-window local R²), shaded
  # green on the matching panel: linrank regions on the linear panels, logrank
  # regions on the log panel.
  lr_lin_reg <- lr_log_reg <- NULL
  if (show_linregions) {
    lr_lin_reg <- .linear_regions(p, "linrank", linreg_win, linreg_r2, linreg_min_len)
    lr_log_reg <- .linear_regions(p, "logrank", linreg_win, linreg_r2, linreg_min_len)
  }
  shade <- function(g, reg) {
    if (is.null(reg)) return(g)
    g + geom_rect(data = reg, inherit.aes = FALSE,
                  aes(xmin = rank_lo, xmax = rank_hi, ymin = -Inf, ymax = Inf),
                  fill = "forestgreen", alpha = 0.13)
  }
  nreg_txt <- function(reg) if (is.null(reg)) "0" else as.character(nrow(reg))

  # base scatter for a panel: shaded linear regions (behind) -> gamma proj -> pts
  mk_base <- function(reg) {
    g <- shade(ggplot(d, aes(rank, log_prom)), reg)
    if (!is.null(gproj)) g <- g + geom_line(data = gproj, colour = "purple", linewidth = 0.8)
    g + geom_point(size = 1.3, colour = "grey40")
  }

  # add one lr fit (+ its used-points ring + drop-lines) to a base plot
  add_fit <- function(g, fit, col, lty) {
    if (is.null(fit)) return(g)
    g +
      geom_point(data = fit$used, aes(rank, log_prom), shape = 21,
                 colour = col, fill = NA, size = 3, stroke = 0.7) +
      geom_line(data = fit$line, colour = col, linewidth = 0.8) +
      geom_segment(data = fit$drop, aes(x = rank, xend = rank, y = y, yend = yend),
                   linetype = lty, colour = col, inherit.aes = FALSE)
  }
  r2txt <- function(fit) if (is.null(fit)) "NA"
                         else sprintf("R²=%.3f, n=%d pts ringed", fit$r2, fit$npts)

  # linear-rank panel as two side-by-side views: full (left) + zoom to peaks
  # 1..zoom_hi (right). The zoom uses its own y-range over the shown ranks.
  zoom_hi <- min(20, max(d$rank))
  yr_zoom <- range(d$log_prom[d$rank <= zoom_hi], na.rm = TRUE)
  lin_base <- function() add_fit(mk_base(lr_lin_reg), lin_fit, "darkorange2", "dashed") +
    geom_point(data = pk, aes(rank, log_prom), colour = "firebrick", size = 2.4) +
    labs(y = "log10(prom) [UMI/µm²]", title = NULL) + theme_minimal()
  g_lin_full <- lin_base() +
    coord_cartesian(xlim = xr, ylim = yr) +
    labs(x = "peak rank (linear)",
         subtitle = sprintf("linear-rank lr fit  %s | fc p1=%.2f p2=%.2f | linregions(green)=%s",
                            r2txt(lin_fit), lr_lfc_lin[1], lr_lfc_lin[2], nreg_txt(lr_lin_reg)))
  g_lin_zoom <- lin_base() +
    coord_cartesian(xlim = c(1, zoom_hi), ylim = yr_zoom) +
    labs(x = "peak rank (linear)", y = NULL,
         subtitle = sprintf("zoom: peaks 1-%d", zoom_hi))
  g_lin_rank <- g_lin_full | g_lin_zoom

  g_log_rank <- add_fit(mk_base(lr_log_reg), log_fit, "steelblue", "dotted") +
    geom_point(data = pk, aes(rank, log_prom), colour = "firebrick", size = 2.4) +
    scale_x_log10() +
    coord_cartesian(xlim = xr, ylim = yr) +
    labs(x = "peak rank (log)", y = "log10(prom) [UMI/µm²]", title = NULL,
         subtitle = sprintf("log-rank lr fit  %s | fc p1=%.2f p2=%.2f | linregions(green)=%s",
                            r2txt(log_fit), lr_lfc_log[1], lr_lfc_log[2], nreg_txt(lr_log_reg))) +
    theme_minimal()
  main_panels <- list(g_lin_rank, g_log_rank)

  # ---- region metric panels (linear-rank only): 1st derivative of the smoothed
  # ---- log-prom curve (with transition minima behind linrank_n_regions), plus
  # ---- the 2nd derivative. Linear rank x-axis (the derivative is wrt rank).
  reg_panels <- list()
  if (show_regions) {
    rf <- .prom_region_fit(p, "linrank", region_min_prom, region_top_exclude)
    if (!is.null(rf)) {
      df <- rf$df; mins <- df[df$is_min, , drop = FALSE]
      p_d1 <- ggplot(df, aes(rank, spline)) +
        geom_hline(yintercept = 0, colour = "grey85") +
        geom_line(aes(y = deriv), colour = "grey70", linewidth = 0.4) +
        geom_line(colour = "forestgreen", linewidth = 0.7)
      if (nrow(mins) > 0)
        p_d1 <- p_d1 + geom_point(data = mins, colour = "firebrick", size = 2.2) +
          geom_vline(data = mins, aes(xintercept = rank),
                     linetype = "dotted", colour = "firebrick")
      p_d1 <- p_d1 + coord_cartesian(xlim = xr) +
        labs(x = NULL, y = "d(smooth logP)/dx",
             subtitle = sprintf("linrank 1st deriv: transitions(refit)=%d  stored n_regions=%s",
                                rf$n_transitions,
                                ifelse(is.finite(n_reg_lin), as.character(round(n_reg_lin)), "NA"))) +
        theme_minimal()

      p_d2 <- ggplot(df, aes(rank, spline2)) +
        geom_hline(yintercept = 0, colour = "grey85") +
        geom_line(colour = "purple", linewidth = 0.7)
      if (nrow(mins) > 0)
        p_d2 <- p_d2 + geom_vline(data = mins, aes(xintercept = rank),
                                  linetype = "dotted", colour = "firebrick")
      p_d2 <- p_d2 + coord_cartesian(xlim = xr) +
        labs(x = "peak rank (linear)", y = "d²(smooth logP)/dx²",
             subtitle = "linrank 2nd derivative") +
        theme_minimal()

      reg_panels <- list(p_d1, p_d2)
    }
  }

  # no gamma fit: show both prom-rank panels (+ region panels) with the lr fits
  if (is.null(glin)) {
    panels  <- c(main_panels, reg_panels)
    heights <- c(2, 2, rep(1, length(reg_panels)))
    return(patchwork::wrap_plots(panels, ncol = 1, heights = heights))
  }

  # deviation z-score (computed above): residual from the gamma line in log
  # space, scaled by the MAD of the background residuals. tight gamma fit
  # (small bg MAD) => larger z for the same fold; loose fit => smaller z.
  r12 <- seq_len(min(2, nrow(glin)))
  zsub <- sprintf("bg resid mad=%.3f   p1 z=%.1f   p2 z=%.1f",
                  mad_bg, glin$z[1], if (length(r12) > 1) glin$z[2] else NA)
  # which model/combo was fitted (e.g. "gamma mix (BIC k=3): w=.. k=..") on top line
  g_title <- if (!is.null(gpar)) gpar else NULL

  pk_lin <- glin[r12, , drop = FALSE]
  g_gam <- ggplot(glin, aes(exp, obs)) +
    geom_abline(slope = 1, intercept = 0, colour = "purple", linewidth = 0.8) +
    geom_point(size = 1, colour = "grey40") +
    geom_point(data = pk_lin, colour = "firebrick", size = 2.4) +
    geom_text(data = pk_lin, aes(label = sprintf("z=%.1f", z)),
              hjust = -0.2, size = 3, colour = "firebrick") +
    scale_x_log10() + scale_y_log10() +
    labs(x = "model-expected prominence", y = "observed prominence",
         title = g_title, subtitle = zsub) +
    theme_minimal() +
    theme(plot.title = element_text(size = 9, face = "plain", colour = "grey20"))

  # distribution of the background residuals that MAD summarises, with the ±MAD
  # band (blue dashed), the gamma line at 0 (purple) and peaks 1/2 (firebrick).
  rd <- data.frame(resid = resid[bg_m])
  p_resid <- ggplot(rd, aes(resid)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40,
                   fill = "grey85", colour = NA) +
    geom_density(colour = "grey30", linewidth = 0.8, adjust = 1.2) +
    geom_vline(xintercept = 0, colour = "purple", linewidth = 0.7) +
    geom_vline(xintercept = c(-mad_bg, mad_bg), colour = "steelblue",
               linetype = "dashed") +
    geom_vline(xintercept = resid[r12], colour = "firebrick", linewidth = 0.7) +
    labs(x = "log10(obs) - log10(gamma exp)", y = "density",
         subtitle = "bg residuals; dashed = ±MAD, red = peaks") +
    theme_minimal()

  panels  <- c(main_panels, reg_panels, list(g_gam, p_resid))
  heights <- c(2, 2, rep(1, length(reg_panels)), 2, 1)
  patchwork::wrap_plots(panels, ncol = 1, heights = heights)
}
