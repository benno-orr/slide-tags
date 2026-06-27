# cumulative UMI within radius r of the c1 centroid (cumumi_r001..r400 from the
# profiles row). normalize = TRUE divides by the value at `norm_r` (default the
# largest radius) so curves are comparable across cells.
# Depends: .get_cell_row().
puck_c1_intden <- function(cell_id, profiles_df, suffix = "-1",
                               normalize = FALSE, norm_r = NULL,
                               bin_width = 3, smooth_spar = 0.55,
                               title = NULL) {
  ttl <- if (missing(title)) cell_id else title
  row <- .get_cell_row(profiles_df, cell_id, suffix)
  if (is.null(row))
    return(ggplot() + ggtitle(paste0(cell_id, " — not in table")) + theme_void())

  cum_cols <- grep("^cumumi_r[0-9]+$", names(row), value = TRUE)
  radius   <- as.integer(sub("^cumumi_r", "", cum_cols))
  ord      <- order(radius)
  radius   <- radius[ord]; cum_cols <- cum_cols[ord]
  y <- as.numeric(row[1, cum_cols])
  if (all(!is.finite(y)))
    return(ggplot() + theme_void())

  ylab <- "cumulative UMI"
  if (normalize) {
    ref_r <- if (is.null(norm_r)) max(radius) else norm_r
    denom <- y[which.min(abs(radius - ref_r))]
    if (is.finite(denom) && denom > 0) { y <- y / denom; ylab <- "fraction of UMI" }
  }

  d <- data.frame(radius = radius, cumumi = y)
  d <- d[is.finite(d$radius) & is.finite(d$cumumi), , drop = FALSE]

  # smooth the points into radius bins of width bin_width (mean cumumi per bin,
  # bin centre as radius) before fitting.
  if (!is.null(bin_width) && bin_width > 0 && nrow(d) > 0) {
    bin_ctr <- (floor(d$radius / bin_width) + 0.5) * bin_width
    d <- aggregate(cumumi ~ bin_ctr, data = d, FUN = mean)
    names(d) <- c("radius", "cumumi")
    d <- d[order(d$radius), , drop = FALSE]
  }

  ymax <- max(d$cumumi[is.finite(d$cumumi)])

  # two independent smoothing splines joined at 30 µm so the heavy inner
  # weighting doesn't starve the outer fit:
  #   inner [0, 30]: origin-anchored (flat start), tightly weighted.
  #   outer [30, max]: its own spline, shifted to meet the inner one at 30.
  # forced always-increasing with cummax.
  cgrad <- function(x, yv) {
    n <- length(yv); g <- numeric(n)
    g[2:(n - 1)] <- (yv[3:n] - yv[1:(n - 2)]) / (x[3:n] - x[1:(n - 2)])
    g[1] <- (yv[2] - yv[1]) / (x[2] - x[1])
    g[n] <- (yv[n] - yv[n - 1]) / (x[n] - x[n - 1])
    g
  }
  split_r <- 30

  fit_df <- deriv_df <- NULL
  if (nrow(d) >= 4) {
    rg <- seq(0, max(d$radius), length.out = 400)

    # inner spline: origin-anchored via reflection about 0, tight
    base <- rbind(data.frame(radius = 0, cumumi = 0),
                  d[d$radius > 0 & d$radius <= split_r, , drop = FALSE])
    refl <- data.frame(radius = -base$radius, cumumi = base$cumumi)
    refl <- refl[refl$radius < 0, , drop = FALSE]
    d_in <- rbind(refl, base)
    d_in <- d_in[order(d_in$radius), , drop = FALSE]
    w_in <- ifelse(d_in$radius == 0, 1e4, 1e5)
    ss_in <- stats::smooth.spline(d_in$radius, d_in$cumumi, w = w_in,
                                  spar = smooth_spar)
    lev_in <- predict(ss_in, rg)$y

    # outer spline: separate fit on points past the split (include the split
    # point for continuity), shifted to meet the inner spline at split_r
    outsub <- d[d$radius >= split_r, , drop = FALSE]
    if (length(unique(outsub$radius)) >= 4) {
      ss_out  <- stats::smooth.spline(outsub$radius, outsub$cumumi, spar = smooth_spar)
      off     <- predict(ss_in, split_r)$y - predict(ss_out, split_r)$y
      lev_out <- predict(ss_out, rg)$y + off
      lev <- ifelse(rg <= split_r, lev_in, lev_out)
    } else {
      lev <- lev_in
    }
    lev <- cummax(lev)

    fit_df   <- data.frame(radius = rg, cumumi = lev)
    deriv_df <- data.frame(radius = rg, d1 = cgrad(rg, lev),
                           d2 = cgrad(rg, cgrad(rg, lev)))
  }

  p <- ggplot(d, aes(radius, cumumi)) +
    geom_point(size = 0.7, colour = "grey30")
  if (!is.null(fit_df))
    p <- p + geom_line(data = fit_df, colour = "grey75", linewidth = 0.8)
  else
    p <- p + geom_line(colour = "grey75", linewidth = 0.7)

  p <- p +
    scale_y_continuous(
      sec.axis = sec_axis(~ . / ymax, name = "fraction of max",
                          breaks = seq(0, 1, 0.25))
    ) +
    labs(x = NULL, y = ylab, title = ttl) +
    theme_minimal()

  if (is.null(deriv_df)) return(p)

  # local maxima of the radial density (1st deriv). Score the strongest maximum
  # within the first 50 µm against the rest: ratio = peak/median(rest), and a
  # robust z = (peak - median(rest)) / MAD(rest).
  d1v <- deriv_df$d1; n1 <- length(d1v)
  is_max <- rep(FALSE, n1)
  if (n1 >= 3) {
    ii <- 2:(n1 - 1)
    is_max[ii] <- d1v[ii] > d1v[ii - 1] & d1v[ii] > d1v[ii + 1]
  }
  lm_df <- deriv_df[is_max, , drop = FALSE]
  in50  <- lm_df[lm_df$radius <= 50, , drop = FALSE]
  rest  <- lm_df[lm_df$radius >  50, , drop = FALSE]
  msub <- NULL; peak_pt <- NULL
  if (nrow(in50) > 0 && nrow(rest) >= 2) {
    peak_pt <- in50[which.max(in50$d1), , drop = FALSE]
    med <- median(rest$d1); mad <- mad(rest$d1)
    ratio <- peak_pt$d1 / med
    z <- if (mad > 0) (peak_pt$d1 - med) / mad else NA_real_
    msub <- sprintf("peak<=50µm / median(rest) = %.2f   z = %.2f", ratio, z)
  }

  # 1st derivative (radial UMI density) of the fit, with local maxima marked
  p_d1 <- ggplot(deriv_df, aes(radius, d1)) +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_line(colour = "grey40", linewidth = 0.7) +
    geom_point(data = lm_df, colour = "grey55", size = 1.2)
  if (!is.null(peak_pt))
    p_d1 <- p_d1 + geom_point(data = peak_pt, colour = "firebrick", size = 2.4)
  p_d1 <- p_d1 +
    labs(x = "radius from c1 centroid (µm)", y = "1st deriv", subtitle = msub) +
    theme_minimal()

  patchwork::wrap_plots(p, p_d1, ncol = 1, heights = c(2, 1))
}
