# cumulative UMI within radius r of the c1 centroid (cumumi_r001..r400 from the
# profiles row). No curve fitting: the lower panels show the RAW first derivative
# d(cumUMI)/dr and the difference in fraction-of-max Δ(cumUMI / max(cumUMI)).
# Each is a `step`-lagged slope — between point n and n+step (default 2; use 3 for
# more smoothing) — placed at the midpoint radius, rather than an adjacent
# (lag-1) difference. normalize = TRUE divides the cumulative curve by the value
# at `norm_r` (default the largest radius) for cross-cell comparison; the
# fraction-of-max difference always uses each cell's own max.
# Depends: .get_cell_row().
puck_c1_intden <- function(cell_id, profiles_df, suffix = "-1",
                               normalize = FALSE, norm_r = NULL,
                               step = 2, roll = 1, title = NULL) {
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

  d <- data.frame(radius = radius, cumumi = y)
  d <- d[is.finite(d$radius) & is.finite(d$cumumi), , drop = FALSE]
  d <- d[order(d$radius), , drop = FALSE]
  ymax <- max(d$cumumi[is.finite(d$cumumi)])

  ylab <- "cumulative UMI"
  cumdisp <- d$cumumi
  if (normalize) {
    ref_r <- if (is.null(norm_r)) max(d$radius) else norm_r
    denom <- d$cumumi[which.min(abs(d$radius - ref_r))]
    if (is.finite(denom) && denom > 0) { cumdisp <- d$cumumi / denom; ylab <- "fraction of UMI" }
  }

  # step-lagged slope between point i and i+step (placed at the midpoint radius):
  #   d1   = (cumUMI[i+step]-cumUMI[i]) / (r[i+step]-r[i])  -> 1st derivative (UMI/µm)
  #   dfom = frac[i+step]-frac[i]                           -> diff in fraction of max
  frac <- d$cumumi / ymax
  n    <- nrow(d)
  step <- max(1L, as.integer(step))
  adj  <- if (n > step) {
            i <- seq_len(n - step); j <- i + step
            data.frame(
              radius = (d$radius[i] + d$radius[j]) / 2,
              d1     = (d$cumumi[j] - d$cumumi[i]) / (d$radius[j] - d$radius[i]),
              dfom   = frac[j] - frac[i])
          } else NULL

  # optional light centred rolling mean (roll > 1) for legibility; raw by default
  if (!is.null(adj) && roll > 1) {
    k <- as.integer(roll)
    rmean <- function(v) as.numeric(stats::filter(v, rep(1 / k, k), sides = 2))
    adj$d1   <- ifelse(is.na(rmean(adj$d1)),   adj$d1,   rmean(adj$d1))
    adj$dfom <- ifelse(is.na(rmean(adj$dfom)), adj$dfom, rmean(adj$dfom))
  }

  # top: cumulative UMI vs radius (raw, no fit), with fraction-of-max on the right
  p <- ggplot(data.frame(radius = d$radius, cumumi = cumdisp),
              aes(radius, cumumi)) +
    geom_line(colour = "grey55", linewidth = 0.5) +
    geom_point(size = 0.6, colour = "grey30") +
    scale_y_continuous(
      sec.axis = sec_axis(~ . / max(cumdisp), name = "fraction of max",
                          breaks = seq(0, 1, 0.25))
    ) +
    labs(x = NULL, y = ylab, title = ttl) +
    theme_minimal()

  if (is.null(adj)) return(p)

  # radius of the strongest raw derivative (core radius) for annotation
  rmax <- adj$radius[which.max(adj$d1)]

  p_d1 <- ggplot(adj, aes(radius, d1)) +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_line(colour = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = rmax, linetype = "dotted", colour = "firebrick") +
    labs(x = NULL, y = "d(cumUMI)/dr  [UMI/µm]",
         subtitle = sprintf("1st derivative (slope over %d-pt step); max at r=%g µm",
                            step, rmax)) +
    theme_minimal()

  p_df <- ggplot(adj, aes(radius, dfom)) +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_line(colour = "steelblue", linewidth = 0.5) +
    geom_vline(xintercept = rmax, linetype = "dotted", colour = "firebrick") +
    labs(x = "radius from c1 centroid (µm)", y = "Δ fraction of max",
         subtitle = sprintf("difference in fraction of max (over %d-pt step)", step)) +
    theme_minimal()

  patchwork::wrap_plots(p, p_d1, p_df, ncol = 1, heights = c(2, 1, 1))
}
