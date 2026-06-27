# region counting on the log10(prominence)-vs-rank curve — R reimplementation of
# map_cells.py `_count_regions`, so the {linrank,logrank}_n_regions / n_transitions
# metrics can be inspected per cell. Steps: smooth log-prom with a [.25,.5,.25]
# kernel, take d/dx in the chosen rank space (linrank: x=rank; logrank: x=log10
# rank), spline-smooth the derivative, then count local minima whose dip-prominence
# (drop from the higher flanking ridge) >= min_prom. transitions = those minima;
# regions = transitions + 1. Returns the per-rank fit frame + the two counts.
.prom_region_fit <- function(prom, space = c("linrank", "logrank"),
                             min_prom = 0.05) {
  space <- match.arg(space)
  lp <- log10(prom)
  n  <- length(lp)
  if (n < 5) return(NULL)

  # [0.25, 0.5, 0.25] smoothing with edge padding (preserves length)
  padded <- c(lp[1], lp, lp[n])
  sm <- as.numeric(stats::filter(padded, c(0.25, 0.5, 0.25)))[2:(n + 1)]

  rank <- seq_len(n)
  x <- if (space == "linrank") rank else log10(rank)

  # central-difference derivative wrt x
  g <- numeric(n)
  g[2:(n - 1)] <- (sm[3:n] - sm[1:(n - 2)]) / (x[3:n] - x[1:(n - 2)])
  g[1] <- (sm[2] - sm[1]) / (x[2] - x[1])
  g[n] <- (sm[n] - sm[n - 1]) / (x[n] - x[n - 1])

  # spline-smooth the derivative (scipy UnivariateSpline analogue)
  ss <- tryCatch(stats::smooth.spline(x, g), error = function(e) NULL)
  sp <- if (is.null(ss)) g else stats::predict(ss, x)$y

  # local minima with dip-prominence >= min_prom
  is_min <- rep(FALSE, n)
  if (n >= 3) for (i in 2:(n - 1)) {
    if (sp[i] < sp[i - 1] && sp[i] <= sp[i + 1]) {
      dip <- min(max(sp[1:(i - 1)]), max(sp[i:n])) - sp[i]
      if (is.finite(dip) && dip >= min_prom) is_min[i] <- TRUE
    }
  }

  list(df = data.frame(rank = rank, x = x, smooth = sm, deriv = g,
                       spline = sp, is_min = is_min),
       n_transitions = sum(is_min),
       n_regions     = sum(is_min) + 1L)
}
