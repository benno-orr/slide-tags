# Identify contiguous rank ranges over which log10(prom) is ~linear in the chosen
# rank space (linrank: x = rank; logrank: x = log10 rank). Slides a window of
# `win` points, fits a line, and flags ranks whose local fit R² >= r2_thr; runs
# of >= min_len consecutive flagged ranks become "linear regions". Returns a
# data.frame of regions (rank_lo, rank_hi, npts, r2 = mean local R² in the run),
# or NULL if none. Used to shade the prom-rank panels.
.linear_regions <- function(prom, space = c("linrank", "logrank"),
                            win = 11, r2_thr = 0.97, min_len = 8) {
  space <- match.arg(space)
  p <- prom[is.finite(prom) & prom > 0]
  n <- length(p)
  if (n < max(win, min_len)) return(NULL)
  y <- log10(p)
  rank <- seq_len(n)
  x <- if (space == "linrank") rank else log10(rank)

  h <- win %/% 2
  loc_r2 <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    lo <- max(1, i - h); hi <- min(n, i + h)
    if (hi - lo + 1 < 4) next
    xs <- x[lo:hi]; ys <- y[lo:hi]
    if (length(unique(xs)) < 2) next
    ss_tot <- sum((ys - mean(ys))^2)
    if (ss_tot <= 0) { loc_r2[i] <- 1; next }
    fit <- stats::lm.fit(cbind(1, xs), ys)
    loc_r2[i] <- 1 - sum(fit$residuals^2) / ss_tot
  }

  flag <- !is.na(loc_r2) & loc_r2 >= r2_thr
  # contiguous runs of TRUE
  r <- rle(flag)
  ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  keep <- which(r$values & r$lengths >= min_len)
  if (length(keep) == 0) return(NULL)
  do.call(rbind, lapply(keep, function(k) {
    idx <- starts[k]:ends[k]
    data.frame(rank_lo = rank[starts[k]], rank_hi = rank[ends[k]],
               npts = length(idx), r2 = mean(loc_r2[idx]))
  }))
}
