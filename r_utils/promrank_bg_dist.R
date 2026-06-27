# exploratory: what distribution do the BACKGROUND peak prominences follow?
# pools the prominence values across the given cells (excluding ranks 1-2 and
# any prominence > bg_prom_max), fits candidate distributions per cell with
# MASS::fitdistr, and ranks them by AIC. returns the per-cell AIC table and
# prints which distribution wins (lowest mean dAIC / most often best).
# Depends: .get_cell_row().
promrank_bg_dist <- function(cell_ids, profiles_df, suffix = "-1",
                             min_rank = 3, bg_prom_max = 0.1,
                             dists = c("exponential", "lognormal",
                                       "gamma", "weibull"),
                             plot = TRUE) {
  require(MASS)

  # background prominence values for one cell
  get_vals <- function(cid) {
    row <- .get_cell_row(profiles_df, cid, suffix)
    if (is.null(row)) return(numeric(0))
    pc <- grep("^peak[0-9]+_prom$", names(row), value = TRUE)
    pc <- pc[order(as.integer(gsub("[^0-9]", "", pc)))]
    pv <- as.numeric(row[1, pc])
    pv <- pv[is.finite(pv) & pv > 0]
    pv[seq_along(pv) >= min_rank & pv <= bg_prom_max]
  }
  per_cell <- setNames(lapply(cell_ids, get_vals), cell_ids)

  fit_aic <- function(v) {
    if (length(v) < 5) return(NULL)
    out <- sapply(dists, function(dn) {
      f <- tryCatch(suppressWarnings(MASS::fitdistr(v, dn)),
                    error = function(e) NULL)
      if (is.null(f)) NA_real_ else stats::AIC(f)
    })
    out
  }

  rows <- list()
  for (cid in cell_ids) {
    a <- fit_aic(per_cell[[cid]])
    if (is.null(a)) next
    rows[[length(rows) + 1]] <- data.frame(
      cell = cid, dist = names(a), aic = as.numeric(a),
      n = length(per_cell[[cid]]), row.names = NULL)
  }
  aic_df <- do.call(rbind, rows)
  if (is.null(aic_df)) { message("no cell had >= 5 background peaks"); return(invisible(NULL)) }

  # dAIC = AIC - best AIC within each cell (0 = best fit for that cell)
  aic_df$dAIC <- ave(aic_df$aic, aic_df$cell,
                     FUN = function(x) x - min(x, na.rm = TRUE))

  best <- do.call(rbind, lapply(split(aic_df, aic_df$cell),
                  function(s) s[which.min(s$aic), ]))

  summary_df <- aggregate(dAIC ~ dist, aic_df, function(x) mean(x, na.rm = TRUE))
  summary_df$n_best <- as.integer(table(best$dist)[summary_df$dist])
  summary_df$n_best[is.na(summary_df$n_best)] <- 0L
  summary_df <- summary_df[order(summary_df$dAIC), ]

  cat("Background prominence distribution fit (",
      length(unique(aic_df$cell)), " cells):\n", sep = "")
  print(summary_df, row.names = FALSE)
  cat("\nBest overall (lowest mean dAIC): ", summary_df$dist[1], "\n", sep = "")

  if (plot) {
    pp <- ggplot(aic_df, aes(reorder(dist, dAIC), dAIC)) +
      geom_boxplot(outlier.size = 0.6, fill = "grey90") +
      geom_jitter(width = 0.12, size = 1, colour = "grey40") +
      labs(x = NULL, y = "ΔAIC vs best (per cell)",
           title = "Background prominence: distribution fit") +
      theme_minimal()
    print(pp)
  }

  invisible(list(aic = aic_df, best_per_cell = best,
                 summary = summary_df, vals = per_cell))
}
