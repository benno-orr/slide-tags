puck_cb_scatter = function(cell_id, bead_df, profiles_df = NULL,
                           mode = c("raw", "norm"),
                           suffix = "-1", eps = 50, title = NULL,
                           min_samples = c(10, 9, 8), norm_to = 1000) {
  require(dbscan)
  mode <- match.arg(mode)
  # bead_df: pair step's cell-barcode_coords.csv with columns cell_bc_10x,
  # bead_bc, nUMI, x_um, y_um — x/y already joined from the puck.
  # mode "raw"  : per-bead nUMI weights, descending min_samples escalation.
  # mode "norm" : w = nUMI/Σ nUMI × norm_to, single DBSCAN at the stored
  #               norm_minpts (from profiles_df), else a high→low scan picking
  #               the highest minPts giving ≥2 clusters — matching map_cells.py.

  cb        <- sub(paste0(suffix, "$"), "", cell_id)

  cb_col <- if ("cell_bc" %in% names(bead_df)) "cell_bc" else "cell_bc_10x"
  df <- bead_df[bead_df[[cb_col]] == cb, , drop = FALSE]
  # drop beads with no coords (matches map_cells.py, which never positions them)
  df <- df[is.finite(df$x_um) & is.finite(df$y_um), , drop = FALSE]
  if (nrow(df) == 0) return(ggplot() + ggtitle(paste0(cell_id, " — no beads")))

  xy <- as.matrix(df[, c("x_um", "y_um")])
  minpts_used <- NA_real_

  if (mode == "raw") {
    # cluster 0 = noise; first min_samples tier with ≥1 cluster wins.
    lbls <- rep(0L, nrow(df))
    for (ms in min_samples) {
      res <- dbscan::dbscan(xy, eps = eps, minPts = ms, weights = df$nUMI)
      if (any(res$cluster > 0)) { lbls <- res$cluster; minpts_used <- ms; break }
    }
  } else {
    w <- df$nUMI / sum(df$nUMI) * norm_to
    # prefer the pipeline's stored norm_minpts for an exact match
    ms <- NA_real_
    if (!is.null(profiles_df)) {
      prow <- .get_cell_row(profiles_df, cell_id, suffix)
      if (!is.null(prow) && "norm_minpts" %in% names(prow))
        ms <- suppressWarnings(as.numeric(prow[["norm_minpts"]][1]))
    }
    if (!is.finite(ms)) {                       # fallback: high→low scan
      for (cand in seq(200, 1)) {
        res <- dbscan::dbscan(xy, eps = eps, minPts = cand, weights = w)
        if (length(unique(res$cluster[res$cluster > 0])) >= 2) { ms <- cand; break }
      }
    }
    lbls <- if (is.finite(ms))
      dbscan::dbscan(xy, eps = eps, minPts = ms, weights = w)$cluster
    else rep(0L, nrow(df))
    minpts_used <- ms
  }

  # c1 = largest cluster by summed UMI, c2+ = remaining clusters, c0 = noise
  df$dbscan_id <- lbls
  non_noise <- lbls[lbls > 0]
  if (length(non_noise) > 0) {
    umi_by_clust <- tapply(df$nUMI[lbls > 0], lbls[lbls > 0], sum)
    c1_id <- as.integer(names(which.max(umi_by_clust)))
    df$cluster_group <- ifelse(lbls == 0, "c0",
                         ifelse(lbls == c1_id, "c1", "c2"))
  } else {
    df$cluster_group <- "c0"
  }

  # per non-noise cluster: total SB UMI + centroid, to annotate on the plot.
  # labelled c1 (largest) then c2, c3, ... by descending UMI.
  clust_lab <- NULL
  if (length(non_noise) > 0) {
    ids <- unique(lbls[lbls > 0])
    clust_lab <- do.call(rbind, lapply(ids, function(id) {
      sub_df <- df[df$dbscan_id == id, , drop = FALSE]
      data.frame(
        dbscan_id = id,
        umi       = sum(sub_df$nUMI),
        x_um      = mean(sub_df$x_um),
        y_um      = mean(sub_df$y_um)
      )
    }))
    clust_lab <- clust_lab[order(-clust_lab$umi), , drop = FALSE]
    clust_lab$name <- NA_character_
    clust_lab$name[clust_lab$dbscan_id == c1_id] <- "c1"
    rest <- which(clust_lab$dbscan_id != c1_id)
    clust_lab$name[rest] <- paste0("c", seq_along(rest) + 1)
    # nudge label to the right of the cluster (~4% of the plotted x-span)
    x_off <- 0.04 * diff(range(df$x_um))
    clust_lab$x_lab <- clust_lab$x_um + x_off
  }

  require(ggnewscale)

  # integer-only breaks for the (log10) UMI colourbars
  int_log_breaks <- function(lims) {
    b <- scales::breaks_log(5)(lims)
    b <- unique(round(b))
    b[b >= 1 & is.finite(b)]
  }

  # within each cluster group draw low-UMI beads first so the brightest sit on top
  ord_umi <- function(d) d[order(d$nUMI), , drop = FALSE]
  d0 <- ord_umi(df[df$cluster_group == "c0", , drop = FALSE])
  d2 <- ord_umi(df[df$cluster_group == "c2", , drop = FALSE])
  d1 <- ord_umi(df[df$cluster_group == "c1", , drop = FALSE])

  # layer order: c0 background → c2 → c1 foreground; each group its own colourbar
  p <- ggplot()

  if (nrow(d0) > 0)
    p <- p +
      geom_point(data = d0, aes(x_um, y_um, colour = nUMI), size = 1.8) +
      scale_colour_gradient(low = "grey85", high = "grey25", trans = "log10",
                            name = "c0 nUMI",
                            breaks = int_log_breaks,
                            labels = scales::label_number(accuracy = 1)) +
      ggnewscale::new_scale_colour()

  if (nrow(d2) > 0)
    p <- p +
      geom_point(data = d2, aes(x_um, y_um, colour = nUMI), size = 1.8) +
      viridis::scale_colour_viridis(option = "plasma", trans = "log10",
                                    name = "c2+ nUMI",
                                    breaks = int_log_breaks,
                                    labels = scales::label_number(accuracy = 1)) +
      ggnewscale::new_scale_colour()

  if (nrow(d1) > 0)
    p <- p +
      geom_point(data = d1, aes(x_um, y_um, colour = nUMI), size = 1.8) +
      viridis::scale_colour_viridis(option = "viridis", trans = "log10",
                                    name = "c1 nUMI",
                                    breaks = int_log_breaks,
                                    labels = scales::label_number(accuracy = 1))

  # annotate each non-noise cluster with its total SB UMI, just to the right
  if (!is.null(clust_lab))
    p <- p + geom_text(
      data = clust_lab,
      aes(x = x_lab, y = y_um, label = paste0(name, ": ", umi, " UMI")),
      hjust = 0, size = 3.2, colour = "black"
    )

  title = if (is.null(title))
    sprintf("%s — %s dbscan (minPts=%s)", cell_id, mode,
            ifelse(is.finite(minpts_used), as.character(round(minpts_used, 1)), "NA"))
  else title

  p +
    coord_fixed() +
    theme(
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()) +

    ggtitle(title)
}
