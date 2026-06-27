puck_cb_scatter = function(cell_id, bead_df, suffix = "-1", eps = 50, title = NULL, min_samples = c(10, 9, 8)) {
  require(dbscan)
  # bead_df: pair step's cell-barcode_coords.csv with columns cell_bc_10x,
  # bead_bc, nUMI, x_um, y_um — x/y already joined from the puck.

  cb        <- sub(paste0(suffix, "$"), "", cell_id)

  cb_col <- if ("cell_bc" %in% names(bead_df)) "cell_bc" else "cell_bc_10x"
  df <- bead_df[bead_df[[cb_col]] == cb, , drop = FALSE]
  # drop beads with no coords (matches map_cells.py, which never positions them)
  df <- df[is.finite(df$x_um) & is.finite(df$y_um), , drop = FALSE]
  if (nrow(df) == 0) return(ggplot() + ggtitle(paste0(cell_id, " — no beads")))

  xy <- as.matrix(df[, c("x_um", "y_um")])

  # DBSCAN with descending min_samples escalation; cluster 0 = noise.
  # UMI-weighted to match the mapping pipeline: minPts is compared to the summed
  # nUMI in each neighbourhood, so a few high-UMI beads can form a cluster.
  lbls <- rep(0L, nrow(df))
  for (ms in min_samples) {
    res <- dbscan::dbscan(xy, eps = eps, minPts = ms, weights = df$nUMI)
    if (any(res$cluster > 0)) { lbls <- res$cluster; break }
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

  title = if (is.null(title)) cell_id else title

  p +
    coord_fixed() +
    theme(
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()) +

    ggtitle(title)
}
