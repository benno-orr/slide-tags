# 3D persp surface of the σ-KDE bead density for one cell (matches map_cells.py's
# Gaussian-sum KDE). bead_df = pair step's cell-barcode_coords.csv (cell_bc_10x,
# nUMI, x_um, y_um). Rendered via ggplotify so it composes in wrap_plots().
# If profiles_df (cell_coords.csv) is given, peaks 1-9 (peakK_x/peakK_y) are
# marked on the surface as circles coloured by prominence rank (plasma; rank 1 =
# highest prom = brightest). Depends: .get_cell_row() when profiles_df is used.
puck_cb_kde_3d <- function(cell_id, bead_df,
                           profiles_df = NULL,
                           suffix = "-1",
                           sigma = 40,
                           grid = 500,
                           theta = 35,
                           phi = 25,
                           clip_bottom = 0.05,
                           wire = FALSE, wire_n = 30, wire_col = "grey30",
                           n_peaks = 9, peak_cex = 1.6, peak_labels = TRUE,
                           title = NULL) {
  require(viridisLite)

  cb <- sub(paste0(suffix, "$"), "", cell_id)
  cb_col <- if ("cell_bc" %in% names(bead_df)) "cell_bc" else "cell_bc_10x"
  u_col  <- if ("nUMI" %in% names(bead_df)) "nUMI" else "u"

  df <- bead_df[bead_df[[cb_col]] == cb, , drop = FALSE]
  df <- df[is.finite(df$x_um) & is.finite(df$y_um), , drop = FALSE]
  # drop (near-)homopolymer beads so the surface matches map_cells.py's KDE peaks
  df <- .drop_homopolymer_beads(df)

  x <- df$x_um
  y <- df$y_um
  u <- as.numeric(df[[u_col]])
  u[!is.finite(u)] <- 0

  pad <- 3 * sigma
  gx <- seq(min(x) - pad, max(x) + pad, length.out = grid)
  gy <- seq(min(y) - pad, max(y) + pad, length.out = grid)

  Ax <- exp(-outer(gx, x, "-")^2 / (2 * sigma^2))
  Ay <- exp(-outer(gy, y, "-")^2 / (2 * sigma^2))
  Z <- (Ax %*% t(sweep(Ay, 2, u, "*"))) / (2 * pi * sigma^2)
  Z <- as.matrix(Z)                          # Z[i, j] = density at (gx[i], gy[j])

  Zfacet <- (
    Z[-1, -1] +
      Z[-1, -ncol(Z)] +
      Z[-nrow(Z), -1] +
      Z[-nrow(Z), -ncol(Z)]
  ) / 4

  z_floor <- quantile(Zfacet, clip_bottom, na.rm = TRUE)
  Zcolor <- pmax(Zfacet, z_floor)

  # use only the top 90% of the viridis scale (skip the darkest 10%)
  cols <- viridisLite::viridis(256, begin = 0.1, end = 1)[
    cut(Zcolor, breaks = 256, include.lowest = TRUE)
  ]

  # peaks 1..n_peaks from the profiles table, positioned on the surface (z via
  # bilinear interpolation of Z at the peak xy), coloured by prominence rank.
  peaks <- NULL
  if (!is.null(profiles_df)) {
    row <- .get_cell_row(profiles_df, cell_id, suffix)
    if (!is.null(row)) {
      K <- seq_len(n_peaks)
      cx <- paste0("peak", K, "_x"); cy <- paste0("peak", K, "_y")
      have <- cx %in% names(row) & cy %in% names(row)
      K <- K[have]
      if (length(K) > 0) {
        px <- suppressWarnings(as.numeric(row[1, paste0("peak", K, "_x")]))
        py <- suppressWarnings(as.numeric(row[1, paste0("peak", K, "_y")]))
        ok <- is.finite(px) & is.finite(py) &
              px >= min(gx) & px <= max(gx) & py >= min(gy) & py <= max(gy)
        K <- K[ok]; px <- px[ok]; py <- py[ok]
        if (length(K) > 0) {
          interp_z <- function(xq, yq) {
            ix <- pmin(pmax(findInterval(xq, gx), 1), length(gx) - 1)
            iy <- pmin(pmax(findInterval(yq, gy), 1), length(gy) - 1)
            tx <- (xq - gx[ix]) / (gx[ix + 1] - gx[ix])
            ty <- (yq - gy[iy]) / (gy[iy + 1] - gy[iy])
            Z[cbind(ix, iy)]         * (1 - tx) * (1 - ty) +
            Z[cbind(ix + 1, iy)]     *      tx  * (1 - ty) +
            Z[cbind(ix, iy + 1)]     * (1 - tx) *      ty  +
            Z[cbind(ix + 1, iy + 1)] *      tx  *      ty
          }
          # plasma by prom rank: rank 1 (highest prom) = brightest
          pal <- rev(viridisLite::plasma(n_peaks))
          peaks <- data.frame(rank = K, x = px, y = py,
                              z = interp_z(px, py), col = pal[K])
        }
      }
    }
  }

  # device-independent 3D→2D projection (avoids depending on trans3d's location)
  proj3d <- function(x, y, z, pmat) {
    tr <- cbind(x, y, z, 1) %*% pmat
    list(x = tr[, 1] / tr[, 4], y = tr[, 2] / tr[, 4])
  }

  # render the base-graphics persp into a grob/ggplot so it COMPOSES (e.g. in
  # wrap_plots) instead of painting over the active device. as.ggplot runs the
  # draw call on its own offscreen device, so global par() is never touched.
  draw <- function() {
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    par(
      bg = "grey95",
      mar = c(0, 0, 0, 0),
      mai = c(0, 0, 0, 0),
      oma = c(0, 0, 0, 0),
      xaxs = "i",
      yaxs = "i",
      plt = c(0, 1, 0, 1)
    )
    # zoom in ~20% by shrinking the plotted x/y box about its centre
    zoom <- 1.21
    cx <- mean(range(gx)); hx <- diff(range(gx)) / 2 / zoom
    cy <- mean(range(gy)); hy <- diff(range(gy)) / 2 / zoom
    # shared zlim so the filled surface and the wireframe use identical vertical
    # scaling (otherwise persp auto-scales each to its own z-range and they drift)
    pargs <- list(
      theta = theta, phi = phi, expand = 0.7,
      xlim = cx + c(-hx, hx), ylim = cy + c(-hy, hy), zlim = range(Z, finite = TRUE),
      box = FALSE, axes = FALSE, xlab = "", ylab = "", zlab = "", main = ""
    )
    pmat <- do.call(persp, c(list(gx, gy, Z, col = cols, border = NA, shade = 0.35), pargs))

    # optional coarse wireframe overlaid on the filled surface
    if (wire) {
      st <- max(1L, floor(grid / wire_n))
      ix <- seq(1L, length(gx), by = st); iy <- seq(1L, length(gy), by = st)
      par(new = TRUE)
      pmat <- do.call(persp, c(list(gx[ix], gy[iy], Z[ix, iy], col = NA,
                                    border = wire_col, lwd = 0.4), pargs))
    }

    # mark peaks 1..n on the surface (lifted slightly so they sit on top)
    if (!is.null(peaks) && nrow(peaks) > 0) {
      lift <- 0.03 * diff(range(Z, finite = TRUE))
      tp <- proj3d(peaks$x, peaks$y, peaks$z + lift, pmat)
      points(tp$x, tp$y, pch = 21, bg = peaks$col, col = "black",
             lwd = 0.7, cex = peak_cex)
      if (peak_labels)
        text(tp$x, tp$y, labels = peaks$rank, pos = 3, offset = 0.25,
             cex = 0.55, col = "grey15")
    }
  }

  ggplotify::as.ggplot(draw) +
    ggtitle(title) +
    theme(
      aspect.ratio = 1,        # square panel (coord_fixed equivalent)
      plot.margin = margin(0, 0, 0, 0),
      plot.background  = element_rect(fill = "grey95", colour = NA),
      panel.background = element_rect(fill = "grey95", colour = NA)
    )
}
