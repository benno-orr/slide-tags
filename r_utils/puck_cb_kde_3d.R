# 3D persp surface of the σ-KDE bead density for one cell (matches map_cells.py's
# Gaussian-sum KDE). bead_df = pair step's cell-barcode_coords.csv (cell_bc_10x,
# nUMI, x_um, y_um). Rendered via ggplotify so it composes in wrap_plots().
puck_cb_kde_3d <- function(cell_id, bead_df,
                           suffix = "-1",
                           sigma = 40,
                           grid = 300,
                           theta = 35,
                           phi = 25,
                           clip_bottom = 0.05,
                           wire = FALSE, wire_n = 30, wire_col = "grey30",
                           title = NULL) {
  require(viridisLite)

  cb <- sub(paste0(suffix, "$"), "", cell_id)
  cb_col <- if ("cell_bc" %in% names(bead_df)) "cell_bc" else "cell_bc_10x"
  u_col  <- if ("nUMI" %in% names(bead_df)) "nUMI" else "u"

  df <- bead_df[bead_df[[cb_col]] == cb, , drop = FALSE]
  df <- df[is.finite(df$x_um) & is.finite(df$y_um), , drop = FALSE]

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
  Z <- as.matrix(Z)

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
    do.call(persp, c(list(gx, gy, Z, col = cols, border = NA, shade = 0.35), pargs))

    # optional coarse wireframe overlaid on the filled surface
    if (wire) {
      st <- max(1L, floor(grid / wire_n))
      ix <- seq(1L, length(gx), by = st); iy <- seq(1L, length(gy), by = st)
      par(new = TRUE)
      do.call(persp, c(list(gx[ix], gy[iy], Z[ix, iy], col = NA,
                            border = wire_col, lwd = 0.4), pargs))
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
