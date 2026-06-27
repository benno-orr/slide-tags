# (near-)homopolymer bead barcodes — TRUE if the barcode is within `max_mismatch`
# substitutions of a pure homopolymer (i.e. some base occupies >= L-max_mismatch
# of the L positions). Mirrors map_cells.py's HOMOPOLYMER_MAX_MISMATCH = 2 filter,
# which drops these beads before DBSCAN/KDE (poly-G etc. are SB artifacts).
.is_near_homopolymer <- function(bc, max_mismatch = 2) {
  bc <- as.character(bc)
  vapply(bc, function(s) {
    ch <- strsplit(s, "")[[1]]
    mx <- max(tabulate(match(ch, c("A", "C", "G", "T")), nbins = 4))
    (length(ch) - mx) <= max_mismatch
  }, logical(1), USE.NAMES = FALSE)
}

# drop near-homopolymer beads from a bead table (matches map_cells.py).
.drop_homopolymer_beads <- function(df, bc_col = "bead_bc", max_mismatch = 2) {
  if (!bc_col %in% names(df)) return(df)
  df[!.is_near_homopolymer(df[[bc_col]], max_mismatch), , drop = FALSE]
}
