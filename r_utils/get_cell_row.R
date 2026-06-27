# Locate one cell's row in a profiles table by cell id (rownames or any of the
# common id columns), trying both with and without the sample suffix.
.get_cell_row <- function(profiles_df, cell_id, suffix = "-1") {
  cb1 <- cell_id
  cb2 <- sub(paste0(suffix, "$"), "", cell_id)

  # possible ID columns
  id_cols <- intersect(
    c("cell_id", "cell_bc", "cell_bc_10x", "barcode", "cb"),
    names(profiles_df)
  )

  # try rownames first
  rn <- rownames(profiles_df)
  hit <- which(rn %in% c(cb1, cb2))

  if (length(hit) == 0 && length(id_cols) > 0) {
    for (id_col in id_cols) {
      vals <- as.character(profiles_df[[id_col]])

      hit <- which(vals %in% c(cb1, cb2))

      # also try removing suffix from table IDs
      if (length(hit) == 0) {
        vals2 <- sub(paste0(suffix, "$"), "", vals)
        hit <- which(vals2 %in% c(cb1, cb2))
      }

      if (length(hit) > 0) break
    }
  }

  if (length(hit) == 0) return(NULL)

  profiles_df[hit[1], , drop = FALSE]
}
