
`%not in%` <- function(x, y) {
  !('%in%'(x, y))
}


fill_na <- function(mat, val) {
  mat[is.na(mat)] <- val
  mat
}


drop_duplicated <- function(tbl, column) {
  duplicated_values <- tbl[duplicated(tbl[[column]]), ] %>%
    dplyr::pull(!!dplyr::sym(column))
  dplyr::filter(tbl, !!dplyr::sym(column) %not in% duplicated_values)
}
