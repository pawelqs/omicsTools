
#' Count DEGs
#'
#' @param rnaseq rna_analysis object
#'
#' @return Summary tibble
#' @export
count_DEGs <- function(rnaseq) {
  UseMethod("count_DEGs")
}


#' add_tmod_mset
#'
#' @param rnaseq rna_analysis object
#' @param ... name = value pairs of mset objects
#'
#' @return rnaseq_analysis object.
#' @export
add_tmod_mset <- function(rnaseq, ...) {
  UseMethod("add_tmod_mset")
}


#' Runs tmod CERNO test for each tmod mset
#'
#' @param rnaseq rnaseq_analysis object
#'
#' @return rnaseq_analysis object
#' @export
run_tmod_cerno <- function(rnaseq, ...) {
  UseMethod("run_tmod_cerno")
}


#' Plots tmod pamelPlot
#'
#' @param rnaseq rnaseq_analysis object
#'
#' @return rnaseq_analysis object
#' @export
plot_tmod_panelplot <- function(rnaseq, mset_name = NULL, tests = rnaseq$tests,
                                pie.style = "pie",
                                filter.rows.pval = 0.001, ...) {
  UseMethod("plot_tmod_panelplot")
}


#' Plots tmod evidence plots
#'
#' @param rnaseq rnaseq_analysis object
#' @param modules_to_detail list(mset_name = vec_of_modules)
#'
#' @export
save_evidence_plots <- function(rnaseq, modules_to_detail, out_dir = "images/") {
  UseMethod("save_evidence_plots")
}
