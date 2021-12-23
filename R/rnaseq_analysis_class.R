
new_rnaseq_analysis <- function(x) {
  stopifnot(is.list(x))
  structure(x, class = "rnaseq_analysis")
}


#' Initialise RNAseq analysis
#'
#' @param deseq2 Named list of deseq2 results
#'
#' @return rnaseq_analysis object.
#' @export
#' @importFrom  magrittr %>%
init_rnaseq_analysis <- function(deseq2) {
  deseq2 <- calc_MSD(deseq2)
  x <- list(
    tests = names(deseq2),
    de_res = deseq2,
    tmod_msets = list(),
    cerno_res = list()
  )
  new_rnaseq_analysis(x)
}


calc_MSD <- function(de_res) {
  de_res %>%
    map(
      ~ mutate(.x,
        lower = log2FoldChange - lfcSE,
        higher = log2FoldChange + lfcSE,
        MSD = abs(log2FoldChange) - lfcSE
      )
    )
}


#' @export
#' @importFrom  magrittr %>%
add_tmod_mset.rnaseq_analysis <- function(rnaseq, ...) {
  new_msets <- list(...)
  to_add <- setdiff(names(new_msets), names(rnaseq$tmod_msets))
  to_replace <- intersect(names(new_msets), names(rnaseq$tmod_msets))

  if (length(to_replace) > 0)
    rnaseq$tmod_msets <- new_msets[to_replace]
  if (length(to_add) > 0)
    rnaseq$tmod_msets <- c(rnaseq$tmod_msets, new_msets[to_add])
  rnaseq
}


#' @export
#' @importFrom tmod tmodCERNOtest
run_tmod_cerno.rnaseq_analysis <- function(rnaseq) {
  rnaseq$cerno_res <- rnaseq$tmod_msets %>%
    map(function(mset) {
      arranged_gene_lists <- map(rnaseq$de_res, extract_gene_list)
      map(arranged_gene_lists, tmodCERNOtest, mset = mset)
    }
  )
  rnaseq
}


extract_gene_list <- function(deseq_table) {
  deseq_table %>%
    filter(symbol != "") %>%
    arrange(desc(MSD)) %>%
    pull(symbol)
}


#' @export
#' @importFrom  magrittr %>%
count_DEGs.rnaseq_analysis <- function(rnaseq, padj_treshold = 0.05, logFC_treshold = 0) {
  rnaseq$de_res %>%
    map(~filter(.x, padj < padj_treshold, abs(log2FoldChange) > logFC_treshold)) %>%
    map_int(nrow) %>%
    enframe(name = "test", value = "N_DEGs")
}


#' @export
#' @importFrom tmod tmodCERNOtest
plot_tmod_panelplot.rnaseq_analysis <- function(rnaseq,
                            mset_name = names(rnaseq$tmod_msets)[[1]],
                            tests = rnaseq$tests,
                            pie.style = "pie",
                            filter.rows.pval = 0.001, ...) {

  message(sprintf("Plotting tmodPanelPlot for "), mset_name)
  mset <- rnaseq$tmod_msets[[mset_name]]
  pie <- rnaseq$de_res %>%
    prepare_de_res() %>%
    prepare_pie(mset)
  tmodPanelPlot(
    rnaseq$cerno_res[[mset_name]][tests],
    pie = pie[tests],
    pie.style = pie.style,
    filter.rows.pval = filter.rows.pval,
    ...
  )
}


prepare_de_res <- function(de_res) {
  de_res <- de_res %>%
    map(drop_duplicated, "symbol")
  common_genes <- de_res %>%
    map(~pull(.x, symbol)) %>%
    reduce(intersect) %>%
    unique()
  de_res %>%
    map(~filter(.x, symbol %in% common_genes)) %>%
    map(~arrange(.x, symbol))
}


prepare_pie <- function(deseq2, mset, lfc_treshold = 0.5) {
  pie_lfc <- deseq2 %>%
    map_dfc("log2FoldChange") %>%
    map_df(replace_na, 0)
  pie_p <- deseq2 %>%
    map("padj") %>%
    map_df(replace_na, 1)
  tmod::tmodDecideTests(
    g = deseq2[[1]][["symbol"]],
    lfc = pie_lfc,
    pval = pie_p,
    mset = mset,
    lfc.thr = 0.5
  )
}


#' @export
save_evidence_plots.rnaseq_analysis <- function(rnaseq,
                                                modules_to_detail,
                                                out_dir = "images/") {
  to_save <- expand_grid(
    mset_name = names(rnaseq$tmod_msets),
    test_name = rnaseq$tests
  )

  pwalk(to_save, function(mset_name, test_name) {
    cli::cat_line("Saving ", test_name, " for ", mset_name)
    file <- str_c(out_dir, test_name, ".", mset_name, ".evidence.pdf")
    pdf(file)

    plot_evidence_plots(
      de_res = rnaseq$de_res[[test_name]],
      cerno_res = rnaseq$cerno_res[[mset_name]][[test_name]],
      modules = modules_to_detail[[mset_name]],
      mset = rnaseq$tmod_msets[[mset_name]]
    )

    dev.off()
  })
}


plot_evidence_plots <- function(de_res, cerno_res, modules, mset) {

  walk(modules, function(module_id) {

    module_genes <- mset$MODULES2GENES[[module_id]]
    gene_data <- de_res %>%
      filter(symbol %in% module_genes) %>%
      arrange(desc(MSD)) %>%
      mutate(
        symbol = parse_factor(symbol, levels = symbol),
        change = if_else(log2FoldChange > 0, "up", "down")
      )

    stats <- filter(cerno_res, ID == module_id)
    stats <- if (nrow(stats) == 0) {
      "not enriched"
    } else {
      sprintf("AUC: %.3f, adj.p: %.1e", stats$AUC, stats$adj.P.Val)
    }
    N = length(module_genes)
    n = nrow(gene_data)
    up = filter(gene_data, change == "up", padj < 0.05) %>% nrow()
    down = filter(gene_data, change == "down", padj < 0.05) %>% nrow()
    subtitle <- sprintf("%s, genes kept %d/%d, %d up, %d down", stats, n, N, up, down)

    p1 <- as.ggplot(function() {
      par(cex = 0.4)
      arranged_genes <- extract_gene_list(de_res)
      evidencePlot(arranged_genes, m = module_id, mset = mset, gene.labels = TRUE)
    })

    p2 <- gene_data %>%
      head(100) %>%
      ggplot(mapping = aes(symbol, log2FoldChange, color = change, alpha = log10(baseMean))) +
      geom_point() +
      geom_errorbar(aes(ymin = lower, ymax = higher), width = 0.1) +
      scale_color_manual(values = c(up = "red", down = "blue")) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(size = 4, angle = 90),
        legend.position = "top"
      ) +
      labs(caption = "(max 100 genes shown)")

    plot(
      p1 / p2 +
        plot_annotation(
          title = filter(mset$MODULES, ID == module_id)$Title,
          subtitle = subtitle
        ) +
        plot_layout(width = 10, heights = 5)
    )
  })
}


#' @export
#' @importFrom  cli cat_line
print.rnaseq_analysis <- function(x, ...) {
  cat_line("Object of rnaseq_analysis class\n", col = "gray")
  cat_line(length(x$de_res), " tests found:")
  cat_line(" ", names(x$de_res), col = "lavenderblush2")
  cat_line()

  cat_line(length(x$tmod_msets), " tmod msets found:")
  cat_line(" ", names(x$tmod_msets), col = "lavenderblush2")
  cat_line()

  cerno_not_run <- setdiff(
    names(x$tmod_msets),
    names(x$cerno_res)
  )
  if (length(cerno_not_run) != 0) {
    cat_line("Cerno test not run for", length(cerno_not_run), " tmod msets:")
    cat_line(" ", cerno_not_run, col = "tomato")
  } else {
    cat_line("Cerno run for all tmod msets", col = "chartreuse4")
  }

  invisible(x)
}
