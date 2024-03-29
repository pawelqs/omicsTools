#' R function to create a genom-wide heatmap for one meta column for a list of GRanges
#'
#' @param granges .
#' @param meta_field .
#' @param row_groups .
#' @param keep_sites_present_in .
#' @param color_breaks .
#' @param colors .
#' @param window_width .
#' @param upper_limit .
#' @param cluster_rows .
#' @param show_row_names .
#' @param show_column_names .
#' @param use_raster .
#' @param cluster_columns .
#' @param border .
#' @param legend_params .
#' @param ...
#'
#' @return ComplexHeatmap heatmap
#'
#' @import dplyr
#' @import tibble
#' @import purrr
#' @import ComplexHeatmap
#'
#' @export
#'
#' @examples
#' .
heatmap_granges <- function(granges, meta_field,
                            row_groups = NULL,
                            keep_sites_present_in = floor(0.8*length(granges)),
                            color_breaks = c(0, 2, 6),
                            colors = c("dodgerblue3", "white", "firebrick3"),
                            window_width = 1000000, upper_limit = 6,
                            cluster_rows = FALSE,
                            show_row_names = TRUE,
                            show_column_names = FALSE,
                            use_raster = TRUE,
                            cluster_columns = FALSE,
                            border = TRUE,
                            legend_params = NULL, ...) {

  ########## Get common ranges
  ranges_cov <- granges %>%
    GenomicRanges::GRangesList() %>%
    unlist %>%
    plyranges::compute_coverage()

  range_fractions <- ranges_cov %>%
    as_tibble() %>%
    group_by(score) %>%
    summarise(width = sum(width)) %>%
    arrange(desc(score)) %>%
    mutate(
      total_len = sum(width),
      frct = width / total_len,
      cum = cumsum(frct)
    ) %>%
    select(score, cum) %>%
    deframe()

  print("Proportion of ranges present in N GRanges:")
  print(range_fractions)
  print(sprintf("Kept sites present in %d/%d GRanges", keep_sites_present_in, length(granges)))
  print(sprintf("Kept %f of all sites", range_fractions[as.character(keep_sites_present_in)]))

  common_ranges <- ranges_cov %>%
    filter(score > keep_sites_present_in) %>%
    plyranges::reduce_ranges()

  ########### List of GRanges to matrix
  windows <- common_ranges %>%
    tile(width = window_width) %>%
    unlist()

  dt <- map(granges, ~plyranges::join_overlap_left(windows, .x))
  mat <- dt %>%
    map(as_tibble) %>%
    # Summarise values for window ranges covering more than one interval
    map(~select(.x, seqnames, start, end, !!sym(meta_field))) %>%
    map(~group_by(.x, seqnames, start, end)) %>%
    map(~summarise(.x, val = mean(!!sym(meta_field)))) %>%
    # Create matrix
    map("val") %>%
    as.data.frame() %>%
    as.matrix() %>%
    t()
  mat[mat > upper_limit] <- upper_limit

  ######### Prepare heatmap chromosome annotation
  chr_means <- as_tibble(windows) %>%
    mutate(i = row_number()) %>%
    group_by(seqnames) %>%
    summarise(m = mean(i)) %>%
    deframe()
  chr_labels <- vector(mode = "character",length = length(windows))
  chr_labels[chr_means] <- names(chr_means)

  seqname_replacements <- names(chr_means) %>%
    as_tibble_col(column_name = "seqname") %>%
    mutate(code = row_number() %% 2)
  chr_bin <- tibble(seqname = as.character(seqnames(windows))) %>%
    left_join(seqname_replacements) %>%
    pull(code)


  chr_bar <- HeatmapAnnotation(
    chr_text = anno_text(chr_labels, gp = gpar(fontsize = 8)),
    chr = chr_bin,
    show_legend = FALSE,
    which = "column",
    col = list(chr = c("0" = "grey88", "1" = "black")
    ))

  ########### Plot Heatmap
  Heatmap(mat,
          col = circlize::colorRamp2(breaks = color_breaks, colors),
          row_split = if(is.null(row_groups)) NULL else factor(row_groups, levels = unique(row_groups)),
          cluster_row_slices = FALSE,
          top_annotation = chr_bar,
          heatmap_legend_param = c(list(title = meta_field), legend_params),
          cluster_rows = cluster_rows,
          show_row_names = show_row_names,
          show_column_names = show_column_names,
          use_raster = use_raster,
          cluster_columns = cluster_columns,
          border = border,
          ...
  )
}
