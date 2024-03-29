
#' Plot Volcano plot
#' @importFrom ggrepel geom_label_repel
#' @import dplyr
#' @import tibble
#' @import ggplot2
#' @importFrom rlang sym
#' @import tidyr
#' @importFrom stringr str_c
#' @export
plot_volcano <- function(df, p_adj_col = "padj", lfc_col = "log2FoldChange",
                         label_col = "mgi_symbol",
                         p_adj_treshold = 0.05, lfc_treshold = 1,
                         labels_n = 20, labels = NULL, labels_fun = ggrepel::geom_label_repel) {
  df <- mutate(df,
    significant = !!sym(p_adj_col) < p_adj_treshold,
    change = case_when(
      significant & !!sym(lfc_col) > lfc_treshold ~ "up",
      significant & !!sym(lfc_col) < -lfc_treshold ~ "down",
      TRUE ~ "not sig."
    )
  )

  x_max <- df[[lfc_col]] %>%
    abs() %>%
    max()
  x_max <- 1.05 * x_max

  annot_positions <- tibble(
    change = c("down", "not sig.", "up"),
    x = c(-Inf, 0, Inf),
    x_just = c(0, 0.5, 1),
    y = Inf,
    y_just = 1)

  numbers <- df %>%
    group_by(change) %>%
    summarize(n = n()) %>%
    full_join(annot_positions) %>%
    mutate(
      n = replace_na(n, 0),
      text = str_c(change, ": ", n))

  labels_data <- list()
  labels_data$empty <- slice(df, 0)
  labels_data$top <- if (labels_n > 0) {
    n_up <- filter(numbers, change == "up")$n
    n_down <- filter(numbers, change == "down")$n
    n_changed <- n_up + n_down
    label_n_up <- if (n_changed > 0) round((n_up * labels_n) / n_changed) else 0
    label_n_down <- if (n_changed > 0) round((n_down * labels_n) / n_changed) else 0

    labels_up <- df %>%
      filter(change == "up") %>%
      arrange(!!sym(p_adj_col)) %>%
      slice(1:label_n_up)
    labels_down <- df %>%
      filter(change == "down") %>%
      arrange(!!sym(p_adj_col)) %>%
      slice(1:label_n_down)
    bind_rows(labels_up, labels_down)
  }
  labels_data$chosen <- if (is.character(labels)) {
    filter(df, !!sym(label_col) %in% labels)
  }
  labels_data <- bind_rows(labels_data)

  ggplot(df, aes(!!sym(lfc_col), -log10(!!sym(p_adj_col)), color = change)) +
    geom_point() +
    scale_color_manual(
      values = c("steelblue2", "gray", "orangered3"),
      breaks = c("down", "not sig.", "up")
    ) +
    labels_fun(
      aes(label = !!sym(label_col)),
      data = labels_data, show.legend = FALSE, max.overlaps = 30
    ) +
    xlim(-x_max, x_max) +
    geom_hline(yintercept = -log10(p_adj_treshold), linetype = "dashed") +
    geom_vline(xintercept = lfc_treshold, linetype = "dashed") +
    geom_vline(xintercept = -lfc_treshold, linetype = "dashed") +
    geom_label(
      aes(x, y, vjust = y_just, hjust = x_just, label = text),
      data = numbers, color = "black", size = 5
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}
