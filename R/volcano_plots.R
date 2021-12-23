
#' Plot Volcano plot
#' @export
plot_volcano <- function(df, p_adj_col = "padj", lfc_col = "log2FoldChange",
                         label_col = "mgi_symbol",
                         p_adj_treshold = 0.05, lfc_treshold = 1) {

  df <- mutate(df,
               significant = !!sym(p_adj_col) < p_adj_treshold,
               change = case_when(
                 significant & !!sym(lfc_col) > lfc_treshold ~ "up",
                 significant & !!sym(lfc_col) < lfc_treshold ~ "down",
                 TRUE ~ "not sig."
               )
  )

  x_max <- df %>%
    pull(!!sym(lfc_col)) %>%
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
    summarise(n = n()) %>%
    full_join(annot_positions) %>%
    mutate(
      n = replace_na(n, 0),
      text = str_c(change, ": ", n))

  n_up <- filter(numbers, change == "up")$n
  n_down <- filter(numbers, change == "down")$n
  n_changed <- n_up + n_down
  label_n_up <- round((n_up * 60) / n_changed)
  label_n_down <- round((n_down * 60) / n_changed)

  labels_up <- df %>%
    filter(change == "up") %>%
    arrange(!!sym(p_adj_col)) %>%
    slice(1:label_n_up)
  labels_down <- df %>%
    filter(change == "down") %>%
    arrange(!!sym(p_adj_col)) %>%
    slice(1:label_n_down)
  labels <- bind_rows(labels_up, labels_down)

  ggplot(df, aes(!!sym(lfc_col), -log10(!!sym(p_adj_col)), color = change)) +
    geom_point() +
    scale_color_manual(
      values = c("steelblue2", "gray", "orangered3"),
      breaks = c("down", "not sig.", "up")
    ) +
    geom_label_repel(aes(label = !!sym(label_col)),
      data = labels, show.legend = FALSE
    ) +
    xlim(-x_max, x_max) +
    geom_hline(yintercept = -log10(p_adj_treshold), linetype = "dashed") +
    geom_vline(xintercept = lfc_treshold, linetype = "dashed") +
    geom_vline(xintercept = -lfc_treshold, linetype = "dashed") +
    geom_label(aes(x, y, vjust = y_just, hjust = x_just, label = text),
               data = numbers, color = "black", size = 5
    ) +
    theme(legend.position = "bottom")
}
