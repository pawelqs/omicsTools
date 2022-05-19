
test_that("Volcano plots does not throw error", {
  p <- plot_volcano(airway_deseq_res, label_col = "gene_id", labels = "ENSG00000179593")
  expect_identical(class(p), c("gg", "ggplot"))
})
