
library(airway)
library(DESeq2)
library(tidyverse)

data("airway")

se <- airway
dds <- DESeqDataSet(se, design = ~ cell + dex)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds)
airway_deseq_res <- res %>%
  as.data.frame() %>%
  rownames_to_column("gene_id")

usethis::use_data(airway_deseq_res)
