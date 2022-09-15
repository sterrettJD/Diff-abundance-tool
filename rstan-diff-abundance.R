library(rstan)
library(tidyverse)
library(qiime2R)

getwd()
options(mc.cores = parallel::detectCores()-1)

# read in qiime count table
gt_data <- read_qza("../../Reisdorph/GT-micro-metabo/microbiome/table_GT.qza")$data %>%
    as.data.frame() %>%
    t() %>% as.data.frame()

# read in metadata
metadata <- data.table::fread("../../Reisdorph/GT-micro-metabo/metadata/mapping.tsv") %>%
    as.data.frame()

metadata$SampleID <- metadata$`#SampleID`
rownames(metadata) <- metadata$SampleID

# do all row names match?
length(match(rownames(metadata),
             rownames(gt_data))
) == nrow(metadata)

# make sure metadata and the count table are ordered the same
metadata <- metadata[rownames(gt_data),]




       