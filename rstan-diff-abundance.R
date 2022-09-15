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
metadata$Germ_free <- metadata$PID == "Control"

# do all row names match?
length(match(rownames(metadata),
             rownames(gt_data))
) == nrow(metadata)

# make sure metadata and the count table are ordered the same
metadata <- metadata[rownames(gt_data),]


# model inputs
formula <- "Germ_free"
N <- nrow(gt_data)
D <- ncol(gt_data)

# Create the predictor vector
X <- model.matrix(~Germ_free, data=cbind(gt_data, metadata)) %>% 
    as.data.frame()
p <- ncol(X-1) #subtract 1 bc of the intercept

# calculate sequencing depth for each sample
depth <- gt_data %>% rowSums()

# The count table is our outcome y
y <- gt_data

#concatenate to a named list for stan
stan_data <- list(
    N=N,
    D=D,
    p=p,
    depth=depth,
    X=X,
    y=y,
)

