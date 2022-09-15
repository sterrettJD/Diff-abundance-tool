library(rstan)
library(tidyverse)
library(qiime2R)

getwd()
options(mc.cores = parallel::detectCores()-1)

# read in qiime count table
gt_data <- read_qza("../../Reisdorph/GT-micro-metabo/microbiome/table_GT.qza")$data %>%
    as.data.frame() %>%
    t() %>% as.data.frame()

# only grab a few of the columns since we're testing stan
total_reads <- sum(gt_data)
microbe_read_sums <- colSums(gt_data)
sum(microbe_read_sums > 0.005*total_reads) # 29 microbes have > 0.5% of the total reads
gt_data <- gt_data[,which(microbe_read_sums > 0.005*total_reads)]


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
    y=y
)


fit1 <- stan(
    file = "diff-abund-NB.stan",  # Stan program
    data = stan_data,    # named list of data
    chains = 4,             # number of Markov chains
    warmup = 500,          # number of warmup iterations per chain
    iter = 1000,            # total number of iterations per chain
    cores = 4,              # number of cores (could use one per chain)
    refresh = 1             # progress shown
)
