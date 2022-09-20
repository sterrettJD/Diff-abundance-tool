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
# total_reads <- sum(gt_data)
# microbe_read_sums <- colSums(gt_data)
# sum(microbe_read_sums > 0.005*total_reads) # 29 microbes have > 0.5% of the total reads
# gt_data <- gt_data[,which(microbe_read_sums > 0.005*total_reads)]


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
x <- model.matrix(~Germ_free, data=cbind(gt_data, metadata)) %>% 
    as.data.frame()
p <- ncol(x-1) #subtract 1 bc of the intercept

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
    x=x,
    y=y
)


fit1 <- stan(
    file = "diff-abund-NB.stan",  # Stan program
    data = stan_data,    # named list of data
    chains = 4,# number of Markov chains
    warmup = 500,          # number of warmup iterations per chain
    iter = 2000,            # total number of iterations per chain
    cores = 4,# number of cores (could use one per chain)
    refresh = 1, # progress shown
)

# assess sampler params
sampler_params <- get_sampler_params(fit1, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)

# print some of the output params, phi and beta
#print(fit1, pars = c("phi", "beta"))
#plot(fit1, pars=c("beta"))

# extract betas from the model
model_betas <- rstan::extract(fit1, permuted=T)$beta[,2,]
dim(model_betas)

# create a function to convert from alr to clr for interpretation
alr2clr <- function(x){
    d <- nrow(x)
    zeros <- matrix(0, nrow=d, ncol=1)
    mat <- cbind(x, zeros)
    x_clr <- mat - rowMeans(mat)
    return(x_clr)
}

model_betas_clr <- alr2clr(model_betas)

# check that this does the same as using alrInv to clr
# either of these options should work
model_betas_clr2 <- compositions::alrInv(model_betas) %>% compositions::clr()
rbind(head(model_betas_clr[1,]),
      head(model_betas_clr2[1,]))

model_betas_clr <- model_betas_clr %>% as.data.frame()
colnames(model_betas_clr) <- colnames(gt_data)

# rename the columns with taxonomy
tax <- read_qza("../../Reisdorph/GT-micro-metabo/microbiome/taxonomy_GT.qza")$data %>%
    as.data.frame()

model_betas_clr_asv <- model_betas_clr
colnames(model_betas_clr) <- sapply(
    colnames(model_betas_clr),
    function(x){tax[tax$Feature.ID==x,"Taxon"]}
    )

microbe_beta_means_asv <- colMeans(model_betas_clr_asv)
microbe_beta_sd_asv <- apply(model_betas_clr_asv, 2, sd)

beta_plot_data <- data.frame(
    microbe_beta_means_asv, microbe_beta_sd_asv, 
    row.names=colnames(model_betas_clr_asv))

beta_plot_data$tax <- colnames(model_betas_clr)
beta_plot_data <- beta_plot_data %>% 
    rename(mean=microbe_beta_means_asv,
                          sd=microbe_beta_sd_asv)
beta_plot_data$nameshort <- beta_plot_data$tax %>% 
    sapply(function(x){str_split(x, pattern=";", n=4)[[1]][4]})
beta_plot_data$asv <- rownames(beta_plot_data)

beta_plot_data <- beta_plot_data[order(beta_plot_data$mean),]

pd <- position_dodge(0.1)
ggplot(beta_plot_data,
       aes(y=reorder(asv,mean), x=mean)) +
    geom_point(position=pd, stat="identity",
             colour='red') +
    geom_errorbar(aes(xmin=mean-sd, xmax=mean+sd),
                  position=pd, alpha=0.1) +
    theme(axis.text.y=element_blank()) +
    ylab("ASV") +
    xlab("Log(Germ Free/Control) + K")

head(beta_plot_data)
tail(beta_plot_data)


fit1_summary <- summary(fit1, probs=c(.025, .975), digits=3)

### - genus table ###########

# read in qiime count table
gt_data2 <- read_qza("../../Reisdorph/GT-micro-metabo/microbiome/L6_table_GT.qza")$data %>%
    as.data.frame() %>%
    t() %>% as.data.frame()

# model inputs
formula <- "Germ_free"
N <- nrow(gt_data2)
D <- ncol(gt_data2)

# Create the predictor vector
x <- model.matrix(~Germ_free, data=cbind(gt_data2, metadata)) %>% 
    as.data.frame()
p <- ncol(x-1) #subtract 1 bc of the intercept

# calculate sequencing depth for each sample
depth <- gt_data2 %>% rowSums()

# The count table is our outcome y
y <- gt_data2

#concatenate to a named list for stan
stan_data2 <- list(
    N=N,
    D=D,
    p=p,
    depth=depth,
    x=x,
    y=y
)


fit2 <- stan(
    file = "diff-abund-NB.stan",  # Stan program
    data = stan_data2,    # named list of data
    chains = 4,# number of Markov chains
    warmup = 500,          # number of warmup iterations per chain
    iter = 2000,            # total number of iterations per chain
    cores = 4,# number of cores (could use one per chain)
    refresh = 1, # progress shown
)


# assess sampler params
sampler_params2 <- get_sampler_params(fit2, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params2), digits = 2)

# extract betas from the model
model_betas2 <- rstan::extract(fit2, permuted=T)$beta[,2,]
dim(model_betas2)

model_betas_clr2 <- alr2clr(model_betas2) %>% as.data.frame()
colnames(model_betas_clr2) <- colnames(gt_data2)


microbe_beta_means2 <- colMeans(model_betas_clr2)
microbe_beta_sd2 <- apply(model_betas_clr2, 2, sd)

beta_plot_data2 <- data.frame(
    mean=microbe_beta_means2, sd=microbe_beta_sd2, 
    row.names=colnames(model_betas_clr2))

beta_plot_data2$tax <- colnames(model_betas_clr2)
beta_plot_data2 <- beta_plot_data2[order(beta_plot_data2$mean),]

pd <- position_dodge(0.1)
ggplot(beta_plot_data2,
       aes(y=reorder(tax,mean), x=mean)) +
    geom_point(position=pd, stat="identity",
               colour='red') +
    geom_errorbar(aes(xmin=mean-sd, xmax=mean+sd),
                  position=pd, alpha=0.2) +
    theme(axis.text.y=element_blank()) +
    ylab("Genus") +
    xlab("Log(Germ Free/Control) + K")
