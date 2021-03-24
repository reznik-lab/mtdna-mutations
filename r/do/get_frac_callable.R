

source(here::here('r/prerequisites.R'))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate the callable-fraction of each patient's tumor and matched normal 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

frac_callable <- function(cov) {
    x <- as.integer(cov[1,])
    frac <- sum(x==1) / length(x)
    list(frac=frac)
}

cov_t <- fread(here('data/original_data/coverage_matrix_tcga_tumor_only.txt.gz'),header=T)
cov_n <- fread(here('data/original_data/coverage_matrix_tcga_normal_only.txt.gz'),header=T)
cov_tn <- fread(here('data/original_data/coverage_matrix_tcga.txt.gz'),header=T)

callable_t <- cov_t[,frac_callable(.SD),by=Tumor_Sample_Barcode]
callable_n <- cov_n[,frac_callable(.SD),by=Tumor_Sample_Barcode]
callable_tn <- cov_tn[,frac_callable(.SD),by=Tumor_Sample_Barcode]

callable_t$sample <- 'Tumor'
callable_n$sample <- 'Normal'
callable_tn$sample <- 'Tumor+Normal'
callable <- rbind(callable_t, callable_n, callable_tn)
callable$min_alt <- 'Min depth = 5'

cov20_t <- fread(here('data/original_data/coverage_matrix_tcga_tumor_only20.txt.gz'),header=T)
cov20_n <- fread(here('data/original_data/coverage_matrix_tcga_normal_only20.txt.gz'),header=T)
cov20_tn <- fread(here('data/original_data/coverage_matrix_tcga20.txt.gz'),header=T)

callable20_t <- cov20_t[,frac_callable(.SD),by=Tumor_Sample_Barcode]
callable20_n <- cov20_n[,frac_callable(.SD),by=Tumor_Sample_Barcode]
callable20_tn <- cov20_tn[,frac_callable(.SD),by=Tumor_Sample_Barcode]

callable20_t$sample <- 'Tumor'
callable20_n$sample <- 'Normal'
callable20_tn$sample <- 'Tumor+Normal'
callable20 <- rbind(callable20_t, callable20_n, callable20_tn)
callable20$min_alt <- 'Min depth = 20'

callable_dat <- rbind(callable, callable20)
write.tsv(callable_dat,here('data/processed_data/frac_callable.txt'))


