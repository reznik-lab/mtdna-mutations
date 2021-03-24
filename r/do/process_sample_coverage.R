# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate misc tables used in other scripts based on sample coverage data
# - table of number of samples covered per position (tumor normal and tumor-only data)
# - matrix of effective-gene-lengths for each sample (tumor normal and tumor-only data)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source(here::here('r/prerequisites.R'))


## get samples callable per position
process_coverage <- function(coverage_file, coverage_tumor_only_file, coverage_normal_only_file, chrM_annotated_file, samples_callable_per_pos_file, sample_effective_gene_lengths_file) {

    get_callable_per_pos <- function(covg) {
        m <- covg[,(valid_positions),with=F]
        N <- colSums(m,na.rm=T)
        samples_callable_per_pos <- data.table(pos=as.integer(names(N)),N=N)
        samples_callable_per_pos
    }
    chrm <- fread(chrM_annotated_file)
    valid_positions <- as.character(chrm$pos)
    covg <- fread(coverage_file,header=T)
    covg_tumoronly <- fread(coverage_tumor_only_file,header=T)
    covg_normalonly <- fread(coverage_normal_only_file,header=T)
    callable_tumornormal <- get_callable_per_pos(covg)
    callable_tumoronly <- get_callable_per_pos(covg_tumoronly)
    callable_normalonly <- get_callable_per_pos(covg_normalonly)
    callable <- merge(callable_tumornormal, callable_tumoronly, by='pos', all=T)
    names(callable) <- c('pos','callable_tn','callable_t')
    callable <- merge(callable, callable_normalonly, by='pos', all=T)
    setnames(callable,'N','callable_n')
    write.tsv(callable,samples_callable_per_pos_file)

    ## get effective gene length for each sample (length of gene for which mutations are callable in each sample)
    get_effective_gene_lengths_per_sample <- function(covg) {
        m <- covg[,(valid_positions),with=F]
        m <- t(m)
        colnames(m) <- covg$Tumor_Sample_Barcode
        m <- cbind(pos=as.integer(rownames(m)), adt(m))
        m <- merge(chrm[,c('symbol','pos'),with=F], m, by='pos', all.x=T)
        m[,pos:=NULL]
        sum_coverage_per_gene <- function(tmp) {
            out <- colSums(tmp, na.rm=T)    
            out <- data.table(Tumor_Sample_Barcode=names(out), effective_length=out)
            out
        }
        res <- m[,sum_coverage_per_gene(.SD),by=symbol]
        res <- reshape(res,idvar='Tumor_Sample_Barcode',timevar='symbol', direction='wide')
        res[is.na(res)] <- 0
        names(res) <- gsub('effective_length[.]','',names(res))
        fields <- c('Tumor_Sample_Barcode',sort(names(res)[2:ncol(res)]))
        res <- res[,(fields),with=F]
        res
    }

    effective_gene_lengths_tumornormal <- get_effective_gene_lengths_per_sample(covg)
    effective_gene_lengths_tumornormal$sample_type <- 'tumor+normal'
    effective_gene_lengths_tumoronly <- get_effective_gene_lengths_per_sample(covg_tumoronly)
    effective_gene_lengths_tumoronly$sample_type <- 'tumor only'
    effective_gene_lengths_normalonly <- get_effective_gene_lengths_per_sample(covg_normalonly)
    effective_gene_lengths_normalonly$sample_type <- 'normal only'
    effective_gene_lengths <- rbind(effective_gene_lengths_tumornormal, 
                                    effective_gene_lengths_tumoronly,
                                    effective_gene_lengths_normalonly)
    fields <- c('Tumor_Sample_Barcode','sample_type',grep('MT-',names(effective_gene_lengths),value=T))
    effective_gene_lengths <- effective_gene_lengths[,(fields),with=F]
    write.tsv(effective_gene_lengths, sample_effective_gene_lengths_file)
    system(paste('gzip',sample_effective_gene_lengths_file),wait=T)
}



## process sample coverage with min coverage of 5 reads
# input
coverage_file <- here('data/original_data/coverage_matrix_tcga.txt.gz')
coverage_tumor_only_file <- here('data/original_data/coverage_matrix_tcga_tumor_only.txt.gz')
coverage_normal_only_file <- here('data/original_data/coverage_matrix_tcga_normal_only.txt.gz')
chrM_annotated_file <- here('data/original_data/chrM_annotated.txt')

# output
samples_callable_per_pos_file <- here('data/processed_data/tcga_samples_callable_per_position.txt')
sample_effective_gene_lengths_file <- here('data/processed_data/tcga_sample_gene_lengths.txt')

process_coverage(coverage_file, coverage_tumor_only_file, coverage_normal_only_file, chrM_annotated_file, samples_callable_per_pos_file, sample_effective_gene_lengths_file)




## process sample coverage with min coverage of 20 reads
# input
coverage_file <- here('data/original_data/coverage_matrix_tcga20.txt.gz')
coverage_tumor_only_file <- here('data/original_data/coverage_matrix_tcga_tumor_only20.txt.gz')
coverage_normal_only_file <- here('data/original_data/coverage_matrix_tcga_normal_only20.txt.gz')
chrM_annotated_file <- here('data/original_data/chrM_annotated.txt')

# output
samples_callable_per_pos_file <- here('data/processed_data/tcga_samples_callable_per_position20.txt')
sample_effective_gene_lengths_file <- here('data/processed_data/tcga_sample_gene_lengths20.txt')

process_coverage(coverage_file, coverage_tumor_only_file, coverage_normal_only_file, chrM_annotated_file, samples_callable_per_pos_file, sample_effective_gene_lengths_file)






