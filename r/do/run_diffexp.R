## is there a difference in gene expression between truncating +/- samples?
## if so, what about for hotspot samples? is there a difference?
## do noncanonical cancertypes with mt-loss phenocopy (renal, thyroid, colorectal)


source(here::here('r/prerequisites.R'))
require(DESeq2)
require(fgsea)
require(apeglm)

rsem_file <- here('data/ext/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv.gz')
hallmark_geneset_file <- here('data/ext/h.all.v7.1.symbols.gmt')

## create the expected directories for output
dirs <- here(c(
          'data/processed_data/rnaseq',
          'data/processed_data/rnaseq/truncating',
          'data/processed_data/rnaseq/vus',
          'data/processed_data/rnaseq/r25q')

setup_dir <- function(directory) {
    if(!dir.exists(directory)) {
        message('Creating directory: ',directory)
        dir.create(directory)
    } else {
        message(directory,' already exists.')
    }
}
hide <- lapply(dirs, setup_dir)


## load prerequisite data
sample_hotspot_status_data <- fread(here('data/processed_data/sample_hotspot_status_tcga_precalculated.txt'))
clin_data <- fread(here('data/original_data/data_samples_tcga.txt.gz'),select=c('PATIENT_ID','rna_available'))

## load/format the RSEM rna-seq data (download from GDC, see: github wiki)
message('loading/formating rnaseq RSEM as counts ...')
cnt_data <- fread(rsem_file)
cnt_data <- cnt_data[grepl('[?]',gene_id)==F,]
genes <- cnt_data$gene_id
cnt_data <- as.matrix(round(cnt_data[,2:ncol(cnt_data),with=F]))
rownames(cnt_data) <- genes
colnames(cnt_data) <- strtrim (colnames(cnt_data),12)
cnt_data[cnt_data < 0] <- NA
rowmeans <- rowMeans(cnt_data)
excluded_genes <- names(rowmeans[is.na(rowmeans)])
cnt_data <- cnt_data[rownames(cnt_data) %nin% excluded_genes,]
cnt_data0 <- rowSums(cnt_data==0)
genes_always_with_counts <- names(cnt_data0)[cnt_data0==0]
cnt_data <- cnt_data[rownames(cnt_data) %in% genes_always_with_counts,]



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define functions used in multiple analyses
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## define main function to run deseq/fgsea for a combo of variant/cancertype
run_type_with_hallmarks <- function(cancer.type,cnt,coldata,rnaseqdir,formula=~affected,rerun=F) { 
    message(paste0('\n',cancer.type))
    msg <- tryCatch({
        outfile <- paste0(rnaseqdir,gsub(' ','_',tolower(cancer.type)),'.rds')
        if(file.exists(outfile) & rerun==F) {
            l <- readRDS(outfile)
            res <- l$data
            cancer.type <- l$type
            affected <- l$affected
            n_affected <- l$n_affected
            n_wt <- l$n_wt
        } else {
            if(cancer.type!='') {
                coldata2 <- coldata[coldata$maintype %in% cancer.type,]
            } else {
                coldata2 <- coldata
            }
            coldata2$affected <- coldata2$status=='affected'
            affected <- rownames(coldata2[coldata2$affected==T,])
            n_affected <- length(affected)
            n_wt <- nrow(coldata2) - n_affected
            cnt2 <- cnt[,rownames(coldata2)]
            dds <- DESeqDataSetFromMatrix(countData = cnt2, colData = coldata2, design = formula)
            dds <- DESeq(dds)
            resultsNames(dds) # lists the coefficients
            res <- lfcShrink(dds, coef="affectedTRUE", type="apeglm")
            res <- cbind(gene=rownames(res),adt(res))
            ## include all padj despite 'independentFiltering'
            res <- res[!is.na(pvalue),]
            res$padj <- p.adjust(res$pvalue,method='BH')
            res <- res[order(padj,decreasing=F),]
        } 
        res$stat <- -log10(res$padj)*sign(res$log2FoldChange)
        res$type <- cancer.type
        f <- function(s) strsplit(s,'[|]')[[1]][1]
        res$symbol <- sapply(res$gene, f)
        stat <- res$stat
        names(stat) <- res$symbol
        stat <- stat[names(stat)!='?']

        ## GSEA
        gsea <- run_fgsea(hallmark_geneset_file, stat, pathways, n_perm=1e5)
        gsea$type <- cancer.type
        out <- list(data=res,gsea=gsea,type=cancer.type,affected=affected, n_affected=n_affected,n_wt=n_wt)
        saveRDS(out,file=outfile) 
        'no problem'
    },error=function(e) {
        message(e)
        as.character(e)
    })
}


run_fgsea <- function(geneset_file,stat,pathways,n_perm) { 
    set.seed(42)
    message(geneset_file)
    pathways <- gmtPathways(geneset_file)
    gsea <- fgsea(pathways = pathways, stats = stat, minSize=10, maxSize=500, nperm=n_perm)
    gsea <- adt(gsea) 
    gsea <- gsea[order(padj,decreasing=F),]
    gsea$geneset <- geneset_file
    gsea
}

load_gsea_results <- function(f,my.geneset='h.all') {
    x <- readRDS(f)
    g <- x$gsea
    g$geneset <- my.geneset
    edgetostring <- function(l) paste(l,collapse=',') 
    g$genes <- sapply(g$leadingEdge, edgetostring)
    g[,leadingEdge:=NULL]
    g
}

load_data_results <- function(f) {
    x <- readRDS(f)
    d <- x$data
    d
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run rnaseq between wt and truncating
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# prepare coldata and cnts 
## - load matrix with genotypes and wild-type status for all tcga samples
info <- copy(sample_hotspot_status_data)
info$PATIENT_ID <- strtrim(info$Tumor_Sample_Barcode,12)
clin <- copy(clin_data)
info <- merge(clin, info, by='PATIENT_ID', all.x=T)
info <- info[rna_available==T,]
info$status <- as.character(NA)
info[annotation_postfilter=='Wildtype',status:='wildtype']
info[annotation_postfilter=='Truncating',status:='affected']

#info <- info[!is.na(status),]
coldata_test <- info[status %in% c('wildtype','affected'),]
coldata_test <- coldata_test[,c('PATIENT_ID','maintype','status'),with=F]
coldata_test$status <- factor(coldata_test$status, levels=c('wildtype','affected'))

## subset cnt and coldata_test for samples in both
cnt_test <- copy(cnt_data)
valid_samples <- sort(intersect(colnames(cnt_test),coldata_test$PATIENT_ID))
cnt_test <- cnt_test[,valid_samples]
coldata_test <- coldata_test[PATIENT_ID %in% valid_samples,]
coldata_test <- as.data.frame(coldata_test) 
rownames(coldata_test) <- coldata_test$PATIENT_ID
coldata_test$PATIENT_ID <- NULL

## get valid types for testing (require 5+ wt and affected per cancer type)
tbl <- xtabs(~maintype+status ,data=coldata_test)
tbl <- cbind(type=rownames(tbl),adt(as.data.frame.matrix(tbl)))
names(tbl) <- c('type','wt','affected')
tbl <- tbl[affected >= 5 & wt >= 5]
tbl$N <- tbl$wt + tbl$affected 
tbl <- tbl[order(N,decreasing=F),]
types <- tbl$type
write.tsv(tbl,here('data/processed_data/rnaseq_truncating_histogram.txt'))

## test each cancer type with hallmarks geneset
l <- lapply(types, run_type_with_hallmarks, cnt_test, coldata_test, rnaseqdir=here('data/processed_data/rnaseq/truncating/'),rerun=T)

## load/extract DESeq results
rds_files <- dir(here('data/processed_data/rnaseq/truncating'),full.names=T)
rds_files <- grep('rds',rds_files,value=T)
l <- lapply(rds_files, load_data_results)
l <- rbindlist(l)
write.tsv(l,here('data/processed_data/rnaseq/truncating/deseq_cancertype_results.txt'))

## load/extract gsea results
rds_files <- dir(here('data/processed_data/rnaseq/truncating'),full.names=T)
rds_files <- grep('rds',rds_files,value=T)
l <- lapply(rds_files, load_gsea_results) 
l <- rbindlist(l)
write.tsv(l,here('data/processed_data/rnaseq_truncating_gsea_hallmark.txt'))




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# repeat for VUSs (protein-coding nontrunc, rRNA, tRNA)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## - load matrix with genotypes and wild-type status for all tcga samples
info <- copy(sample_hotspot_status_data)
info$PATIENT_ID <- strtrim(info$Tumor_Sample_Barcode,12)
clin <- copy(clin_data)
info <- merge(clin, info, by='PATIENT_ID', all.x=T)
info <- info[rna_available==T,]
info$status <- as.character(NA)
info[annotation_postfilter=='Wildtype',status:='wildtype']
info[annotation_postfilter=='VUS',status:='affected']
info <- info[!is.na(status),]

coldata_test <- info[status %in% c('wildtype','affected'),]
coldata_test <- coldata_test[,c('PATIENT_ID','maintype','status'),with=F]
coldata_test$status <- factor(coldata_test$status, levels=c('wildtype','affected'))

## subset cnt and coldata_test for samples in both
cnt_test <- copy(cnt_data)
valid_samples <- sort(intersect(colnames(cnt_test),coldata_test$PATIENT_ID))
cnt_test <- cnt_test[,valid_samples]

coldata_test <- coldata_test[PATIENT_ID %in% valid_samples,]
coldata_test <- as.data.frame(coldata_test) 
rownames(coldata_test) <- coldata_test$PATIENT_ID
coldata_test$PATIENT_ID <- NULL

## get valid types for testing (require 5+ wt and affected per cancer type)
tbl <- xtabs(~maintype+status ,data=coldata_test)
tbl <- cbind(type=rownames(tbl),adt(as.data.frame.matrix(tbl)))
names(tbl) <- c('type','wt','affected')
tbl <- tbl[affected >= 5 & wt >= 5]
tbl$N <- tbl$wt + tbl$affected 
tbl <- tbl[order(N,decreasing=F),]
types <- tbl$type
write.tsv(tbl,here('data/processed_data/rnaseq_vus_histogram.txt'))

## test each cancer type with hallmarks geneset
l <- lapply(types, run_type_with_hallmarks, cnt_test, coldata_test, rnaseqdir=here('data/processed_data/rnaseq/vus/'),rerun=T)

## load/extract DESeq results
rds_files <- dir(here('data/processed_data/rnaseq/vus'),full.names=T)
rds_files <- grep('rds',rds_files,value=T)
l <- lapply(rds_files, load_data_results)
l <- rbindlist(l)
write.tsv(l,here('data/processed_data/rnaseq/vus/deseq_cancertype_results.txt'))

## load/extract gsea results
rds_files <- dir(here('data/processed_data/rnaseq/vus'),full.names=T)
rds_files <- grep('rds',rds_files,value=T)
l <- lapply(rds_files, load_gsea_results) 
l <- rbindlist(l)
write.tsv(l,here('data/processed_data/rnaseq_vus_gsea_hallmark.txt'))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test R25Q vs WT in colorectal cancer
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cancer.type <- 'Colorectal Cancer'

## - load matrix with genotypes and wild-type status for all tcga samples
info <- copy(sample_hotspot_status_data)
info$PATIENT_ID <- strtrim(info$Tumor_Sample_Barcode,12)
clin <- copy(clin_data)
info <- merge(clin, info, by='PATIENT_ID', all.x=T)
info <- info[rna_available==T,]
info <- info[maintype==cancer.type,]
info$status <- as.character(NA)
info[annotation_postfilter=='Wildtype',status:='wildtype']
info[truncating==F & trna==F & rrna==F & `SNV MT-ND1 p.R25`==1,status:='affected']
info <- info[!is.na(status),]

coldata_test <- info[status %in% c('wildtype','affected'),]
coldata_test <- coldata_test[,c('PATIENT_ID','maintype','status'),with=F]
coldata_test$status <- factor(coldata_test$status, levels=c('wildtype','affected'))

## subset cnt and coldata_test for samples in both
cnt_test <- copy(cnt_data)
valid_samples <- sort(intersect(colnames(cnt_test),coldata_test$PATIENT_ID))
cnt_test <- cnt_test[,valid_samples]

coldata_test <- coldata_test[PATIENT_ID %in% valid_samples,]
coldata_test <- as.data.frame(coldata_test) 
rownames(coldata_test) <- coldata_test$PATIENT_ID
coldata_test$PATIENT_ID <- NULL

## test each cancer type with hallmarks geneset
l <- run_type_with_hallmarks(cancer.type, cnt_test, coldata_test, rnaseqdir=here('data/processed_data/rnaseq/r25q/'),rerun=T)

## load/extract DESeq results
rds_files <- dir(here('data/processed_data/rnaseq/r25q'),full.names=T)
rds_files <- grep('rds',rds_files,value=T)
l <- lapply(rds_files, load_data_results)
l <- rbindlist(l)
write.tsv(l,here('data/processed_data/rnaseq/r25q/deseq_cancertype_results.txt'))

## load/extract gsea results
rds_files <- dir(here('data/processed_data/rnaseq/r25q'),full.names=T)
rds_files <- grep('rds',rds_files,value=T)
l <- lapply(rds_files, load_gsea_results)
l <- rbindlist(l)
write.tsv(l,here('data/processed_data/rnaseq_r25q_gsea_hallmark.txt'))


