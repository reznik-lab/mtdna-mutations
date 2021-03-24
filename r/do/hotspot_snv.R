## implementation of mt-hotspots algorithm 

source(here::here('r/prerequisites.R'))

## input files
mutation_file <- here('data/original_data/data_mutations_tcga.txt.gz')
clinical_file <- here('data/original_data/data_samples_tcga.txt.gz')
fa_file <- here('data/original_data/resources/Homo_sapiens.GRCh38v95.dna.chromosome.MT.fa')
chrM_annotated_file <- here('data/original_data/resources/chrM_annotated.txt')
covg_per_pos_file <- here('data/processed_data/misc/samples_callable_per_position_precalculated.txt')

## output files
hotspots_output <- here('data/processed_data/hotspots/hotspots_snv.txt') ## results 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# tabulate the 64 unique possible trinucleotides for all retained positions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fa <- fread(fa_file)[[1]]
fa <- paste(fa,collapse='')
dat <- data.table(pos=1:16569)
dat$nt <- strsplit(fa,'')[[1]]
dat$ntMinus1 <- c('-',dat$nt[1:(nrow(dat)-1)])
dat$ntPlus1 <- c(dat$nt[2:(nrow(dat))],'-')
dat$trinuc <- paste0(dat$ntMinus1,dat$nt,dat$ntPlus1)
dat <- dat[,c('pos','trinuc'),with=F]

d <- fread(chrM_annotated_file)
d <- merge(d, dat, by='pos', all.x=T)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate trinucleotide mutabilities
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

maf <- fread(mutation_file)
maf <- maf[Variant_Type=='SNP' & class=='somatic'] # & Variant_Classification!='Silent']
maf <- maf[helix_freq < 0.1 & genbank_freq < 0.1,]
tbl <- table.freq(maf$flanking_bps)
setnames(tbl,'N','s_trinuc')
tbl$s_total <- sum(tbl$s_trinuc)
setnames(tbl,'value','trinuc')
d <- merge(d, tbl, by='trinuc', all.x=T)
d[is.na(s_trinuc),s_trinuc:=0]




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add number of samples with coverage at position and overall
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

covg_per_pos <- fread(covg_per_pos_file)
covg_per_pos[,callable_t:=NULL]
setnames(covg_per_pos,'callable_tn','covd_samples')
d <- merge(d, covg_per_pos, by='pos', all.x=T)
samples <- fread(clinical_file,select='Tumor_Sample_Barcode')[[1]]
n_samples <- length(samples)
d$total_samples <- n_samples


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add number of samples mutated per gene
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## get mutations per position
mutation_counts <- table.freq(maf$Start_Position)
mutation_counts$value <- as.integer(mutation_counts$value)
d <- merge(d, mutation_counts, by.x='pos', by.y='value', all.x=T)
d[is.na(N),N:=0]
setnames(d,'N','mutations_at_pos')


## add mutations per gene
tbl <- table.freq(maf$Hugo_Symbol)
setnames(tbl,'N','mutations_in_gene')
d <- merge(d, tbl, by.x='symbol', by.y='value', all.x=T)
d[is.na(mutations_in_gene),mutations_in_gene:=0]
d$mu_pos <- (d$s_trinuc / d$s_total) * (d$covd_samples / d$total_samples)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# normalize the position mutability by the gene's total trinucleotide-mutability
# this is the binomial parameter
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get_binomial_parameter <- function(res) {
    res$binomial_probability <- res$mu_pos / sum(res$mu_pos)
    res
}
d <- d[,get_binomial_parameter(.SD),by=symbol]
d <- d[order(pos),]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test each position
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test_position <- function(dd) {
    successes <- dd$mutations_at_pos
    trials <- dd$mutations_in_gene
    probability <- dd$binomial_probability
    p.value <- binom.test(x=successes, n=trials, p=probability, alternative='two.sided')$p.value
    out <- list(p.value=p.value)
    out
}
hs <- d[, test_position(.SD), by=pos]
res <- merge(d, hs, by='pos', all=F)
res$q.value <- p.adjust(res$p.value,method='BH')
res <- res[order(p.value, decreasing=F),]



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add additional info about each hotspot position
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

toadd <- maf[,c('Tumor_Sample_Barcode','maintype','Start_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','HGVSp_Short'),with=F]
info <- HGVSp_Short_parse(toadd$HGVSp_Short)
info[,HGVSp_Short:=NULL]
toadd <- cbind(toadd, info)
toadd[Reference_Amino_Acid=='',Reference_Amino_Acid:=NA]
toadd$start <- toadd$Start_Position

collapse_info <- function(s) {
    s[is.na(s)] <- 'NA'
    tbl <- table.freq(s)
    include.n <- nrow(tbl) > 1
    if(include.n) {
        label <- paste(paste0(tbl$value,'(',tbl$N,')'), collapse='; ')
    } else {
        label <- paste(tbl$value,collapse='; ')
    }
    label
}

collapse_pos <- function(x) {
    n_samples <- length(unique(x$Tumor_Sample_Barcode))
    cancertypes <- collapse_info(x$maintype)
    ref <- x$Reference_Allele[1]
    alt <- x$Tumor_Seq_Allele2
    pos <- x$Start_Position[1]
    dna_change <- collapse_info(paste0(ref,pos,alt))
    variant_classification <- collapse_info(x$Variant_Classification)
    reference_amino_acid <- x$Reference_Amino_Acid[1]
    amino_acid_position <- x$Amino_Acid_Position[1]
    variant_amino_acid <- x$Variant_Amino_Acid
    protein_change <- collapse_info(paste0(reference_amino_acid,amino_acid_position,variant_amino_acid))

    if(!is.na(amino_acid_position)) {
        label <- paste0('p.',reference_amino_acid,amino_acid_position)
    } else {
        protein_change <- as.character(NA)
        label <- paste0('c.',pos)
    }
    
    list(n_samples=n_samples, cancertypes=cancertypes, 
         dna_change=dna_change,
         protien_change=protein_change, label=label)
}

position_info <- toadd[,collapse_pos(.SD),by=start]

out <- merge(res, position_info, by.x='pos', by.y='start', all.x=T)
out <- out[order(p.value,decreasing=F),]
out$label <- paste(out$symbol,out$label,sep=' ')
out <- out[mutations_at_pos >= 1,]
setnames(out,'pos','Start_Position')
setnames(out,'symbol','Hugo_Symbol')
write.tsv(out,hotspots_output)


