## implementation of mt-hotspots algorithm for indels

source(here::here('r/prerequisites.R'))

## input files
tumor_coverage_file <- here('data/original_data/coverage_matrix_tcga_tumor_only.txt.gz')
homopolymer_track_file <- here('data/original_data/resources/chrM_homopolymers.txt')
common_variants_file <- here('data/original_data/resources/mitomap_common_variants_20191015.txt')
data_mutations <- here('data/original_data/data_mutations_tcga.txt.gz')
chrM_annotated_file <- here('data/original_data/resources/chrM_annotated.txt')

## output files
binomial_test_results_file <- here('data/processed_data/hotspots/hotspots_indel.txt')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create table of sample coverage per unique homopolymer repeat locus (tumor-only coverage)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## load padded homopolymer regions 
homopolymer <- fread(homopolymer_track_file)
setnames(homopolymer,'homopolymer','id')
homopolymer$homopolymer <- paste0(homopolymer$start,':',homopolymer$end,',',homopolymer$id)
homopolymer$start <- homopolymer$start-1; homopolymer$end <- homopolymer$end+1 ## pad around the homopolymer by +/- 1 bp
homopolymer$Chromosome <- 'MT'
setnames(homopolymer,'start','homopolymer_start')
setnames(homopolymer,'end','homopolymer_end')
setkey(homopolymer,'Chromosome','homopolymer_start','homopolymer_end')
expand <- function(homopolymer) {
    pos <- seq(homopolymer$homopolymer_start,homopolymer$homopolymer_end)
    list(pos=pos)
}
homopolymer_expanded <- homopolymer[,expand(.SD),by=c('homopolymer','id','sequence')]


collapse <- function(d) {
    out <- d[1,]
    out$homopolymer_start <- min(d$pos)
    out$homopolymer_end <- max(d$pos)
    out$rep <- unique(d$id)
    out
}
collapsed_homopolymer <- homopolymer_expanded[,collapse(.SD),by=homopolymer]
collapsed_homopolymer[,c('pos','id'):=NULL]
collapsed_homopolymer$Chromosome <- 'MT'
setkey(collapsed_homopolymer,'Chromosome','homopolymer_start','homopolymer_end')
collapsed_homopolymer$homopolymer_width <- nchar(collapsed_homopolymer$sequence)
setnames(collapsed_homopolymer,'rep','homopolymer_element')


## for each sample in each homopolymer, get the number of basepairs where is it covered
regions <- fread(tumor_coverage_file,header=T) 
ids <- regions$Tumor_Sample_Barcode
regions[,Tumor_Sample_Barcode:=NULL]
regions <- t(regions)
regions <- cbind(pos=rownames(regions),adt(regions))
names(regions)[2:ncol(regions)] <- ids
regions$pos <- suppressWarnings(as.integer(regions$pos))

region_annotations <- fread(chrM_annotated_file)
region_annotations[,strand:=NULL]
setnames(region_annotations,'pos','start')
region_annotations$end <- region_annotations$start + 1
region_annotations$chr <- 'MT'
regions <- merge(region_annotations, regions, by.x='start', by.y='pos', all.x=T)


setkey(regions,'chr','start','end')
regions_annotated <- foverlaps(regions, collapsed_homopolymer, type='any')
regions_annotated <- regions_annotated[!is.na(homopolymer)]
regions_annotated <- regions_annotated[start >= homopolymer_start & start <= homopolymer_end]
sample_fields <- grep('TCGA-',names(regions),value=T)
regions_annotated <- regions_annotated[,c('homopolymer','homopolymer_start','homopolymer_end','homopolymer_width','homopolymer_element',sample_fields),with=F]
collapse_homopolymer <- function(d) {
    tmp <- d[,(sample_fields),with=F]
    nucleotides_covered <- colSums(tmp)
    out <- as.list(nucleotides_covered)
    out
}
coverage <- regions_annotated[,collapse_homopolymer(.SD),
                              by=c('homopolymer','homopolymer_start','homopolymer_end',
                                   'homopolymer_element','homopolymer_width')]
coverage_melted <- melt(coverage,id.vars=c('homopolymer','homopolymer_start','homopolymer_end','homopolymer_element','homopolymer_width'))
setnames(coverage_melted,'variable','Tumor_Sample_Barcode')
coverage_melted <- coverage_melted[value==homopolymer_width+2,] ## subset the sample/homopolymer combos for those where we had sufficient coverage
coverage_melted[,value:=NULL]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# for each homopolymer, get the number of sample covered and 
# number of samples with frame-shift indels
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## load indels, merge with homopolymers (note some samples/mutations may arise in several homopolymers at once)
d <- fread(data_mutations)
d <- d[Variant_Classification %in% c('Frame_Shift_Ins','Frame_Shift_Del'),] 
setkey(d,'Chromosome','Start_Position','End_Position')
d <- foverlaps(d, collapsed_homopolymer, type='any')
d <- d[!is.na(homopolymer)]
d <- d[Start_Position >= homopolymer_start & End_Position <= homopolymer_end]
d <- d[,c('Tumor_Sample_Barcode','maintype','Start_Position','Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','HGVSp_Short','homopolymer','ShortVariantID'),with=F]
d <- d[!is.na(homopolymer),]

## annotate the coverage per homopolymer with mutations
coverage_melted <- merge(coverage_melted, d, by=c('Tumor_Sample_Barcode','homopolymer'), all.x=T)
x <- coverage_melted[,c('Tumor_Sample_Barcode','ShortVariantID','homopolymer','homopolymer_element','homopolymer_width'),with=F]
x$affected <- !is.na(coverage_melted$ShortVariantID)
x$nuc_repeat <- x$homopolymer_element %in% c('A','C','T','G')


## prep data for various potential enrichment tests
get_number_affected_per_homopolymer <- function(x) {
    mutant <- length(unique(x$Tumor_Sample_Barcode[x$affected==T]))
    wildtype <- length(unique(x$Tumor_Sample_Barcode[x$affected==F]))
    list(mutant=mutant,wildtype=wildtype)
}
xx <- x[,get_number_affected_per_homopolymer(.SD),by=c('homopolymer_element','homopolymer')]
xx$covered <- xx$wildtype+xx$mutant 
xx$prop <- xx$mutant/xx$covered
xx <- merge(xx, homopolymer, by='homopolymer', all.x=T)
setnames(xx,'width','homopolymer_element_width')
xx$homopolymer_width <- nchar(xx$sequence)


## add gene annotations to data
genes <- fread(chrM_annotated_file)
collapse_gene <- function(genes) {
    gene_start <- min(genes$pos)
    gene_end <- max(genes$pos)
    list(gene_start=gene_start, gene_end=gene_end)
}
genes <- genes[,collapse_gene(.SD),by=c('symbol','strand')]
genes$chr <- 'MT'
setkey(genes,'chr','gene_start','gene_end')
setkey(xx,'Chromosome','homopolymer_start','homopolymer_end')
xx <- foverlaps(xx, genes, type='within')
xx <- xx[grepl('MT-R',symbol)==F & grepl('MT-T',symbol)==F,] ## exclude RNAs from this
Ns <- length(sample_fields)
xx$prop_callable <- xx$covered/Ns


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run binomial test
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tmp <- xx
total_homopolymer_indels <- sum(tmp$mutant)
total_homopolymer_regions_length <- sum(tmp$homopolymer_width) ## this is the padded width
tmp <- tmp[order(homopolymer_start)]

test <- function(tmp) {
    cnt_affected <- tmp$mutant
    homopolymer_length <- tmp$homopolymer_width
    probability <- homopolymer_length/total_homopolymer_regions_length
    p <- tryCatch({
        binom.test(cnt_affected,total_homopolymer_indels,probability,alternative='greater')$p.value
    },error=function(e) {
        as.numeric(NA)
    })
    list(x=cnt_affected,probability=probability,p=p)
}
l <- tmp[,test(.SD),by=homopolymer]
l$q <- p.adjust(l$p,method='BH')
l <- l[order(p,decreasing=F),]
tmp <- merge(tmp, l, by='homopolymer', all.x=T)


## format output. We will only include nuc_repeats.
setnames(tmp,'id','homopolymer_element')
out <- tmp[,c('homopolymer','homopolymer_start','homopolymer_end','homopolymer_element',
              'symbol','Chromosome','strand','gene_start','gene_end',
              'mutant','covered','prop','sequence',
              'homopolymer_width','prop_callable','p','q'),with=F]
setnames(out,'p','p.value')
setnames(out,'q','q.value')
setnames(out,'symbol','gene')
out <- out[order(p.value,decreasing=F)]
out$labels <- paste(out$gene,out$homopolymer) 
out$start <- out$homopolymer_start
out$end <- out$homopolymer_end
out[,c('homopolymer_start','homopolymer_end'):=NULL]
result <- out


## get data for observed mutations at each homopolymer
toadd <- coverage_melted[!is.na(ShortVariantID),c('Tumor_Sample_Barcode','maintype','Start_Position','Reference_Allele',
                             'Tumor_Seq_Allele2','Variant_Classification','HGVSp_Short','homopolymer'),with=F]
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
    ref <- x$Reference_Allele
    alt <- x$Tumor_Seq_Allele2
    pos <- x$Start_Position
    dna_change <- collapse_info(paste0(ref,pos,alt))
    variant_classification <- collapse_info(x$Variant_Classification)
    reference_amino_acid <- x$Reference_Amino_Acid
    amino_acid_position <- x$Amino_Acid_Position
    variant_amino_acid <- x$Variant_Amino_Acid
    protein_change <- collapse_info(paste0(reference_amino_acid,amino_acid_position,variant_amino_acid))

    list(n_samples=n_samples, cancertypes=cancertypes, 
         dna_change=dna_change,
         protien_change=protein_change)
         #Variant_Classification=variant_classification,
         #Reference_Amino_Acid=reference_amino_acid, 
         #Amino_Acid_Position=amino_acid_position, 
         #Variant_Amino_Acid=variant_amino_acid)
}

position_info <- toadd[,collapse_pos(.SD),by=homopolymer]
out <- merge(result, position_info, by='homopolymer', all.x=T)
out[is.na(n_samples),n_samples:=0]
setnames(out,'gene','Hugo_Symbol')
out <- out[order(homopolymer,Hugo_Symbol),]
out$hotspot_class <- 'homopolymer-indel'
setnames(out,'start','Start_Position')
setnames(out,'end','End_Position')
out[,c('n_samples'):=NULL]
out$label <- paste(out$Hugo_Symbol, out$homopolymer, 'homopolymer-indel')
out[q.value >= 0.01, hotspot_class:='']
out <- out[order(p.value,Start_Position,decreasing=F),]
out[,labels:=NULL]
write.tsv(out, binomial_test_results_file)




