## find any positions in aligned trnas with positive selection based exact permutation tests

source(here::here('r/prerequisites.R'))


test_trna_hotspots <- function(d, output_file, samples_callable_per_position_data) {
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # for all tRNA positions, get their cloverleaf 
    # regions/positions, & total samples covered and mutated
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    d <- d[class=='somatic' & Variant_Type=='SNP' & Variant_Classification=='tRNA',]
    if(all(c('genbank_freq','helix_freq') %in% names(d))) d <- d[genbank_freq < 0.1 & helix_freq < 0.1,] 

    regions <- fread(chrM_annotations_file)
    strands <- regions[!duplicated(pos),c('symbol','strand'),with=F]
    strands <- strands[!duplicated(symbol),]

    al <- fread(alignment_file)
    al <- al[Position!='',]

    position_order <- unique(al$Position)
    al <- melt(al, id.vars=c('Position','index','Region'))
    setnames(al,'variable','Hugo_Symbol')
    setnames(al,'value','Start_Position')
    al$Hugo_Symbol <- gsub('TRN','MT-T',al$Hugo_Symbol)
    al$Position <- factor(al$Position, levels=position_order)

    ## resolve any positions not included in study
    m <- fread(chrM_annotations_file)
    m$id <- paste(m$symbol,m$pos)
    al$id <- paste(al$Hugo_Symbol,al$Start_Position)
    al <- al[id %in% m$id,]

    ## add coverage per position
    samples_callable_per_position_data[,callable_t:=NULL]
    setnames(samples_callable_per_position_data,'callable_tn','N')
    al <- merge(al, samples_callable_per_position_data, by.x='Start_Position', by.y='pos', all.x=T)
    al <- al[!is.na(Start_Position) & !is.na(N),]

    count <- function(d) {
        N <- length(unique(d$Tumor_Sample_Barcode))
        list(N=N)
    }
    cnt <- d[,count(.SD),by=Start_Position]

    al <- merge(al, cnt, by='Start_Position', all.x=T)
    setnames(al,'N.x','samples_covered')
    setnames(al,'N.y','samples_mutated')
    al[is.na(samples_mutated),samples_mutated:=0]
    al <- merge(al, strands, by.x='Hugo_Symbol', by.y='symbol',all.x=T)
    al <- al[order(Start_Position,decreasing=F),]

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # test for recurrent mutations at aligned position
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    collapse_position <- function(al) {   
        total_bps_callable_at_position <- sum(al$samples_covered,na.rm=T) 
        total_mutations_at_position <- sum(al$samples_mutated,na.rm=T)
        list(total_bps_callable_at_position=total_bps_callable_at_position, 
             total_mutations_at_position=total_mutations_at_position)
    }
    pos <- al[,collapse_position(.SD),by=Position]
    total_bps_callable_all_positions <- sum(pos$total_bps_callable_at_position)
    total_mutations_all_positions <- sum(pos$total_mutations_at_position)
    pos$prob_expected <- pos$total_bps_callable_at_position/total_bps_callable_all_positions

    test_position <- function(i,pos) {
        tst <- binom.test(pos$total_mutations_at_position[i],total_mutations_all_positions,
                          p=pos$prob_expected[i],alternative='greater')
        out <- pos[i,]
        out$p.value=tst$p.value
        out
    } 
    result <- rbindlist(lapply(1:nrow(pos), test_position, pos))
    result <- result[order(p.value,decreasing=F),]
    result$q.value <- p.adjust(result$p.value,method='BH') 
    result$nlog10q <- -log10(result$q.value)



    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # format/write hotspots output
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    result2 <- merge(al, result, by='Position', all.x=T)
    result2 <- result2[Position!='' & !is.na(Position),]
    result2[,id:=NULL]

    names(result2) <- c('trna_pos','Hugo_Symbol','Start_Position','index','trna_region',
                        'samples_covered',
                        'samples_mutated','strand',
                        'trna_bps_callable_at_aligned_pos',
                        'trna_samples_mutated_at_aligned_pos',
                        'prob_expected','p.value','q.value','nlog10q')


    result2$Hugo_Symbol <- gsub('TRN','MT-T',as.character(result2$Hugo_Symbol))
    result2 <- result2[order(Start_Position),]
    result2 <- merge(result2, trna_labels, by='Hugo_Symbol', all.x=T)
    result2[,aa:=NULL]
    setnames(result2,'label','AA')
    out <- result2

    ## add mutation counts per gene/pos
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

    count_samples <- function(d) {
        N <- length(unique(d$Tumor_Sample_Barcode))
        tumortypes <- collapse_info(d$maintype)
        ref <- d$Reference_Allele
        alt <- d$Tumor_Seq_Allele2
        pos <- d$pos
        dna_change <- collapse_info(paste0(ref,pos,alt))
        list(N=N,tumortypes=tumortypes,dna_change=dna_change)
    }
    tmp <- d[class=='somatic' & Variant_Classification=='tRNA' & Variant_Type=='SNP',
             c('Tumor_Sample_Barcode','maintype','Start_Position','Variant_Type',
               'Reference_Allele','Tumor_Seq_Allele2'),with=F]
    tmp$pos <- tmp$Start_Position
    tbl <- tmp[,count_samples(.SD),by=Start_Position]

    result2 <- merge(result2, tbl, by='Start_Position', all.x=T)
    result2[is.na(N),N:=0]
    result2[is.na(tumortypes),tumortypes:='']
    result2$hotspot_class <- 'tRNA'
    result2$End_Position <- result2$Start_Position
    result2$label <- paste(result2$Hugo_Symbol, result2$Start_Position, 'tRNA')
    result2 <- result2[order(Start_Position,Hugo_Symbol),]
    result2[,N:=NULL]
    result2[,End_Position:=NULL]
    setnames(result2,'Start_Position','pos')
    write.tsv(result2, output_file)

    out
}


alignment_file <- here('data/original_data/resources/trna_alignment.txt')
chrM_annotations_file <- here('data/original_data/resources/chrM_annotated.txt')
samples_callable_per_position_file <- here('data/processed_data/misc/samples_callable_per_position_precalculated.txt')

## run for TCGA
maf_file <- here('data/original_data/data_mutations_tcga.txt.gz')
output_file <- here('data/processed_data/hotspots/hotspots_trna.txt')

maf_data <- fread(maf_file)
samples_callable_per_position_data <- fread(samples_callable_per_position_file)
valid_positions <- samples_callable_per_position_data$pos
result_tcga <- test_trna_hotspots(maf_data, output_file, samples_callable_per_position_data)



## run for PCAWG for validation
maf_file <- here('data/original_data/data_mutations_pcawg.txt.gz')
output_file <- here('data/processed_data/hotspots/hotspots_trna_pcawg.txt')
clin_file <- here('data/original_data/data_samples_pcawg.txt.gz')

maf_data <- fread(maf_file)
clin_data <- fread(clin_file)
clin_data <- clin_data[grepl('TCGA',Tumor_Sample_Barcode)==F,]
setnames(clin_data,'dcc_project_code','maintype')
maf_data <- maf_data[Tumor_Sample_Barcode %in% clin_data$Tumor_Sample_Barcode,]
maf_data$class <- 'somatic'
maf_data <- merge(maf_data, clin_data[,c('Tumor_Sample_Barcode','maintype'),with=F], by='Tumor_Sample_Barcode', all.x=T)
n_pcawg <- nrow(clin_data)
samples_callable_per_position_data <- data.table(pos=samples_callable_per_position_data$pos, callable_tn=n_pcawg, callable_t=n_pcawg)
result_pcawg <- test_trna_hotspots(maf_data, output_file, samples_callable_per_position_data)





