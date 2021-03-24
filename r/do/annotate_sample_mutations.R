## for every sample, annotate which hotspots it has if any, and any other truncating, trna, rrna, snv SNV mutations

source(here::here('r/prerequisites.R'))

## input files:
snv_hotspots_file <- here('data/processed_data/hotspots_snv_precalculated.txt')
trna_hotspots_file <- here('data/processed_data/hotspots_trna_precalculated.txt')
indel_hotspots_file <- here('data/processed_data/hotspots_indel_precalculated.txt')

tcga_maf_file <- here('data/original_data/data_mutations_tcga.txt.gz')
tcga_samples_file <- here('data/original_data/data_samples_tcga.txt.gz')
tcga_tumornormal_coverage_file <- here('data/original_data/coverage_matrix_tcga.txt.gz')
tcga_tumoronly_coverage_file <- here('data/original_data/coverage_matrix_tcga_tumor_only.txt.gz')

impact_maf_file <- here('data/original_data/data_mutations_impact.txt.gz')
impact_samples_file <- here('data/original_data/data_samples_impact.txt.gz')
impact_tumornormal_coverage_file <- here('data/original_data/coverage_matrix_impact.txt.gz')
impact_tumoronly_coverage_file <- here('data/original_data/coverage_matrix_impact_tumor_only.txt.gz')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load results from 3 hotspots tests and get positions of their loci
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hotspots_regions_map_file <- here('data/processed_data/hotspot_region_map.txt')

if(!file.exists(hotspots_regions_map_file)) {
    message('Generating hotspot regions map ...')
    ## SNP hotspots
    snv_hotspots <- fread(snv_hotspots_file)
    snv_hotspots <- snv_hotspots[q.value < 0.01,c('label','Start_Position'),with=F]
    snv_hotspots <- snv_hotspots[order(Start_Position),]
    snv_hotspots$type <- 'SNV'

    ## tRNA hotspots
    trna_hotspots <- fread(trna_hotspots_file)
    trna_hotspots$label <- paste0(trna_hotspots$Hugo_Symbol,' pos',trna_hotspots$trna_pos,'.',trna_hotspots$pos)
    trna_hotspots <- trna_hotspots[q.value < 0.01,c('label','pos'),with=F]
    setnames(trna_hotspots,'pos','Start_Position')
    trna_hotspots <- trna_hotspots[Start_Position %nin% snv_hotspots$Start_Position]
    trna_hotspots$type <- 'tRNA'

    ## INDEL hotspots
    indel_hotspots <- fread(indel_hotspots_file)
    indel_hotspots <- indel_hotspots[q.value < 0.01]
    indel_hotspots$homopolymer <- paste(indel_hotspots$Hugo_Symbol,indel_hotspots$homopolymer)
    indel_hotspots <- indel_hotspots[,c('homopolymer','Start_Position','End_Position'),with=F]
    get_positions <- function(map) {
        ## extract the positions for each hotspot loci
        pos <- map$Start_Position:map$End_Position
        data.table(Start_Position=pos)
    }
    indel_hotspot_positions <- indel_hotspots[,get_positions(.SD),by=homopolymer]
    setnames(indel_hotspot_positions,'homopolymer','label')
    indel_hotspot_positions$type <- 'indel'

    ## merge togethe 
    hotspot_positions <- rbind(snv_hotspots, trna_hotspots, indel_hotspot_positions)
    valid_positions <- hotspot_positions$Start_Position
    write.tsv(hotspot_positions, hotspots_regions_map_file)

} else {
    message('Loading existing hotspot regions map ...')
    hotspot_positions <- fread(hotspots_regions_map_file)
    valid_positions <- hotspot_positions$Start_Position
}



get_sample_mtdna_status <- function(maf_file, samples_file, coverage_tumoronly_file, coverage_tumornormal_file, output_file, min_coverage, minvaf=NA, maxvaf=NA, somaticonly=F, tumoronly=F) { 

    hotspot_positions$label <- paste(hotspot_positions$type, hotspot_positions$label) 
    data_coverage_tumornormal <- fread(coverage_tumornormal_file,header=T)
    data_coverage_tumoronly <- fread(coverage_tumoronly_file,header=T)
    data_samples <- fread(samples_file)
    d_data <- fread(maf_file)

    d_data$censor <- F
    if(!is.na(minvaf)) d_data[TumorVAF < minvaf,censor:=T]
    if(!is.na(maxvaf)) d_data[TumorVAF >= maxvaf,censor:=T]
    if(somaticonly==T) d_data[class!='somatic',censor:=T]

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # for all samples, get the coverage at each hotspot loci, 0=not covered, 1=covered
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    data_coverage_tumoronly <- data_coverage_tumoronly[,c('Tumor_Sample_Barcode', as.character(valid_positions)),with=F]
    tonly_cov <- melt(data_coverage_tumoronly,id.var='Tumor_Sample_Barcode')
    data_coverage_tumornormal <- data_coverage_tumornormal[,c('Tumor_Sample_Barcode', as.character(valid_positions)),with=F]
    tnormal_cov <- melt(data_coverage_tumornormal,id.var='Tumor_Sample_Barcode')
    coverage <- merge(tonly_cov, tnormal_cov, by=c('Tumor_Sample_Barcode','variable'), all=T)
    coverage$id <- paste(coverage$Tumor_Sample_Barcode,coverage$variable)
    coverage <- coverage[!duplicated(id),]
    coverage[,id:=NULL]     

    names(coverage) <- c('Tumor_Sample_Barcode','pos','covered_tumor_only','covered_tumor_normal')   
    coverage[is.na(coverage)] <- 0   
    coverage$pos <- as.integer(as.character(coverage$pos))
    coverage <- merge(coverage, hotspot_positions, by.x='pos', by.y='Start_Position', all.x=T)
    coverage$covered <- coverage$covered_tumor_normal
    if(tumoronly==F) {
        coverage$covered <- coverage$covered_tumor_normal
    } else {
        coverage$covered <- coverage$covered_tumor_only
    }
    coverage[type=='indel',covered:=covered_tumor_only]
    coverage <- coverage[order(Tumor_Sample_Barcode,pos),]

    collapse_sample_hotspot_coverage <- function(coverage) {
        ## collapse coverage within hotspot positions, either entirely covered, or not covered
        covered <- all(coverage$covered==1)
        list(covered=covered)
    }
    coverage <- coverage[,collapse_sample_hotspot_coverage(.SD),by=c('Tumor_Sample_Barcode','label')]
    coverage$covered <- as.integer(coverage$covered)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # annotate the coverage data with if the mutation was actually observed
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    dat <- merge(d_data[,c('Tumor_Sample_Barcode','Start_Position','Variant_Classification','Variant_Type','TumorVAF','ShortVariantID','censor'), with=F], hotspot_positions,by='Start_Position', all.x=T)
    dat <- dat[(Variant_Classification %in% c('Frame_Shift_Ins','Frame_Shift_Del') & type=='indel') | (Variant_Type=='SNP' & type!='indel' & !is.na(type))]
    dat <- dat[order(Tumor_Sample_Barcode,Start_Position),]
    dat <- dat[,c('Tumor_Sample_Barcode','label','type','TumorVAF','ShortVariantID','censor'),with=F]
    dat$id <- paste(dat$Tumor_Sample_Barcode, dat$label)
    dat <- dat[order(id,TumorVAF,decreasing=T),]
    dat <- dat[!duplicated(id),]
    dat[,ShortVariantID:=NULL]

    coverage <- merge(coverage, dat, by=c('Tumor_Sample_Barcode','label'), all.x=T)
    coverage$affected <- as.integer(!is.na(coverage$type))
    coverage[covered==0,affected:=NA]

    coverage_matrix <- dcast(coverage, Tumor_Sample_Barcode ~ label, value.var='affected')
    coverage_matrix$N_hotspots <- rowSums(coverage_matrix[,(hotspot_positions$label),with=F])
    coverage_matrix$N_hotspots_indel <- rowSums(coverage_matrix[,(unique(hotspot_positions$label[hotspot_positions$type=='indel'])), with=F])
    coverage_matrix$N_hotspots_snv <- rowSums(coverage_matrix[,(unique(hotspot_positions$label[hotspot_positions$type=='SNV'])), with=F])
    coverage_matrix$N_hotspots_trna <- rowSums(coverage_matrix[,(unique(hotspot_positions$label[hotspot_positions$type=='tRNA'])), with=F])

    fields <- c('Tumor_Sample_Barcode','N_hotspots','N_hotspots_indel','N_hotspots_snv','N_hotspots_trna')
    tmp <- hotspot_positions[!duplicated(label),]
    tmp$type <- factor(tmp$type, levels=c('indel','SNV','tRNA'))
    tmp <- tmp[order(type,Start_Position),]
    fields <- c(fields, tmp$label)
    coverage_matrix <- coverage_matrix[,(fields),with=F]


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # annotate clinical data with mt-status and variants
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    dt <- d_data
    summarize_sample <- function(tmp) {
        censor <- any(tmp$Variant_Classification %in% c('Frame_Shift_Ins','Frame_Shift_Del','Nonsense_Mutation') & tmp$censor==T)
        truncating <- any(tmp$Variant_Classification %in% c('Frame_Shift_Ins','Frame_Shift_Del','Nonsense_Mutation'))
        rrna <- any(tmp$Variant_Classification %in% c('rRNA'))
        trna <- any(tmp$Variant_Classification %in% c('tRNA'))
        missense <- any(tmp$Variant_Classification %in% c('Missense_Mutation','In_Frame_Del','In_Frame_Ins','Translation_Start_Site','Nonstop_Mutation'))
        silent <- any(tmp$Variant_Classification %in% 'Silent')
        list(silent=silent,missense=missense,trna=trna,rrna=rrna,truncating=truncating, censor=censor)
    }
    l2 <- dt[,summarize_sample(.SD),by=Tumor_Sample_Barcode]
    clin <- data_samples
    if('CANCER_TYPE' %nin% names(clin)) clin$CANCER_TYPE <- clin$maintype

    if(tumoronly==F) {
        clin$covered <- clin$coverage_tumornormal >= min_coverage
    } else {
        clin$covered <- clin$coverage_tumor >= min_coverage
    }

    l2 <- merge(l2, clin[,c('CANCER_TYPE','maintype','Tumor_Sample_Barcode','covered','coverage_tumornormal','coverage_tumor'),with=F], by='Tumor_Sample_Barcode',all.y=T)
    l2 <- merge(l2, coverage_matrix, by='Tumor_Sample_Barcode', all.x=T)
    l2[is.na(silent),silent:=F]
    l2[is.na(missense),missense:=F]
    l2[is.na(trna),trna:=F]
    l2[is.na(rrna),rrna:=F]
    l2[is.na(truncating),truncating:=F]
    l2[is.na(censor),censor:=F]

    ## prefilter is for the barplot showing incidence of variants across cancer types
    ## we only consider samples with adequate coverage. Also this breakdown includes Silent
    l2$annotation_detailed_prefilter <- as.character(NA)
    l2[is.na(annotation_detailed_prefilter) & (covered==F | is.na(N_hotspots_indel)), annotation_detailed_prefilter:='Unknown']
    l2[is.na(annotation_detailed_prefilter) & censor==T, annotation_detailed_prefilter:='Censor']
    l2[is.na(annotation_detailed_prefilter) & truncating==T, annotation_detailed_prefilter:='Truncating']
    l2[is.na(annotation_detailed_prefilter) & (trna + rrna + missense) > 1, annotation_detailed_prefilter:='2+ types, non-truncating/silent']
    l2[is.na(annotation_detailed_prefilter) & trna==T, annotation_detailed_prefilter:='tRNA']
    l2[is.na(annotation_detailed_prefilter) & rrna==T, annotation_detailed_prefilter:='rRNA']
    l2[is.na(annotation_detailed_prefilter) & missense==T, annotation_detailed_prefilter:='Missense']
    l2[is.na(annotation_detailed_prefilter) & silent==T, annotation_detailed_prefilter:='Silent']
    l2[is.na(annotation_detailed_prefilter), annotation_detailed_prefilter:='Wildtype']

    ## postfilter is for the RNA-Seq and Survival data. Samples with ANY truncating variants
    ## are called Truncating. Then samples with inadequate coverage or missing coverage at indel
    ## hoptspot loci are removed as Unknown. Silent is combined with Wild-type
    l2$annotation_detailed_postfilter <- as.character(NA)
    l2[is.na(annotation_detailed_postfilter) & censor==T, annotation_detailed_postfilter:='Censor']
    l2[is.na(annotation_detailed_postfilter) & truncating==T, annotation_detailed_postfilter:='Truncating']
    l2[is.na(annotation_detailed_postfilter) & (covered==F | is.na(N_hotspots_indel)), annotation_detailed_postfilter:='Unknown']
    l2[is.na(annotation_detailed_postfilter) & (trna + rrna + missense) > 1, annotation_detailed_postfilter:='2+ types, non-truncating/silent']
    l2[is.na(annotation_detailed_postfilter) & trna==T, annotation_detailed_postfilter:='tRNA']
    l2[is.na(annotation_detailed_postfilter) & rrna==T, annotation_detailed_postfilter:='rRNA']
    l2[is.na(annotation_detailed_postfilter) & missense==T, annotation_detailed_postfilter:='Missense']
    l2[is.na(annotation_detailed_postfilter), annotation_detailed_postfilter:='Wildtype']

    ## add broader annotation groupings: Truncating, SNV, WT for both pre/post-filters
    l2$annotation_prefilter <- l2$annotation_detailed_prefilter
    l2[annotation_prefilter %in% c('Missense','tRNA','rRNA','2+ types, non-truncating/silent'), annotation_prefilter:='VUS']
    l2[annotation_prefilter %in% c('Silent','Wildtype'), annotation_prefilter:='Wildtype']
    l2$annotation_postfilter <- l2$annotation_detailed_postfilter
    l2[annotation_postfilter %in% c('Missense','tRNA','rRNA','2+ types, non-truncating/silent'), annotation_postfilter:='VUS']
    l2[annotation_postfilter %in% c('Silent','Wildtype'), annotation_postfilter:='Wildtype']
    
    message('Writing result: ',output_file)
    write.tsv(l2,output_file)
    l2
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate sample-level mtDNA-status classification
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Run for TCGA
qc <- get_sample_mtdna_status(tcga_maf_file, tcga_samples_file, tcga_tumoronly_coverage_file, tcga_tumornormal_coverage_file, here('data/processed_data/sample_hotspot_status_tcga.txt'), min_coverage=0.9, tumoronly=F, somaticonly=F)

## Run for IMPACT. We will use likely-somatic, somatic and truncating variants:
## - likely somatic are those which are heteroplasmic and never seen impact, tcga or helix germline
## - we also require only 60% coverage in tumor
qc <- get_sample_mtdna_status(impact_maf_file, impact_samples_file, impact_tumoronly_coverage_file, impact_tumornormal_coverage_file, here('data/processed_data/sample_hotspot_status_impact.txt'), min_coverage=0.60, tumoronly=T, somaticonly=F)


