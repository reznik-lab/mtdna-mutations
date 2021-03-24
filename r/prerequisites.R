## set seed for reproducibility of figures
set.seed(42)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load required packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

required.packages <- c('data.table','ggplot2','cowplot','RColorBrewer','parallel','ggsignif','binom','scales','ggrepel','Hmisc','here','magrittr','knitr','rmarkdown','epitools','survminer','pheatmap')
hide <- suppressMessages(lapply(required.packages, require, character.only = TRUE))
missing.packages <- required.packages[!required.packages %in% (.packages())]
if(length(missing.packages)>0) stop(paste('Could not load required packages:',paste(missing.packages,collapse=', ')))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define commonly used objects and helper functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

coding_classes <- c('Missense_Mutation','Nonsense_Mutation','Splice_Site','In_Frame_Del','In_Frame_Ins','Frame_Shift_Del','Frame_Shift_Ins','Translation_Start_Site','Nonstop_Mutation','Silent')
truncating_classes <- c('Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins','Splice_Site')

mtdna_classes <- unique(c(coding_classes,truncating_classes,'tRNA','rRNA'))
mtdna_classes <- mtdna_classes[mtdna_classes!='Silent']
nontruncating_classes <- mtdna_classes[mtdna_classes %nin% truncating_classes]

trna_labels <- data.table(aa=c('A','R','N','D','C','Q','E','G','H','I','L1','L2','K','M','F','P','S1','S2','T','W','Y','V'),
                          label=c('Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Ile','Leu1','Leu2','Lys','Met','Phe','Pro','Ser1','Ser2','Thr','Trp','Tyr','Val'))
trna_labels$Hugo_Symbol <- paste0('MT-T',trna_labels$aa)


## ggplot theme
theme_std <- function(base_size = 11, base_line_size = base_size/22, base_rect_size = base_size/22) {
    theme_classic(base_size = base_size, base_family = 'ArialMT')  %+replace%
    theme(
          line = element_line(colour = "black", size = base_line_size, linetype = 1, lineend = "round"),
          text = element_text(family = 'ArialMT', face = "plain", colour = "black", size = base_size, lineheight = 0.9, hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), debug=F),
          axis.text = element_text(colour = "black", family='ArialMT', size=rel(0.8)),
          axis.ticks = element_line(colour = "black", size=rel(1)),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black", size = rel(1)),
          legend.key = element_blank(),
          strip.background = element_blank())
}


## function to force axis break in ggplot
break_axis <- function(y, maxlower, minupper=NA, lowerticksize, upperticksize, ratio_lower_to_upper) {
    if(is.na(minupper)) {
        breakpos <- maxlower
        lowerticklabels <- seq(0,breakpos,by=lowerticksize); lowerticklabels
        upperticklabels <- seq(breakpos+upperticksize,max(y)+upperticksize,by=upperticksize); upperticklabels
        ticklabels <- c(lowerticklabels, upperticklabels); ticklabels
        lowertickpos <- lowerticklabels
        uppertickspacing <- ratio_lower_to_upper * lowerticksize
        uppertickpos <- breakpos + ((1:length(upperticklabels))*uppertickspacing)
        tickpos <- c(lowertickpos, uppertickpos)
        newy <- as.numeric(y)
        ind <- newy > breakpos
        newy[ind] <- breakpos + uppertickspacing*((newy[ind]-breakpos) / upperticksize)
        list(newy=newy, breaks=tickpos, labels=ticklabels, limits=range(tickpos))
    } else {
        lowerticklabels <- seq(0,maxlower,by=lowerticksize); lowerticklabels
        upperticklabels <- seq(minupper,max(y)+upperticksize,by=upperticksize); upperticklabels
        ticklabels <- c(lowerticklabels, upperticklabels); ticklabels
        lowertickpos <- lowerticklabels
        uppertickspacing <- ratio_lower_to_upper * lowerticksize
        uppertickpos <- maxlower + 0.5*lowerticksize + ((1:length(upperticklabels))*uppertickspacing)
        tickpos <- c(lowertickpos, uppertickpos)
        newy <- as.numeric(y)
        ind <- newy > maxlower
        newy[ind] <- maxlower + 0.5*lowerticksize + 1*uppertickspacing + uppertickspacing*((newy[ind]-minupper) / upperticksize)
        list(newy=newy, breaks=tickpos, labels=ticklabels, limits=range(tickpos))
    }
}


## function to convert ggplot object into separate objects for plot vs legend
extract_gglegend <- function(p){
    ## extract the legend from a ggplot object
    tmp <- ggplot_gtable(ggplot_build(p))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    if(length(leg) > 0) leg <- tmp$grobs[[leg]]
    else leg <- NULL
    leg

    ## return the legend as a ggplot object
    legend <- cowplot::ggdraw() + cowplot::draw_grob(grid::grobTree(leg))
    plot <- p + theme(legend.position='none')
    list(plot=plot,legend=legend)
}


## function to split HGVSp_Short field in MAF into Amino_Acid_Position, Reference_Amino_Acid, Variant_Amino_Acid, Amino_acids
HGVSp_Short_parse <- function(x,cpus=1) {
    x1 <- gsub('p.','',x)
    m <- gregexpr('[0-9]+',x1)
    nums <- regmatches(x1,m)
    get.aa <- function(num) num[1]
    aas <- sapply(nums,get.aa,USE.NAMES=F)
    tmp <- data.table(HGVSp_Short=x1,aa=aas)
    tmp$i <- 1:nrow(tmp)
    split <- function(d) {
        s <- strsplit(d$HGVSp_Short,d$aa)[[1]]
        s[2] <- gsub('_sice','splice',s[2])
        list(Reference_Amino_Acid=s[1],Variant_Amino_Acid=s[2])
    }
    info <- tmp[,split(.SD),by=i]
    info$HGVSp_Short <- x
    info$Amino_Acid_Position=as.integer(aas)
    info <- info[,c(4,2,5,3),with=F]
    info$Amino_acids <- paste(info$Reference_Amino_Acid,info$Variant_Amino_Acid,sep='/')
    info$Amino_acids[is.na(info$Amino_Acid_Position)] <- NA
    info
}


## creates a data.table with histogram of values in a vector 
table.freq <- function(value) {
    if(is.null(value) | length(value)==0) {
        tbl <- data.table(value=NA,N=NA)
    } else {
        tbl <- adt(table(value))
        tbl <- tbl[order(tbl$N,decreasing=T),]
    }
    tbl
}


## shortcut for sort(unique(...))
sortunique <- function(x,...) sort(unique(x),na.last=T,...)

## shortcut to write tab-delimited data with consistent format 
write.tsv <- function(d, file, sep = "\t", quote = F, row.names = F, ...)  write.table(d, file = file, sep = sep, quote = quote, row.names = row.names, ...)

## shortcut for as.data.table
adt <- function(d) as.data.table(d)

## standardize p-value labels
plabel <- function(p) {
    if(p==0) {
        p <- .Machine$double.xmin
        p <- paste0('P<',prettyNum(p,digits=1))
    } else {
        p <- paste0('P=',prettyNum(p,digits=1))
    }
    p
}





