# mtdna-mutations

![alt text](https://github.com/reznik-lab/mtdna-mutations/blob/master/data/wiki.png "Hi there!")

# Instructions:
The following steps can be used to recreate the essential parts of the figures in the manuscript. All code is in R (>=3.6.0).

Clone this repo:
```shell
git clone https://github.com/reznik-lab/mtdna-mutations
cd mtdna-mutations
```

Install required R packages from CRAN:

```r
## check for missing required packages, install them.
required.packages <- c('data.table','ggplot2','cowplot','RColorBrewer',
                       'parallel','ggsignif','binom','scales',
                       'MASS','ggrepel','Hmisc','Rcpp',
                       'car','here','magrittr','knitr','rmarkdown')

new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)
```

Generate html file with figures in the main text and extended data:
```shell
R -e "rmarkdown::render('html/figures.Rmd', output_file = 'figures.html')"
```

## Optional: 

The above commands use data provided in this repo, including a number of "precalculated" datasets based on methods described in the manuscript. All such precalculated datasets are in the data/processed_data subdirectory, and can be recreated with the following commands:

```r
## If Bioconductor isn't already installed, install it.
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

## Install additional required packages from Bioconductor if not installed.
if(!'DESeq2' %in% installed.packages()) BiocManager::install("DESeq2")
if(!'apeglm' %in% installed.packages()) BiocManager::install("apeglm")
if(!'fgsea' %in% installed.packages()) BiocManager::install("fgsea")
```

1. Generate misc tables used in other scripts based on sample coverage data: 
```shell
## Generate a table of number of samples covered per position, 
## and a matrix of effective-gene-lengths for each sample. Generates: 
## - data/processed_data/misc/samples_callable_per_position.txt
## - data/processed_data/misc/sample_gene_lengths.txt.gz
Rscript r/do/process_sample_coverage.R
```

2. Determine recurrent positions for mtDNA SNV, truncating indels, and tRNA mutations: 
```shell
## Calculate SNP hotspots. Generates: 
## - data/processed_data/hotspots/hotspots_snp.txt
Rscript r/do/hotspot_snp.R

## Calculate frame-shift indel hotspots. Generates:
## - data/processed_data/hotspots/hotspots_indel.txt
Rscript r/do/hotspot_indel.R

## Calculate tRNA alignment hotspots. Generates:
## - data/processed_data/hotspots/hotspots_trna.txt
Rscript r/do/hotspot_trna.R
```

3. Annotate samples with mtDNA variants, % callable, and overall mtDNA-status:
```shell
## Generates: 
## - data/processed_data/hotspots/sample_hotspot_status.txt 
## - data/processed_data/hotspots/hotspot_region_map.txt
Rscript r/do/annotate_sample_mutations.R
```

4. Calculate differentially-expressed genes and run GSEA:
```shell
## Generates: 
## - data/processed_data/rnaseq/
## Download TCGA RNA-Seq RSEM counts (1.8GB file, will take some time):
wget http://api.gdc.cancer.gov/data/3586c0da-64d0-4b74-a449-5ff4d9136611 -O data/original_data/resources/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv

## Run DESeq and fGSEA analyses:
Rscript r/do/run_deseq.R
```

### Citation
- URL: **pending**
- DOI: **pending**

### Contact
E-mail any questions to [gorelica@mskcc.org](mailto:gorelica@mskcc.org?subject=[GitHub]%20mtDNA-Mutations%20paper).
