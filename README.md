# Respiratory complex and tissue lineage drive recurrent mutations in tumour mtDNA

![alt text](https://github.com/reznik-lab/mtdna-mutations/blob/master/data/wiki.png "Title")

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

The generated file `html/figures.html` can be opened in a web browser to view the panels from the main text and extended data figures.


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
## - data/processed_data/samples_callable_per_position.txt
## - data/processed_data/sample_gene_lengths.txt.gz
## Additionally generate similar files using a cutoff of 20 alt-reads rather than the default 5:
## - data/processed_data/samples_callable_per_position20.txt
## - data/processed_data/sample_gene_lengths20.txt.gz
Rscript r/do/process_sample_coverage.R
```

2. Generate table of fraction of mtDNA with callable mutations in tumor, normal, and combined tumor+normal samples, at a threshold of 5+ alt reads, and 20+ alt reads: 
```shell
## Generates: data/processed_data/frac_callable.txt
Rscript r/do/get_frac_callable.R
```

3. Determine recurrent positions for mtDNA SNV, truncating indels, and tRNA mutations: 
```shell
## Calculate SNP hotspots. Generates: 
## - data/processed_data/hotspots_snv.txt
Rscript r/do/hotspot_snv.R

## Calculate frame-shift indel hotspots. Generates:
## - data/processed_data/hotspots_indel.txt
Rscript r/do/hotspot_indel.R

## Calculate tRNA alignment hotspots. Generates:
## - data/processed_data/hotspots_trna.txt
## - data/processed_data/hotspots_trna_pcawg.txt ## for figure EDF7b
Rscript r/do/hotspot_trna.R
```

4. Annotate samples with mtDNA variants, % callable, and overall mtDNA-status:
```shell
## Generates: 
## - data/processed_data/sample_hotspot_status_tcga.txt
## - data/processed_data/sample_hotspot_status_impact.txt
## - data/processed_data/hotspot_region_map.txt
Rscript r/do/annotate_sample_mutations.R
```

5. Calculate differentially-expressed genes and run GSEA:
```shell
## Generates 
## - data/processed_data/rnaseq_truncating_gsea_hallmark_precalculated.txt
## - data/processed_data/rnaseq_truncating_histogram_precalculated.txt
## - data/processed_data/rnaseq_vus_gsea_hallmark_precalculated.txt
## - data/processed_data/rnaseq_vus_histogram_precalculated.txt
## - data/processed_data/rnaseq_r25q_gsea_hallmark_precalculated.txt

## First, download TCGA RNA-Seq RSEM data (1.8GB file):
wget http://api.gdc.cancer.gov/data/3586c0da-64d0-4b74-a449-5ff4d9136611 -O data/original_data/resources/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv

## Now run DESeq and fGSEA analyses:
Rscript r/do/run_deseq.R
```

### Citation
- URL: **pending**
- DOI: https://doi.org/10.1038/s42255-021-00378-8

### Contact
E-mail any questions to [gorelica@stanford.edu](mailto:gorelica@stanford.edu?subject=[GitHub]%20mtDNA-Mutations%20paper).
