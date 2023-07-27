library(data.table)
library(devtools)
document()

# To Do -------------------------------------------------------------------

# make tests
# split by CHROM for "in" functions for speed and memory use
# more features for read.VCF that
  #parses by input type
  #option for long format (for multiple samples)
  #deepvariant

#wrapper function for Strelka2 function
  #directory argument
  #directory structure: dir/NORMAL/TUMOR


source("R/read.bedGraph.R")
source("R/features_in_features.R")
source("R/bedGraph_total.R")


input<-read.bedGraph("~/Dropbox/Research/rice mutation paper/data/GSM4668650_Col0_rep1_H3K4me1_Input_unique_reads.bedGraph.gz")
genes_input<-features_in_features(features = genes, features2 = input, mode = "sumxlength", value="DEPTH")
colnames(genes_input)[2]<-"input"

chip<-read.bedGraph("~/Dropbox/Research/rice mutation paper/data/GSM4668649_Col0_rep1_H3K4me1_ChIP_unique_reads.bedGraph.gz")
genes_chip<-features_in_features(features = genes, features2 = chip, mode = "sumxlength", value="DEPTH")
colnames(genes_chip)[2]<-"chip"
merge<-merge(genes_chip, genes_input)

chip_total<-bedGraph_total(chip)
input_total<-bedGraph_total(input)

merge$enrich<-log2((1+merge$chip)/chip_total) - log2((1+merge$input)/input_total)
return(merge)

source("R/features_chip_enrich.R")

chipfile = "~/Dropbox/Research/rice mutation paper/data/GSM4668649_Col0_rep1_H3K4me1_ChIP_unique_reads.bedGraph.gz"
inputfile = "~/Dropbox/Research/rice mutation paper/data/GSM4668650_Col0_rep1_H3K4me1_Input_unique_reads.bedGraph.gz"


genes_chip<-features_chip_enrich(genes, chipfile = chipfile, inputfile = inputfile)

source("R/feature_windows.R")
source("R/features_in_features.R")
source("R/plot_feature_windows.R")
genes$DIRECTION<-genes$direction
gene_windows<-feature_windows(features = genes, breaks = 10, directed = T, ID="gene")
gene_windows_chip<-features_chip_enrich(gene_windows, chipfile = chipfile, inputfile = inputfile)
gene_windows$enrich<-gene_windows_chip$enrich

source("R/strelka2_filter.R")
source("R/read.VCF.R")

strelka2_filter<-strelka2_filter(directory = "~/Desktop/strelka_out/", output = "~/Desktop/")

gene_windows_mutations<-sites_in_features(features = gene_windows, sites = strelka2_filter, mode = "counts")
gene_windows$mutations<-gene_windows_mutations$counts

source("R/sites_in_features.R")

plot_feature_windows(gene_windows,variable = "enrich",mode="mean")
plot_feature_windows(gene_windows,variable = "mutations",mode="percent")


genes<-fread("~/Dropbox/Research/rice mutation paper/data/A_thal_genes_PDS5_enrich.csv")
#genes<-gff[type=="gene"]
genes$CHROM<-as.character(genes$chr)
genes$ID<-genes$gene
genes$START<-genes$start
genes$STOP<-genes$stop

source("R/read.GFF.R")
GFF<-read.GFF("~/Dropbox/Research/rice mutation paper/data/all.gff3")

mutations<-fread("~/Dropbox/Research/rice mutation paper/data/MSH6_KO_mutations.csv")


source("R/sites_in_features.R")
sites_in_features(features=genes,
                  sites=mutations,
                  mode="counts")
mutations$depth_pct
depth<-sites_in_features(features=genes,
                  sites=mutations,
                  mode="mean",
                  value="depth_pct")



features$CHROM<-as.character(features$CHROM)
sites$CHROM<-as.character(sites$CHROM)
