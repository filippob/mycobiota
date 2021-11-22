#' ---
#' title: "Data analysis on mycobiota data of USEFUL project"
#' author: Francesco Vitali
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---

#' ```{r global_options, echo = FALSE, include = FALSE}
#' options(width = 9999)
#' knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE,
#'                      cache = FALSE, tidy = FALSE, size = "small",
#'                      fig.height = 10, fig.width = 15,fig.align = "center")
#' ```

#' # INTRODUCTION
#' 
#' In this script we will perform classical alpha and beta diversity analysis on the results of mycobiota sequencing 
#' of Taleggio cheese. We will perform analysis on the R1 and R2 separately (i.e. using picking results from non merged reads)
#' then comparing results (i.e. Mantel test or procsutes rotation). Finally, we will explore sthe effect of some experimental metadata
#' 
#' 

#+ echo=FALSE, message = FALSE

# Load libraries
library(phyloseq)
library(ggplot2)
library(reshape2)
library(ape)
library(data.table)
library(ggsci)
library(ggpubr)
library(picante)
library(microbiome)
library(cowplot)
library(metagenomeSeq)
library(pals)
library(dplyr)
library(tidyr)
library(knitr)
library(pairwiseAdonis)
library(stringr)
library(patchwork)
library(ggrepel)

# Loading data 

# Loading R1 data

metadata_Taleggio <- read.table(file = ".", header = T, sep = "\t", row.names = 1)
otutable_Taleggio <- read.table(file = "./MEATIC_WP3WP4_rawdata/16S/MICCA_16S_WP3WP4/otutable.txt", header = T, sep = "", row.names = 1)
taxonomy_Taleggio <- read.table(file = "./MEATIC_WP3WP4_rawdata/16S/MICCA_16S_WP3WP4/taxa_SILVA_R.csv", sep = "\t", row.names = 1)
colnames(taxonomy_bact_WP3WP4) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
row.names(metadata_bact_WP3WP4) <- gsub("-", ".", row.names(metadata_bact_WP3WP4)) # R converts the "-" in the otu table as "."

# reading fasta file
fastaFile <- readDNAStringSet("./MEATIC_WP3WP4_rawdata/16S/MICCA_16S_WP3WP4/otus.fasta")
length(fastaFile)
summary(width(fastaFile)) ## all reads are 400bp, so the whole md5 hash makes sense
seq_name = names(fastaFile)
sequence = paste(fastaFile)
fasta_bact_WP3WP4 <- data.frame(seq_name, sequence)
# obtaining md5 hash of the sequence
library(digest)
fasta_bact_WP3WP4$md5 <- sapply(fasta_bact_WP3WP4$sequence, digest, algo="md5")
# use the md5 hash of the sequence to rename OTU in the different objects (taxonomy and otutable)
summary(row.names(otutable_bact_WP3WP4) == fasta_bact_WP3WP4$seq_name) # check if order is the same
summary(row.names(taxonomy_bact_WP3WP4) == fasta_bact_WP3WP4$seq_name) # check if order is the same
row.names(otutable_bact_WP3WP4) <- fasta_bact_WP3WP4$md5
row.names(taxonomy_bact_WP3WP4) <- fasta_bact_WP3WP4$md5

summary(unique(fasta_bact_WP3WP4$md5))
summary(fasta_bact_WP3WP4$md5)


