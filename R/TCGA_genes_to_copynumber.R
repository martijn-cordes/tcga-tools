library(Biobase)
library(QDNAseq)
library(IRanges)
library(GenomicRanges)
source("~/Dropbox/Werk/TCGA/TCGA_functions.R")

setwd("~/Dropbox/Werk/TCGA/TCGA_CESC_GSNP6noCNV-2015-02-24/")

###Clinical data
copynumber_clinical_data <- read.delim("clinical_data", header=T)
colnames(copynumber_clinical_data) <- gsub("X_","_",colnames(copynumber_clinical_data), fixed=TRUE)

### Segment Data
segment_data <- read.delim("genomicSegment", header=F)
copynumber_samples <- as.vector(unique(segment_data[,1]))

#Create tumor dataset from all copynumber_samples
tumor_segments <- data.frame(Name=segment_data[,1], Chromosome=segment_data[,2], Start=segment_data[,3], End=segment_data[,4], log2=segment_data[,6])

#Make bins
bin_size <- 15000
build = "HG19"
bins <- makeBins(bin_size,build)

#Copynumber frame
copynumbers_frame <- SegmentFrame(copynumber_samples, tumor_segments)

#Expression data
setwd("/Users/martijncordes/Dropbox/Werk/TCGA/TCGA_CESC_exp_HiSeqV2-2015-02-24/")

### Expression Data
expression_data <- read.delim("genomicMatrix", header=T)
colnames(expression_data) <- gsub(".", "-", colnames(expression_data), fixed=TRUE)
genes <- as.vector((expression_data[,1]))
samples <- colnames(expression_data)

### Probe Data
probe_data <- read.delim("probe", stringsAsFactors = F)
gene_list <- probe_data[which(probe_data$gene %in% genes),]
gene_list <- gene_list[,1:6]
input_genes <- as.vector(unique(gene_list$gene))
not_in_set <- which(!genes %in% input_genes)

duplicates <- gene_list$gene[duplicated(gene_list$gene)]
gene_list[gene_list$gene %in% duplicates,]


#gene_copynumber
gene_copynumber <- GeneToCopynumber(gene_list, copynumbers_frame, copynumber_samples)

#Mutation data
setwd("/Users/martijncordes/Dropbox/Werk/TCGA/TCGA_CESC_mutation_bcgsc_gene-2015-02-24/")

### Expression Data
mutation_data <- read.delim("genomicMatrix", header=T)
colnames(mutation_data) <- gsub(".", "-", colnames(mutation_data), fixed=TRUE)
mut_genes <- as.vector((mutation_data[,1]))
mut_samples <- colnames(mutation_data)


### Probe Data
probe_data <- read.delim("probe", stringsAsFactors = F)
gene_list <- probe_data[which(probe_data$gene %in% mut_genes),]
gene_list <- gene_list[,1:6]
input_genes <- as.vector(unique(gene_list$gene))
not_in_set <- which(!mut_genes %in% input_genes)

duplicates <- gene_list$gene[duplicated(gene_list$gene)]
gene_list[gene_list$gene %in% duplicates,]


#gene_copynumber
mut_gene_copynumber <- GeneToCopynumber(gene_list, copynumbers_frame, copynumber_samples)


#Make one genefile for all genes used for expression and mutation
gene_copynumber[which(!gene_copynumber$gene %in% mut_gene_copynumber$gene),]

all_gene_copynumber <- rbind(mut_gene_copynumber, gene_copynumber[which(!gene_copynumber$gene %in% mut_gene_copynumber$gene),] )

saveRDS(copynumber_clinical_data, file="TCGA_CESC_clinical_data.rds")
saveRDS(expression_data, file="TCGA_CESC_expression_data.rds")
saveRDS(mutation_data, file="TCGA_CESC_mutation_data.rds")
saveRDS(all_gene_copynumber, file="TCGA_CESC_copynumber_per_gene.rds")

