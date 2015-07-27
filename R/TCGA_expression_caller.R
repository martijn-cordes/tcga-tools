#library(QDNAseq)
#library(Biobase)
library(CGHcall)
library(IRanges)
library(GenomicRanges)
library(survival)
#source("frequencyPlot.R")
library(biomaRt)

#Get list of genes of a region
ensembl54=useMart(host='may2009.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL')
ensembl <- useDataset("hsapiens_gene_ensembl", ensembl54)
attributes <- c("ensembl_gene_id", "external_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "percentage_gc_content")

#ensembl <- useMart(biomart="ensembl")
#ensembl <- useDataset("hsapiens_gene_ensembl", ensembl)
#attributes <- c("ensembl_gene_id", "hgnc_symbol", "entrezgene", "chromosome_name", "start_position", "end_position", "percentage_gc_content")

#on chromosome
#filters <- c("chromosome_name")
#values <- list(chromosome="10")
#on chromosome part
#filters <- c("chromosome_name","start","end")
#values <- list(chromosome="10",start="112950001",end="135435000")

#on gene name
filters <- c("hgnc_symbol")
values <- list(external_gene_id=duplicates)
geneList <- getBM(attributes=attributes, filters=filters, values=values, mart=ensembl)
geneList
geneListAll <- geneList
#save(geneListAll, file="TCGAllg.Rdata")

# remove all gene without entrezID
geneList <- geneList[!is.na(geneList[,3]),]

# remove all genes not mapping to chr 1 to 22 or X and Y
geneList <- geneList[which(geneList[,4] %in% c(1:22, "X", "Y")),]
#geneList <- geneList[which(geneList[,4] %in% c(values)),]
geneList <- geneList[,-1]
#geneList <- geneList[,-2]
geneList <- unique(geneList)
geneList <- geneList[order(geneList[,3], geneList[,4], geneList[,5]),]

input_genes <- geneList$hgnc_symbol

not_found <- genes[which(genes%in%input_genes != TRUE)]

if(length(which(input_genes%in%genes != TRUE))>0){
  print (paste("gene", input_genes[which(input_genes%in%genes != TRUE)], "not found"))
  input_genes <- input_genes[-which(input_genes%in%genes != TRUE)]
} else {print("all genes found, extracting deletion status")}

hgnc <- read.delim(url("http://www.genenames.org/cgi-bin/download?col=gd_app_sym&col=gd_prev_sym&col=gd_aliases&status=Approved&status=Entry+Withdrawn&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&hgnc_dbtag=on&submit=submit"))

for (i in 1:length(not_found)) {
  print(as.vector(hgnc$Previous.Symbols[grep(not_found[i],hgnc$Previous.Symbols)]))
}

as.vector(hgnc[grep(not_found[1],hgnc$Approved.Symbol),])
as.vector(hgnc$Previous.Symbols[grep(not_found[1],hgnc$Previous.Symbols)])

#set source for files
#setwd("/home/cordes/Dropbox/Werk/TCGA/TCGA_LGG_exp_HiSeqV2-2014-08-28/")
#setwd("/home/cordes/Dropbox/Werk/TCGA/old_data/TCGA_LGG_exp_HiSeqV2-2014-08-28/")
#setwd("/home/cordes/Dropbox/Werk/TCGA/TCGA_COAD_exp_HiSeqV2-2014-08-28/")
#setwd("/home/cordes/Dropbox/Werk/TCGA/TCGA_STAD_exp_HiSeq-2014-08-28/")
#setwd("/home/cordes/Desktop/TCGA_LGG_exp_HiSeqV2_exon-2014-08-28/")
setwd("/home/cordes/Dropbox/Werk/TCGA/TCGA_HNSC_exp_HiSeqV2-2015-02-24/")
####LOAD FILES

### Clinical data
expression_clinical_data <- read.delim("clinical_data", header=T)
colnames(expression_clinical_data) <- gsub("X_","_",colnames(expression_clinical_data), fixed=TRUE)

### Expression Data
expression_data <- read.delim("genomicMatrix", header=T)
genes <- as.vector((expression_data[,1])
samples <- colnames(expression_data)-
samples <- gsub(".", "-", samples, fixed=TRUE)

### Probe Data
probe_data <- read.delim("probe", stringsAsFactors = F)
gene_list <- probe_data[which(probe_data$gene %in% genes),]
gene_list <- probe_data[,1:6]
input_genes <- as.vector(unique(gene_list$gene))
not_in_set <- which(!genes %in% input_genes)

duplicates <- gene_list$gene[duplicated(gene_list$gene)]
gene_list[gene_list$gene %in% duplicates,]

if(length(which(input_genes%in%genes != TRUE))>0){
  print (paste("gene", input_genes[which(input_genes%in%genes != TRUE)], "not found"))
  not_found <- c()
  not_found <- c(not_found, input_genes[which(input_genes%in%genes != TRUE)])
  input_genes <- input_genes[-which(input_genes%in%genes != TRUE)]
} else {print("all genes found, extracting expression status")}

#remove samples from clinical data with no expression status
expression_clinical_data <- expression_clinical_data[which(expression_clinical_data$sampleID %in% samples),]

#remove patients with no OS
if (length(which(is.na(expression_clinical_data[,which(colnames(expression_clinical_data)=="_OS")]))) > 0){
  expression_clinical_data <- expression_clinical_data[-which(is.na(expression_clinical_data[,which(colnames(expression_clinical_data)=="_OS")])),]
}

#normalize expression data
expression_matrix <- data.matrix(expression_data[,2:ncol(expression_data)])
normalized_expression <- expression_matrix - rowMeans(expression_matrix)
normalized_expression <- as.data.frame(normalized_expression)

expression_info <- data.frame(Gene=0,Pvalue=0, Expression=0, No_expression=0)
for (i in 1:length(input_genes)) {
  #expression vs no expression
  no_expression <- colnames(normalized_expression[which(normalized_expression[which(genes==input_genes[i]),]<=0)])
  no_expression <- gsub(".", "-", no_expression, fixed=TRUE)
  
  expression <- colnames(normalized_expression[which(normalized_expression[which(genes==input_genes[i]),]>=0)])
  expression <- gsub(".", "-", expression, fixed=TRUE)
  
  
  #add clinical parameter (expression or no expression) to clinical data
  
  expression_clinical_data[,input_genes[i]] <- 0
  for (j in 1:length(expression)) {
    expression_clinical_data[,input_genes[i]][which(expression_clinical_data$sampleID == expression[j])] <- paste("a_",input_genes[i], "_expression", sep="")
    
  }
  
  for (j in 1:length(no_expression)) {
    expression_clinical_data[,input_genes[i]][which(expression_clinical_data$sampleID == no_expression[j])] <- paste("no_",input_genes[i], "_expression", sep="")
    
  }
  
  #remove zeros
  if (length(which(expression_clinical_data[,input_genes[i]] == 0)) > 0){
    expression_clinical_data <- expression_clinical_data[-which(expression_clinical_data[,input_genes[i]] == 0),]
  }
  
  #KM
  column_number <- which(colnames(expression_clinical_data)==input_genes[i])
  expression_clinical_data$SurvObj <- with(expression_clinical_data, Surv(as.numeric(as.vector(expression_clinical_data[,which(colnames(expression_clinical_data)=="_OS")])) ,as.numeric(as.vector(expression_clinical_data[,which(colnames(expression_clinical_data)=="_EVENT")]) == 1)))
  km.by.gene <- survfit(SurvObj ~ expression_clinical_data[,column_number], data = expression_clinical_data, conf.type = "log-log")
  rate <- (km.by.gene$n[1]/sum(km.by.gene$n))*100
  
  #more than one groupe to test
  if (length(km.by.gene$n) > 1) {
    pval <- survdiff(SurvObj ~ expression_clinical_data[,column_number], data = expression_clinical_data)
    pval.by.gene <- 1-pchisq(pval$chisq,length(pval$n)-1)
    
    #pvalue below 0.05 and expression rate above 5%
    #if (pval.by.gene <= 0.05 || rate > 5) {
    if (pval.by.gene <= 0.05) {      
      print(paste("Creating KM for",input_genes[i]))
      expression_info <- rbind(expression_info, c(input_genes[i],pval.by.gene, km.by.gene$n[1],km.by.gene$n[2]))
      
      #Make KM plot 
      png(paste("plots/",input_genes[i],".png",sep=""), width=297, height=210, units='mm', res=150)
      plot(km.by.gene, mark.time=T, xlab='Days', ylab='Survival', col=c("dodgerblue4","gold2"),lwd=2, lty=1:2, yaxt="n", main=input_genes[i])
      axis(1, at = c(500,1500,2500,3500,4500,5500,6500))
      axis(2, at = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), las=1)
      legend("topright", c(paste(input_genes[i],"expression","(n =",length(which(expression_clinical_data[,input_genes[i]]==paste("a_",input_genes[i], "_expression", sep=""))),")"),paste("no", input_genes[i], "expression","(n =",length(which(expression_clinical_data[,input_genes[i]]==paste("no_",input_genes[i],"_expression",sep=""))),")")), col=c("dodgerblue4","gold2"), lwd=2)
      legend("bottomright",inset = 0.05, paste("pvalue =",pval.by.gene),box.col=NA )
      dev.off()
      
    }
  }
}

expression_info <- expression_info[-1,]


