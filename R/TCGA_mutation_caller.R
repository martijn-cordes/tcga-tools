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
attributes <- c("ensembl_gene_id", "external_gene_id", "entrezgene", "chromosome_name", "start_position", "end_position", "percentage_gc_content")
#on chromosome
#filters <- c("chromosome_name")
#values <- list(chromosome="10")
#on chromosome part
#filters <- c("chromosome_name","start","end")
#values <- list(chromosome="10",start="112950001",end="135435000")
#on gene name
filters <- c("hgnc_symbol")
values <- list(external_gene_id=c("APC","KRAS","BRAF","NRAS","PIK3CA","TP53","SMAD4","FBXW7"))
geneList <- getBM(attributes=attributes, filters=filters, values=values, mart=ensembl)
geneListAll <- geneList
#save(geneListAll, file="TCGAllg.Rdata")

# remove all gene without entrezID
geneList <- geneList[!is.na(geneList[,3]),]

# remove all genes not mapping to chr 1 to 22 or X and Y
geneList <- geneList[which(geneList[,4] %in% c(1:22, "X", "Y")),]
#geneList <- geneList[which(geneList[,4] %in% c(values)),]
geneList <- geneList[,-1]
geneList <- geneList[,-2]
geneList <- unique(geneList)
geneList <- geneList[order(geneList[,2], geneList[,3], geneList[,4]),]

input_genes <- geneList$external_gene_id



### Specific for Mutation
 
#set source for files
setwd("/home/cordes/Dropbox/Werk/TCGA/old_data/TCGA_COAD_mutation_bcm_gene-2014-08-27/")
#setwd("/home/cordes/Dropbox/Werk/TCGA/TCGA_STAD_mutation-2014-05-02/")
#setwd("/home/cordes/Dropbox/Werk/TCGA/TCGA_LGG_mutation_curated_broad_gene-2014-08-27/")
#setwd("/Users/martijncordes/Dropbox/Werk/TCGA/TCGA_LGG_mutation_broad_gene-2014-08-27/")
#setwd("/home/cordes/Dropbox/Werk/TCGA/TCGA_LGG_mutation_broad_gene-2014-08-27/")

####LOAD FILES

### Clinical data
mut_clinical_data <- read.delim("clinical_data", header=T)
colnames(mut_clinical_data) <- gsub("X_","_",colnames(mut_clinical_data), fixed=TRUE)
#mut_clinical_data$sampleID <- substr(mut_clinical_data$sampleID, 1, 12)

#Remove MSI
no_msi <- c(which(mut_clinical_data$MSI_updated_Oct62011 == "MSS"), which(mut_clinical_data$microsatellite_instability == "NO"),which(mut_clinical_data$MSI_updated_Oct62011 == ""))
no_msi <- no_msi[-which(duplicated(no_msi))]

mut_clinical_data <- mut_clinical_data[no_msi,]


grade2 <- as.vector(mut_clinical_data$sampleID[which(mut_clinical_data$tumor_stage %in% c("Stage II","Stage IIA", "Stage IIB", "Stage IIC"))])
grade3 <- as.vector(mut_clinical_data$sampleID[which(mut_clinical_data$tumor_stage %in% c("Stage III","Stage IIIA", "Stage IIIB", "Stage IIIC"))])

os_grade2 <- mut_clinical_data[which(mut_clinical_data$tumor_stage %in% c("Stage II","Stage IIA", "Stage IIB", "Stage IIC")),]
os_grade3 <- mut_clinical_data[which(mut_clinical_data$tumor_stage %in% c("Stage III","Stage IIIA", "Stage IIIB", "Stage IIIC")),]

#mut_clinical_data <- os_grade3


### mutation Data
mutation_data <- read.delim("genomicMatrix", header=T)
genes <- as.vector(unique(mutation_data[,1]))
samples <- colnames(mutation_data)
samples <- gsub(".", "-", samples, fixed=TRUE)
#samples <- substr(samples, 1, 12)

input_genes <- geneList$external_gene_id
#input_genes <- c("KRAS","ARID1A","PTEN")

if(length(which(input_genes%in%genes != TRUE))>0){
  print (paste("gene", input_genes[which(input_genes%in%genes != TRUE)], "not found"))
  input_genes <- input_genes[-which(input_genes%in%genes != TRUE)]
} else {print("all genes found, extracting mutation status")}
  
#remove samples from clinical data with no mutation status
mut_clinical_data <- mut_clinical_data[which(mut_clinical_data$sampleID %in% samples),]

#remove patients with no OS
if (length(which(is.na(mut_clinical_data[,which(colnames(mut_clinical_data)=="_OS")]))) > 0){
  mut_clinical_data <- mut_clinical_data[-which(is.na(mut_clinical_data[,which(colnames(mut_clinical_data)=="_OS")])),]
}

#only StageIV
#samples <- as.vector(mut_clinical_data$sampleID[c(which(mut_clinical_data$pathologic_stage == "Stage IV"), which(mut_clinical_data$pathologic_stage == "Stage IVA"))])
#mut_clinical_data <- mut_clinical_data[c(which(mut_clinical_data$pathologic_stage == "Stage IV"), which(mut_clinical_data$pathologic_stage == "Stage IVA")),]

  
mutation_info <- data.frame(Gene=0,Pvalue=0, Mutation=0, No_Mutation=0)
for (i in 1:length(input_genes)) {
  #Mutation vs no mutation
  no_mutation <- colnames(mutation_data[which(mutation_data[which(genes==input_genes[i]),]==0)])
  no_mutation <- gsub(".", "-", no_mutation, fixed=TRUE)
  
  mutation <- colnames(mutation_data[which(mutation_data[which(genes==input_genes[i]),]==1)])
  mutation <- gsub(".", "-", mutation, fixed=TRUE)
  
  
  #add clinical parameter (mutation or no mutation) to clinical data
  
  mut_clinical_data[,input_genes[i]] <- 0
  for (j in 1:length(mutation)) {
    mut_clinical_data[,input_genes[i]][which(mut_clinical_data$sampleID == mutation[j])] <- paste("a_",input_genes[i], "_mutation", sep="")
    
  }
  
  for (j in 1:length(no_mutation)) {
    mut_clinical_data[,input_genes[i]][which(mut_clinical_data$sampleID == no_mutation[j])] <- paste("no_",input_genes[i], "_mutation", sep="")
    
  }
  
  #remove zeros
  if (length(which(mut_clinical_data[,input_genes[i]] == 0)) > 0){
    mut_clinical_data <- mut_clinical_data[-which(mut_clinical_data[,input_genes[i]] == 0),]
  }
  
  #KM
  column_number <- which(colnames(mut_clinical_data)==input_genes[i])
  mut_clinical_data$SurvObj <- with(mut_clinical_data, Surv(as.numeric(as.vector(mut_clinical_data[,which(colnames(mut_clinical_data)=="_OS")])) ,as.numeric(as.vector(mut_clinical_data[,which(colnames(mut_clinical_data)=="_EVENT")]) == 1)))
  km.by.gene <- survfit(SurvObj ~ mut_clinical_data[,column_number], data = mut_clinical_data, conf.type = "log-log")
  rate <- (km.by.gene$n[1]/sum(km.by.gene$n))*100
  
  #more than one groupe to test
  if (length(km.by.gene$n) > 1) {
    pval <- survdiff(SurvObj ~ mut_clinical_data[,column_number], data = mut_clinical_data)
    pval.by.gene <- 1-pchisq(pval$chisq,length(pval$n)-1)
    
    #pvalue below 0.05 and mutation rate above 5%
    if (pval.by.gene <= 0.9 || rate > 0) {
      print(paste("Creating KM for",input_genes[i]))
      mutation_info <- rbind(mutation_info, c(input_genes[i],pval.by.gene, km.by.gene$n[1],km.by.gene$n[2]))
        
      #Make KM plot 
      png(paste("plots/",input_genes[i],".png",sep=""), width=297, height=210, units='mm', res=150)
      plot(km.by.gene, mark.time=T, xlab='Days', ylab='Survival', col=c("dodgerblue4","gold2"),lwd=2, lty=1:2, yaxt="n", main=paste(input_genes[i], " in TCGA Colon adenocarcinoma's  (n = ", nrow(mut_clinical_data), ")",sep=""))
      #plot(km.by.gene, mark.time=T, xlab='Days', ylab='Survival', col=c("dodgerblue4","gold2"),lwd=2, lty=1:2, yaxt="n", main=paste(input_genes[i], " in TCGA Stage IV Colon adenocarcinoma's  (n = ", nrow(mut_clinical_data), ")",sep=""))
      axis(1, at = c(500,1500,2500,3500,4500,5500,6500))
      axis(2, at = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), las=1)
      legend("topright", c(paste(input_genes[i],"mutation","(n =",length(which(mut_clinical_data[,input_genes[i]]==paste("a_",input_genes[i], "_mutation", sep=""))),")"),paste("no", input_genes[i], "mutation","(n =",length(which(mut_clinical_data[,input_genes[i]]==paste("no_",input_genes[i],"_mutation",sep=""))),")")), col=c("dodgerblue4","gold2"), lwd=2)
      legend("bottomright",inset = 0.05, paste("pvalue =",pval.by.gene),box.col=NA )
      dev.off()
    
    }
  }
}

mutation_info <- mutation_info[-1,]

















