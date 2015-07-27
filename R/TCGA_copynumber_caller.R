library(CGHcall)
library(IRanges)
library(GenomicRanges)
library(survival)
library(biomaRt)
library(Biobase)
#library(QDNAseq)

#Get list of genes of a region
ensembl54=useMart(host='may2009.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL')
ensembl <- useDataset("hsapiens_gene_ensembl", ensembl54)
attributes <- c("ensembl_gene_id", "external_gene_id", "entrezgene", "chromosome_name", "start_position", "end_position", "percentage_gc_content")
#on chromosome
#filters <- c("chromosome_name")
#values <- list(chromosome="10")
#on chromosome part
#filters <- c("chromosome_name","start","end")
#values <- list(chromosome="10",start="112939991",end="135284990")
#on gene name
filters <- c("hgnc_symbol")
values <- list(external_gene_id=genes)

geneList <- getBM(attributes=attributes, filters=filters, values=values, mart=ensembl)
geneListAll <- geneList

#part of genome
#       chr10 135374737
#       chr11 134452384
#geneList[1,] <- c("Distal 10q", "Distal 10q","Distal 10q", "10", "112939991","135284990","1" )
#geneList <- geneList[-c(2:nrow(geneList)),] 
#geneList[2,] <- c("Chromosome 10q", "Chromosome 10q","Chromosome 10q", "10", "42100000","135374737","1" )
#geneList$start_position <- as.numeric(geneList$start_position)
#geneList$end_position <- as.numeric(geneList$end_position)

#chr13
#chr13	114142980
#chr14  106368585
#geneList[1,] <- c("Chromosome 13q", "Chromosome 13q","Chromosome 13q", "13", "18400000","114142980","1" )
#geneList <- geneList[-c(2:nrow(geneList)),]
#geneList$start_position <- as.numeric(geneList$start_position)
#geneList$end_position <- as.numeric(geneList$end_position)

#save(geneListAll, file="TCGAllg.Rdata")

# remove all gene without entrezID
geneList <- geneList[!is.na(geneList[,3]),]

# remove all genes not mapping to chr 1 to 22 or X and Y
geneList <- geneList[which(geneList[,4] %in% c(1:22, "X", "Y")),]
geneList <- geneList[which(geneList[,4] %in% c(values)),]
geneList <- geneList[,-1]
geneList <- geneList[,-2]
geneList <- unique(geneList)
geneList <- geneList[order(geneList[,2], geneList[,3], geneList[,4]),]

input_genes <- geneList$external_gene_id

print_segments = "FALSE"

if(length(which(input_genes%in%genes != TRUE))>0){
  print (paste("gene", input_genes[which(input_genes%in%genes != TRUE)], "not found"))
  input_genes <- input_genes[-which(input_genes%in%genes != TRUE)]
} else {print("all genes found, extracting deletion status")}

### Specific for Copynumber

####LOAD FILES
#setwd("/Users/martijncordes/Dropbox/Werk/TCGA/TCGA_COAD_GSNP6noCNV-2015-02-24/")
#setwd("/home/cordes/Dropbox/Werk/TCGA/TCGA_LUAD_GSNP6noCNV-2015-02-24/")
#setwd("/home/cordes/Dropbox/Werk/TCGA/TCGA_BRCA_GSNP6noCNV-2015-02-24/")
#setwd("/home/cordes/Dropbox/Werk/TCGA/TCGA_LUSC_GSNP6noCNV-2015-02-24/")
#setwd("/home/cordes/Dropbox/Werk/TCGA/TCGA_LGG_GSNP6noCNV-2015-02-24/")
#setwd("/home/cordes/Dropbox/Werk/TCGA/TCGA_HNSC_GSNP6noCNV-2015-02-24/")
#setwd("/home/cordes/Dropbox/Werk/TCGA/TCGA_STAD_GSNP6noCNV-2015-02-24//")
setwd("/Users/martijncordes/Dropbox/Werk/TCGA/TCGA_HNSC_GSNP6noCNV-2015-02-24/")




#setwd("/home/cordes/Dropbox/Werk/TCGA/old_data/TCGA_LGG_GSNP6noCNV-2014-08-22/")

#Tumorscape
#setwd("/home/cordes/Dropbox/Werk/TCGA/Tumorscape_Beatriz/")


source("../TCGA_functions.R")
source("/home/cordes/Dropbox/Werk/TCGA/TCGA_functions.R")

###Clinical data
copynumber_clinical_data <- read.delim("clinical_data", header=T)
colnames(copynumber_clinical_data) <- gsub("X_","_",colnames(copynumber_clinical_data), fixed=TRUE)

#Remove MSI
#no_msi <- c(which(copynumber_clinical_data$MSI_updated_Oct62011 == "MSS"), which(copynumber_clinical_data$microsatellite_instability == "NO"))
#no_msi <- no_msi[-which(duplicated(no_msi))]

#copynumber_clinical_data <- copynumber_clinical_data[no_msi,]

#remove patients with no OS
if (length(which(is.na(copynumber_clinical_data[,which(colnames(copynumber_clinical_data)=="_OS")]))) > 0){
  copynumber_clinical_data <- copynumber_clinical_data[-which(is.na(copynumber_clinical_data[,which(colnames(copynumber_clinical_data)=="_OS")])),]
}

#For LGG set
#Filter called data for 10q region
#hg19: chr10:112950001-135435000/hg18:chr10:112939991-135284990)

#filter for clinical data with survival values
#grade2 <- as.vector(copynumber_clinical_data$sampleID[which(copynumber_clinical_data$neoplasm_histologic_grade == "G2")])
#grade3 <- as.vector(copynumber_clinical_data$sampleID[which(copynumber_clinical_data$neoplasm_histologic_grade == "G3")])

#os_grade2 <- copynumber_clinical_data[which(copynumber_clinical_data$neoplasm_histologic_grade == "G2"),]
#os_grade3 <- copynumber_clinical_data[which(copynumber_clinical_data$neoplasm_histologic_grade == "G3"),]

#copynumber_clinical_data <- os_grade3

### Segment Data
segment_data <- read.delim("genomicSegment", header=F)
copynumber_samples <- as.vector(unique(segment_data[,1]))

#remove CN copynumber_samples without OS status
copynumber_samples <- copynumber_samples[which(copynumber_samples %in% copynumber_clinical_data$sampleID)]

#remove copynumber_samples from clinical data with no CN status
copynumber_clinical_data <- copynumber_clinical_data[which(copynumber_clinical_data$sampleID %in% copynumber_samples),]

#only StageIV
#copynumber_samples <- as.vector(copynumber_clinical_data$sampleID[c(which(copynumber_clinical_data$pathologic_stage == "Stage IV"), which(copynumber_clinical_data$pathologic_stage == "Stage IVA"))])
#copynumber_clinical_data <- copynumber_clinical_data[c(which(copynumber_clinical_data$pathologic_stage == "Stage IV"), which(copynumber_clinical_data$pathologic_stage == "Stage IVA")),]


#Create tumor dataset from all copynumber_samples
tumor_segments <- data.frame(Name=segment_data[,1], Chromosome=segment_data[,2], Start=segment_data[,3], End=segment_data[,4], log2=segment_data[,6])

#correct for log2ratio
#tumor_segments_correct <- tumor_segments_all[which(tumor_segments_all$log2 < 1),]
#tumor_segments_correct <- tumor_segments_correct[which(tumor_segments_all$log2 > -1),]
#tumor_segments <- tumor_segments_correct
#tumor_segments <- tumor_segments_all


#Tumorscape methods
#segment_data <- read.delim("tumorscape_100217.seg", header=T)
#tumor_segments <- data.frame(Name=segment_data[,1], Chromosome=paste("chr",segment_data[,2],sep=""), Start=segment_data[,3], End=segment_data[,4], log2=segment_data[,6])

#copynumber_samples <- read.delim("sample_info_current.txt")
#copynumber_samples <- as.vector(copynumber_samples$Array[which(copynumber_samples$Type == "Colorectal")])


#Divide chromosomes in bins
#10kb  = 10000
#100kb = 100000
#1MB   = 1000000
#10MB  = 10000000    
bin_size <- 10000
build = "HG18"
bins <- makeBins(bin_size,build)

#tumor_percentage_range <- seq(10,30,10)
tumor_percentage_range <- 30

for (percentage in 1:length(tumor_percentage_range)) {
  if (file.exists("copynumbers_called.Rdata") == TRUE) {
    load("copynumbers_called.Rdata")
  }
  
  else {
    tumor_percentage <- tumor_percentage_range[percentage]
    print ( paste("Tumor percentage ",tumor_percentage,"pct",sep =""))
    
    #Call copynumbers
    copynumbers_called <- SegmentCalling(copynumber_samples, tumor_segments, tumor_percentage)
    
    if (bin_size <= 100000) {
      #Calculate overlap of gene regions with bins (only for bin sizes smaller then 100kb)
      for (percent in seq(90,90,10)) {
        copynumber_clinical_data <- GeneCalling(input_genes, copynumbers_called, percent)
        km_info <- KM_per_gene(input_genes,copynumber_clinical_data)
    
      }
    }
  
  }
}

#plots
#make cghRegions object
regions_object <- make_cghRegions(copynumbers_called)

#Make frequency Plot
print("Making frequency plot")
png(paste("plots/frequencyplot_cutoff",tumor_percentage[percentage], "pct.png",sep=""), width=297, height=210, units='mm', res=150)
frequencyPlot(regions_object, build=build, main=paste("Frequency plot of TCGA stomach adenocarcinoma samples (n=", length(copynumber_samples), ")", sep=""))
dev.off()

#Make Correlation plot
cgh.reg <- regions_object

y <- regions(cgh.reg)
y[y == 2] <- 1
y[y == -2] <- -1
m <- cor(t(y), method='pearson',use='complete.obs')

png("plots/correlation_plot.png", width=297, height=210, units='mm', res=150)
par(mar=c(5,4,4,4)+.1)
cor.mat.by.region(cgh.reg, m, build, main=paste('Correlation matrix for Tumorscape Colorectal set (n=', ncol(cgh.reg), ')', sep=''))
dev.off()

cor.mat <- cor(t(y),method='spearman',use='complete.obs')
pdf("plots/spearman_correlation.pdf",width=8,height=7.5)
par(oma=c(5,0,0,0))
heatmap(cor.mat,symm=T,col = maPalette(high="blue",low="white"))
dev.off()

#SV samples
SVsamples <- copynumber_samples[1:length(copynumber_samples)]
interesting_samples <- c()
sample_list <- c()
for (i in 1:length(SVsamples)){
  tumor_dataset <- tumor_segments_all[which(tumor_segments_all$Name == SVsamples[i]),]
  segments <- which(tumor_dataset$log2 > 1)
  if (length(segments) > 5) {
    chroms <- table(tumor_dataset$Chromosome[segments])
      for(chr in 1:length(chroms)) {
        if (chroms[chr] > 10 && !(i %in% sample_list) ){
           interesting_samples <- c(interesting_samples,SVsamples[i] )
           sample_list <- c(sample_list, i)
        }
      } 
   }
}

for (i in 1:length(interesting_samples)) {
  tumor_dataset <- tumor_segments_all[which(tumor_segments_all$Name == interesting_samples[i]),]
  saveRDS(tumor_dataset, "segments.rds")
  
  source("create aberrations for Martijn.R")
  
  colnames(arrayData) <- interesting_samples[i]
  png(paste("",interesting_samples[i],"_2.png",sep=""), width=297, height=210, units='mm', res=150)
  plot(arrayData, dotres=1, ylim=c(-5,5), main="TEST", build=genome_build)
  dev.off()

}

#Segment plot
#Print Segments
if (print_segments == "TRUE") {
  plot_samples <- sample_list
  PrintSegments(plot_samples, "HG18")    
}



#TEST STUFF
#Survival for SV samples 
name <- "unusual amount of segments (>10) per chromosome"
copynumber_clinical_data[,name] <- "no_SV"
for (j in 1:length(interesting_samples)) {
  copynumber_clinical_data[,name][which(copynumber_clinical_data$sampleID == interesting_samples[j])] <- paste("a_SV_deletion", sep="")  
}

input_genes <- "unusual amount of segments (>10) per chromosome"
gene <- 1
km_info <- KM_per_gene(input_genes,copynumber_clinical_data)


