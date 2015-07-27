

##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################



rm(list=ls())
gc()

library(cgdsr)
library(biomaRt)
library(MASS)

##########################################################################################################################################

# get list of human genes
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

# ensembl=useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
#Filter called data for 10q region
#hg19: chr10:112950001-135435000/hg18:chr10:112939991-135284990)

attributes = listAttributes(ensembl)
geneList <- getBM(attributes=c("ensembl_gene_id", "external_gene_id", "entrezgene", "chromosome_name", "start_position", "end_position", "percentage_gc_content"), filters = c("chromosome_name", "start", "end"), values = c("10", "112950001", "135435000"),mart = ensembl)
geneListAll <- geneList
save(geneListAll, file="TCGAllg.Rdata")


library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", ensembl)
#attributes <- c("ensembl_gene_id","start_position","end_position","strand","hgnc_symbol","chromosome_name","entrezgene","ucsc","band")
attributes <- c("ensembl_gene_id", "external_gene_id", "entrezgene", "chromosome_name", "start_position", "end_position", "percentage_gc_content")
filters <- c("chromosome_name","start","end")
values <- list(chromosome="10",start="112950001",end="135435000")
geneList <- getBM(attributes=attributes, filters=filters, values=values, mart=ensembl)
geneListAll <- geneList
save(geneListAll, file="TCGAllg.Rdata")

# remove all gene without entrezID
geneList <- geneList[!is.na(geneList[,3]),]

# remove all genes not mapping to chr 1 to 22 or X and Y
geneList <- geneList[which(geneList[,4] %in% c(1:22, "X", "Y")),]
geneList <- geneList[which(geneList[,4] %in% c(10)),]
geneList <- geneList[,-1]
geneList <- geneList[,-2]
geneList <- unique(geneList)
geneList <- geneList[order(geneList[,2], geneList[,3], geneList[,4]),]

##########################################################################################################################################

# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

# Get list of cancer studies at server
getCancerStudies(mycgds)

# Get available case lists (collection of samples) for a given cancer study
# 21, 22, 23
mycancerstudy = getCancerStudies(mycgds)[65,1]
mycaselist = getCaseLists(mycgds, mycancerstudy)[1,1]

# Get available genetic profiles
mrnaProf = getGeneticProfiles(mycgds,mycancerstudy)[c(3),1]
methProf = getGeneticProfiles(mycgds,mycancerstudy)[c(6),1]

mutProf = getGeneticProfiles(mycgds,mycancerstudy)[c(9),1]
cnProf = getGeneticProfiles(mycgds,mycancerstudy)[c(6),1]

##########################################################################################################################################

assign("last.warning", NULL, envir = baseenv())

# Get data slices for a specified list of genes, genetic profile and case list
setwd("/home/wessel/Consultation/Ylstra_Bauke")
meData <- numeric()
cnData <- numeric()
geData <- numeric()
geneInfo <- numeric()
for (j in 1:nrow(geneList)){
	if (length(warnings()) == 0){
		print((c(j, nrow(geneList), dim(geData))))
		geneName <- as.character(geneList[j,1])
		geneData <- getProfileData(mycgds, geneName, c(cnProf, methProf, mrnaProf), mycaselist)
		if (dim(geneData)[2] == 3 & dim(geneData)[1] > 0){
			for (k in 1:3){
				if (class(geneData[,k]) == "factor"){ geneData[,k] <- as.numeric(levels(geneData[,k])[geneData[,k]]) }
				if (class(geneData[,k]) == "character"){ geneData[,k] <- as.numeric(geneData[,k]) }
			}
			if (any(apply(is.na(geneData), 2, sum) == nrow(geneData)) | sum(apply(is.na(geneData), 1, sum) == 0) < 2){ 
			} else {
				cnData <- rbind(cnData, geneData[,1])
				meData <- rbind(meData, geneData[,2])
				geData <- rbind(geData, geneData[,3])
				geneInfo <- rbind(geneInfo, geneList[j,])
				save(geneListAll, geneInfo, cnData, meData, geData, file="TCGAllg.Rdata")
			}
		}
	}
}
colnames(cnData) <- rownames(geneData)
colnames(meData) <- rownames(geneData)
colnames(geData) <- rownames(geneData)

rownames(cnData) <- geneInfo[,1]
rownames(meData) <- geneInfo[,1]
rownames(geData) <- geneInfo[,1]

clinInfo <- getClinicalData(mycgds, mycaselist)
save(clinInfo, geneListAll, geneInfo, cnData, meData, geData, file="TCGAllg.Rdata")


##########################################################################################################################################
# to do:
# - add mutation data
# - clean-up. 




