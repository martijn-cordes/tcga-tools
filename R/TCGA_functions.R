makeBins <- function(bin_size, build) {
  #load chromosome positions
  if (build == "HG18") {
    chromInfo <- read.delim("../ChromInfoHg18.txt", sep="", header=F)
    colnames(chromInfo) <- c("chromosome","position")
  }
  if (build == "HG19") {
    chromInfo <- read.delim("../ChromInfoHG19.txt", sep="", header=F)
    colnames(chromInfo) <- c("chromosome","position")
  }
  
  #Divide chromosomes in bins
  #10kb = 10000
  #1MB  = 1000000
  #10MB = 10000000        
  bins <- data.frame()
  for (i in 1:22) {
    chrom_length <- chromInfo$position[chromInfo$chromosome == paste("chr",i,sep="")]
    
    start <- seq(1,chrom_length-bin_size, bin_size)
    end <- start+(bin_size-1)
    bins <- rbind(bins, data.frame(chr=paste("chr",i,sep=""), start=start, end=end))
  }
  #X chromosome
  chrom_length <- chromInfo$position[chromInfo$chromosome == paste("chr","X",sep="")]
  
  start <- seq(1,chrom_length-bin_size, bin_size)
  end <- start+(bin_size-1)
  bins <- rbind(bins, data.frame(chr=paste("chr","X",sep=""), start=start, end=end))
  
  #Y chromosome
  chrom_length <- chromInfo$position[chromInfo$chromosome == paste("chr","Y",sep="")]
  
  start <- seq(1,chrom_length-bin_size, bin_size)
  end <- start+(bin_size-1)
  bins <- rbind(bins, data.frame(chr=paste("chr","Y",sep=""), start=start, end=end))
  
  bins
}

#Segmentcalling function
SegmentCalling <- function (copynumber_samples, tumor_segments,tumor_percentage) {
  tumor_perc <- tumor_percentage/100
  loss_cutoff <- log2((2-tumor_perc)/2)
  gain_cutoff <- log2((2+tumor_perc)/2)
  
  copynumbers_called <- data.frame(sample=1:nrow(bins))
  for (sample in 1:length(copynumber_samples)) {
    print (paste("performing calling for sample ", copynumber_samples[sample], " (",sample, "/", length(copynumber_samples),")", sep=""))
    
    temp_dataset <- tumor_segments[which(tumor_segments$Name == copynumber_samples[sample]),]
    
    #Range search
    #search for segments with bin regions
    #ranges of the bins
    bin_range <- with(bins,GRanges(chr,IRanges(bins$start,bins$end)))
    
    #ranges of segements of the sample
    sample_range <- with(temp_dataset,GRanges(Chromosome,IRanges(temp_dataset$Start,temp_dataset$End)))
    
    #overlap of segements in bins
    ov <- findOverlaps(sample_range, bin_range)
    
    #the overlapping hits in a dataframe
    hits <- data.frame(segment=queryHits(ov), bin_row=subjectHits(ov))
    
    #make dataframe for called regions
    bin_frame <- data.frame(sample=1:nrow(bins))
    bin_frame$sample <- 0
    
    hits$log2 <- temp_dataset$log2[hits$segment]
    
    #filter hits on log2 cutoff threshold
    #delete losses above threshold 
    loss_filtered <- hits$log2 < loss_cutoff
    bin_frame[hits$bin_row[loss_filtered],] <- -1
    
    #delete gains below threshold
    gain_filtered <- hits$log2 > gain_cutoff
    bin_frame[hits$bin_row[gain_filtered],] <- 1
    
    copynumbers_called <- cbind(copynumbers_called, sample=bin_frame[,1])
    
  }
  
  copynumbers_called <- copynumbers_called[,-1]
  colnames(copynumbers_called) <- copynumber_samples[1:length(copynumber_samples)]
  copynumbers_called
}

GeneCalling <- function(input_genes, coypnumbers_called, percent) {
  for (gene in 1:length(input_genes)) {
    
    print(paste("Scoring copynumber samples with an overlap of", input_genes[gene] ,"region with", percent, "percent"))         
    
    #start of copynumber region
    start_cn_region <- gene_list$chromStart[gene_list$gene == input_genes[gene]]
    end_cn_region <- gene_list$chromEnd[gene_list$gene == input_genes[gene]]
    
    #calculate bin start and end to copynumber region provided
    bin_start <- which(bins$chr == paste("chr",gene_list$chrom[gene_list$gene == input_genes[gene]], sep="") 
                       & bins$start == round(start_cn_region,-4) + 1)
    bin_end <- which(bins$chr == paste("chr",gene_list$chrom[gene_list$gene == input_genes[gene]], sep="") 
                     & bins$end == round(end_cn_region,-4))
    #bin_end <- bin_start +1
    criterium <- ((bin_end-bin_start)/100)*percent
    
    no_deletion <- c()
    deletion <- c()
    #check for region if all bins are called
    for (i in 1:ncol(copynumbers_called)) {
      calling_check <- sum(copynumbers_called[bin_start:(bin_end),i])
      
      if (calling_check == 0) { no_deletion <- c(no_deletion,colnames(copynumbers_called)[i])}
      if (calling_check > 0) { no_deletion <- c(no_deletion,colnames(copynumbers_called)[i])}
      if (calling_check < 0) { 
        calling_check <- calling_check*-1    
        if (calling_check > criterium) {
          deletion <- c(deletion,colnames(copynumbers_called)[i])
        } else { no_deletion <- c(no_deletion,colnames(copynumbers_called)[i]) }
      }
    }    
    
    #add clinical parameter (deletion or no deletion) to clinical data
    copynumber_clinical_data[,input_genes[gene]] <- 0
    for (j in 1:length(deletion)) {
      copynumber_clinical_data[,input_genes[gene]][which(copynumber_clinical_data$sampleID == deletion[j])] <- paste("a_",input_genes[gene], "_deletion", sep="")  
    }
    
    for (j in 1:length(no_deletion)) {
      copynumber_clinical_data[,input_genes[gene]][which(copynumber_clinical_data$sampleID == no_deletion[j])] <- paste("no_",input_genes[gene], "_deletion", sep="")     
    }

  }
  copynumber_clinical_data
}


#Kaplan Meier
make_KM <- function(copynumber_clinical_data, gene){
  column_number <- which(colnames(copynumber_clinical_data)==input_genes[gene])
  copynumber_clinical_data$SurvObj <- with(copynumber_clinical_data, Surv(as.numeric(as.vector(copynumber_clinical_data[,which(colnames(copynumber_clinical_data)=="_OS")])) ,as.numeric(as.vector(copynumber_clinical_data[,which(colnames(copynumber_clinical_data)=="_EVENT")]) == 1)))
  km.by.gene <- survfit(SurvObj ~ copynumber_clinical_data[,column_number], data = copynumber_clinical_data, conf.type = "log-log")
  rate <- (km.by.gene$n[1]/sum(km.by.gene$n))*100
  
  #more than one group to test
  if (length(km.by.gene$n) > 1) {
    pval <- survdiff(SurvObj ~ copynumber_clinical_data[,column_number], data = copynumber_clinical_data)
    pval.by.gene <- 1-pchisq(pval$chisq,length(pval$n)-1)
    
    #pvalue below 0.05 and deletion rate above 5%
    if (pval.by.gene <= 1 || rate > 5) {
      print(paste("Creating KM for",input_genes[gene]))
      #copynumber_info <- rbind(copynumber_info, c(input_genes[gene],pval.by.gene, km.by.gene$n[1],km.by.gene$n[2]))
      km_info <- c(input_genes[gene],pval.by.gene, km.by.gene$n[1],km.by.gene$n[2])
      
      #Make KM plot 
      png(paste("plots/",input_genes[gene],".png",sep=""), width=297, height=210, units='mm', res=150)
      plot(km.by.gene, mark.time=T, xlab='Days', ylab='Survival', col=c("dodgerblue4","gold2"),lwd=2, lty=1:2, yaxt="n", main=paste(input_genes[gene], " in TCGA Grade 2 and 3 LGG  (n = ", length(copynumber_samples), ")",sep=""))
      #plot(km.by.gene, mark.time=T, xlab='Days', ylab='Survival', col=c("dodgerblue4","gold2"),lwd=2, lty=1:2, yaxt="n", main=paste(input_genes[gene], " in TCGA Stage IV Colon adenocarcinoma's  (n = ",length(copynumber_samples), ")",sep=""))
      axis(1, at = c(500,1500,2500,3500,4500,5500,6500))
      axis(2, at = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), las=1)
      legend("topright", c(paste(input_genes[gene],"deletion","(n =",length(which(copynumber_clinical_data[,input_genes[gene]]==paste("a_",input_genes[gene], "_deletion", sep=""))),")"),paste("no", input_genes[gene], "deletion","(n =",length(which(copynumber_clinical_data[,input_genes[gene]]==paste("no_",input_genes[gene],"_deletion",sep=""))),")")), col=c("dodgerblue4","gold2"), lwd=2)
      #legend("topright", c(paste(input_genes[gene],"(n =",length(which(copynumber_clinical_data[,input_genes[gene]]==paste("a_SV_deletion", sep=""))),")"),paste("no ", input_genes[gene],"(n =",length(which(copynumber_clinical_data[,input_genes[gene]]==paste("no_SV",sep=""))),")")), col=c("dodgerblue4","gold2"), lwd=2)
      legend("bottomright",inset = 0.05, paste("pvalue =",pval.by.gene),box.col=NA )
      mtext(paste(percent,"% overlap with deletion regions and log2 cutoff for tumor percentage of ",tumor_percentage,"%",sep=""))
      
      dev.off()
      
    }
  } 
  km_info
}


KM_per_gene <- function(input_genes, copynumber_clinical_data){
  #KM loop
  km_gene_overview <- data.frame()
  for (gene in 1:length(input_genes)) {
    km_info <- make_KM(copynumber_clinical_data, gene)
    km_gene_overview <- rbind(km_gene_overview, km_info)
    colnames(km_gene_overview) <- c("Gene","Pvalue","Deletion","No_Deletion")
  }
  km_gene_overview
}

#make cghRegions object
make_cghRegions <- function(copynumbers_called) {
  print("Creating cghRegions object")  
  assay_matrix <- as.matrix(copynumbers_called,rownames=T)
  regions_assaydata <- assayDataNew(regions=assay_matrix)
  
  regions_dataframe <- data.frame(Chromosome=as.numeric(substr(bins$chr, 4,5)), Start=as.numeric(bins$start), End=as.numeric(bins$end), Nclone=0, AveDist=0)
  featurdata_object <- new("AnnotatedDataFrame", data=regions_dataframe)
  
  regions_object <- new("cghRegions", assayData=regions_assaydata, featureData=featurdata_object)
  regions_object
}

#Frequency plot
setMethod("frequencyPlot", signature(x="cghRegions", y="missing"),
          function (x, y, main='Frequency Plot', gaincol='blue', losscol='red', misscol=NA, build=build,... )
          {
            if (build=="HG18") { genome_build='GRCh36'}
            if (build=="HG19") { genome_build='GRCh37'}

            chrom <- chromosomes(x)
            pos <- bpstart(x)
            pos2 <- bpend(x)
            uni.chrom <- unique(chrom)
            chrom.lengths <- CGHbase:::.getChromosomeLengths(genome_build)[as.character(uni.chrom)]
            chrom.ends <- integer()
            cumul <- 0
            for (j in uni.chrom) {
              pos[chrom > j] <- pos[chrom > j] + chrom.lengths[as.character(j)]
              pos2[chrom > j] <- pos2[chrom > j] + chrom.lengths[as.character(j)]
              cumul <- cumul + chrom.lengths[as.character(j)]
              chrom.ends <- c(chrom.ends, cumul)
            }
            names(chrom.ends) <- names(chrom.lengths)
            calls <- regions(x)
            loss.freq <- rowMeans(calls < 0)
            gain.freq <- rowMeans(calls > 0)
            plot(NA, xlim=c(0, max(pos2)), ylim=c(-1,1), type='n', xlab='chromosomes', ylab='frequency', xaxs='i', xaxt='n', yaxs='i', yaxt='n', main=main,...)
            if (!is.na(misscol)) {
              rect(0, -1, max(pos2), 1, col=misscol, border=NA)
              rect(pos, -1, pos2, 1, col='white', border=NA)
            }
            rect(pos, 0, pos2, gain.freq, col=gaincol, border=gaincol)
            rect(pos, 0, pos2, -loss.freq, col=losscol, border=losscol)
            box()
            abline(h=0)
            if (length(chrom.ends) > 1)
              for (j in names(chrom.ends)[-length(chrom.ends)])
                abline(v=chrom.ends[j], lty='dashed')
            ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
            axis(side=1,at=ax,labels=uni.chrom,cex=.2,lwd=.5,las=1,cex.axis=1,cex.lab=1)
            axis(side=2, at=c(-1, -0.5, 0, 0.5, 1), labels=c('100 %', ' 50 %', '0 %', '50 %', '100 %'), las=1)
            mtext('gains', side=2, line=3, at=0.5)
            mtext('losses', side=2, line=3, at=-0.5)
            ### number of data points
            mtext(paste(nrow(x), 'regions'), side=3, line=0, adj=0)
          })



#Correlation plot
library(RColorBrewer)
cor.mat.by.region <- function(x, m, build, main='Correlation Matrix', ...) {
  if (build=="HG18") { genome_build='GRCh36'}
  if (build=="HG19") { genome_build='GRCh37'}
  
  image(x=1:(nrow(m)+1), y=1:(nrow(m)+1), z=m, zlim=c(-1,1), col=brewer.pal(10, 'RdBu'), xlab=NA, ylab=NA, useRaster=TRUE, main=main, xaxt='n', yaxt='n', ...)
  
  chr.change <- chromosomes(x) != c(chromosomes(x)[-1], 0)
  chr.change <- which(chr.change) + 1
  abline(h=chr.change[-length(chr.change)], v=chr.change[-length(chr.change)], lty='dashed')
  
  if (nrow(m) < 100) {
    a <- add.cytobands(fData(x), genome.build=genome_build)
    axis(side=1, at=1:nrow(x)+.5, labels=a$cytoband, las=2, line=-.5, tick=FALSE, cex.axis=.4)
    axis(side=2, at=1:nrow(x)+.5, labels=a$cytoband, las=2, line=-.5, tick=FALSE, cex.axis=.1)
    axis(side=4, at=1:nrow(x)+.5, labels=a$cytoband, las=2, line=-.5, tick=FALSE, cex.axis=.4)
  } else {
    ax <- (chr.change + c(0, chr.change[-length(chr.change)]))/2
    axis(side=1, at=ax, labels=unique(chromosomes(x)), lwd=.5, las=1, cex.axis=1)
    axis(side=2, at=ax, labels=unique(chromosomes(x)), lwd=.5, las=1, cex.axis=.5)
    axis(side=4, at=ax, labels=unique(chromosomes(x)), lwd=.5, las=1, cex.axis=.5)
  }
}

cor.mat.by.size <- function(x, m, build, main='Correlation Matrix', ...) {
  if (build=="HG18") { genome_build='GRCh36'}
  if (build=="HG19") { genome_build='GRCh37'}
  
  chrom <- chromosomes(x)
  pos <- bpstart(x)
  pos2 <- bpend(x)
  uni.chrom <- unique(chrom)
  chrom.lengths <- CGHbase:::.getChromosomeLengths(genome_build)[as.character(uni.chrom)]
  chrom.ends <- integer()
  cumul <- 0
  for (j in uni.chrom) {
    pos[chrom > j] <- pos[chrom > j] + chrom.lengths[as.character(j)]
    pos2[chrom > j] <- pos2[chrom > j] + chrom.lengths[as.character(j)]
    cumul <- cumul + chrom.lengths[as.character(j)]
    chrom.ends <- c(chrom.ends, cumul)
  }
  names(chrom.ends) <- names(chrom.lengths)
  
  image(x=c(pos, pos2[length(pos2)]), y=c(pos, pos2[length(pos2)]), z=m, zlim=c(-1,1), col=brewer.pal(10, 'RdBu'), xlab=NA, ylab=NA, main=main, xaxt='n', yaxt='n', ...)
  
  abline(h=chrom.ends[-length(chrom.ends)], v=chrom.ends[-length(chrom.ends)], lty='dashed')
  
  ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
  axis(side=1, at=ax, labels=uni.chrom, lwd=.5, las=1, cex.axis=1)
  axis(side=2, at=ax, labels=uni.chrom, lwd=.5, las=1, cex.axis=.5)
  axis(side=4, at=ax, labels=uni.chrom, lwd=.5, las=1, cex.axis=.5)
}

PrintSegments <- function(plot_samples,build) {
  #load chromosome positions
  if (build == "HG18") {
    chromInfo <- read.delim("ChromInfoHg18.txt", sep="", header=F)
    colnames(chromInfo) <- c("chromosome","position")
    genome_build='GRCh36'
  }
  if (build == "HG19") {
    chromInfo <- read.delim("ChromInfoHG19.txt", sep="", header=F)
    colnames(chromInfo) <- c("chromosome","position")
    genome_build='GRCh37'
  }
  
  SVsamples <- copynumber_samples[plot_samples]
  
  #Chromosome lengths ordered
  chromInfo$position<-as.numeric(chromInfo$position)
  
  chr.list <- paste("chr", c(1:22, "X", "Y"), sep="")
  data <- chromInfo[chromInfo[,1] %in% chr.list,]
  
  substr(data[,1], 4,5) -> cn
  cn[cn == "X"] <- 100
  cn[cn == "Y"] <- 101
  order(as.numeric(cn), as.numeric(data[,2])) -> ord
  data[ord,] -> data2
  cno <- cn[ord]
  
  cnou <- sort(as.numeric(unique(cno)))
  sapply(cnou, function(x) { max(as.numeric(as.vector(data2[cno == x,2])))  }) -> cmax
  clin <- c(0, cumsum(cmax))
  data2$conlin <- NA
  sapply(1:length(cnou), function(x) { as.numeric(data2[cno == cnou[x],2]) + clin[x] -> data2$conlin[cno == cnou[x]]}) -> res
  data2$conlin <- unlist(res)
  
  
  #Plot segment plots of certain samples
  for (i in 1:length(SVsamples)) {
    tumor_dataset <- tumor_segments_all[which(tumor_segments_all$Name == SVsamples[i]),]
    #tumor_dataset <- input_frame
    
    #make tumor_object to plot
    input_frame <- tumor_segments_all[which(tumor_segments_all$Name == SVsamples[i]),]
    colnames(input_frame)[5] <- SVsamples[i]
    input_frame$Name <- paste(input_frame$Start,"-",input_frame$End, sep="")
    if (length(which(duplicated(input_frame$Name))) > 0){
      input_frame <- input_frame[-which(duplicated(input_frame$Name)),]
    }
    tumor_object <- make_cghRaw(input_frame)
    
    png(paste("plots/",SVsamples[i],".png",sep=""), width=297, height=210, units='mm', res=150)
    plot(tumor_object, dotres=1, ylim=c(-5,5), main=samples[i], build=genome_build,dlcol="red")
    #plot(tumor_object[chromosomes(tumor_object)==14], dotres=1, ylim=c(-5,5), main=samples[i], build=genome_build,dlcol="red")
    
    for (i in 1:22) {
      
      start <- tumor_dataset$Start[which(tumor_dataset$Chromosome == paste("chr",i,sep=""))]
      end <- tumor_dataset$End[which(tumor_dataset$Chromosome == paste("chr",i,sep=""))]
      seg <- tumor_dataset$log2[which(tumor_dataset$Chromosome == paste("chr",i,sep=""))]
      
      if (i > 1) {
        
        start <- start + data2$conlin[i-1]
        end <- end + data2$conlin[i-1]
        
        
        for (xi in 1:length(seg)) {
          lines(c(start[xi],end[xi]),c(seg[xi],seg[xi]), col="red",lwd=2)
          #lines(c(2300791338,2401130253),c(0,0), col="red",lwd=2)
          
        }  
      }
      
      else {		
        for (xi in 1:length(seg)) { lines(c(start[xi],end[xi]),c(seg[xi],seg[xi]), col="red", lwd=2) }
      }	
    }
    
    dev.off()
    
  }
  
  
}




#SegmentFrame function
SegmentFrame <- function (copynumber_samples, tumor_segments) {
  copynumbers_called <- data.frame(sample=1:nrow(bins))
  for (sample in 1:length(copynumber_samples)) {
    print (paste("performing calling for sample ", copynumber_samples[sample], " (",sample, "/", length(copynumber_samples),")", sep=""))
    
    temp_dataset <- tumor_segments[which(tumor_segments$Name == copynumber_samples[sample]),]
    
    #Range search
    #search for segments with bin regions
    #ranges of the bins
    bin_range <- with(bins,GRanges(chr,IRanges(bins$start,bins$end)))
    
    #ranges of segements of the sample
    sample_range <- with(temp_dataset,GRanges(Chromosome,IRanges(temp_dataset$Start,temp_dataset$End)))
    
    #overlap of segements in bins
    ov <- findOverlaps(sample_range, bin_range)
    
    #the overlapping hits in a dataframe
    hits <- data.frame(segment=queryHits(ov), bin_row=subjectHits(ov))
    
    #make dataframe for called regions
    bin_frame <- data.frame(sample=1:nrow(bins))
    bin_frame$sample <- NA
    
    hits$log2 <- temp_dataset$log2[hits$segment]
    
    #filter hits on log2 cutoff threshold
    #delete losses above threshold 
    loss_filtered <- hits$log2 < 0
    bin_frame[hits$bin_row[loss_filtered],] <- hits$log2[loss_filtered]
    
    #delete gains below threshold
    gain_filtered <- hits$log2 > 0
    bin_frame[hits$bin_row[gain_filtered],] <- hits$log2[gain_filtered]
    
    copynumbers_called <- cbind(copynumbers_called, sample=bin_frame[,1])
    
  }
  
  if (ncol(copynumbers_called) > 2) {
  copynumbers_called <- copynumbers_called[,-1]
  colnames(copynumbers_called) <- copynumber_samples[1:length(copynumber_samples)]
  }
  copynumbers_called
}



GeneToCopynumberBiomart <- function(input_genes, coypnumbers_frame, copynumber_samples) {
  print(paste("Getting copynumber values of ", length(copynumber_samples), " samples for ", length(input_genes) ," genes (running time: 30 min)",sep=""))
  
  #layout of gene frame
  gene_frame <- data.frame(matrix(ncol = length(copynumber_samples)+1, nrow = length(input_genes)))
  gene_frame[,1] <- input_genes
  colnames(gene_frame)[1] <- "genes"
  colnames(gene_frame)[2:(length(copynumber_samples)+1)] <- copynumber_samples
  

  start_cn_region <- sapply(input_genes, function(x) { gene_list$chromStart[gene_list$gene == x]})
  #start_cn_region <- unlist(start_cn_region)
    
  end_cn_region <- sapply(input_genes,function(x) {  gene_list$chromEnd[gene_list$gene == x]})
  #end_cn_region <- unlist(end_cn_region)
  
  print ("Located start and stop positions of all genes, getting copynumber bins for gene positions")
  
  bin_start <- sapply(start_cn_region, function(x) {which(bins$chr == gene_list$chrom[gene_list$chromStart == x] & bins$start == round(as.numeric(x),-4) + 1)} )
  
  #Correct for the same end positions
  for (i in 1:length(bin_start)) {
    if(length(bin_start[[i]])>1) {
      print(bin_start[[i]])
      bin_start[i] <- which(bins$chr == input_genes[i] 
                            & bins$start == round(start_cn_region[i],-4))+1
    }
  }
  
  bin_end <- sapply(end_cn_region, function(x) {which(bins$chr == gene_list$chrom[gene_list$chromEnd == x] & bins$end == round(x,-4))+1})
  #Correct for the same end positions
  for (i in 1:length(bin_end)) {
    if(length(bin_end[[i]])>1) {
      print(bin_end[[i]])
      bin_end[i] <- which(bins$chr == input_genes[i]
                          & bins$end == round(end_cn_region[i],-4))+1
    }
  }
  
  bin_start <- as.vector(unlist(bin_start))
  bin_end <- as.vector(unlist(bin_end))
  
  #fill gene frame
  for (sample in 2:(length(copynumber_samples)+1)) {
    #get NAs
    nas <- which(is.na(copynumbers_frame[bin_start,sample] == copynumbers_frame[bin_end,sample]))
    gene_frame[nas,sample] <- NA
    #if copynumber level is same at beginning and end of gene position then put in copynumber leel
    equal <- which(copynumbers_frame[bin_start,sample] == copynumbers_frame[bin_end,sample])
    gene_frame[equal,sample] <- copynumbers_frame[equal,sample]
    #if copynumber differs from beginning of gene to end, breakpoint in gene
    not_equal <- which(copynumbers_frame[bin_start,sample] != copynumbers_frame[bin_end,sample])
    gene_frame[not_equal,sample] <- "breakpoint in gene"
  }
  
  gene_frame
}



GeneToCopynumber <- function(gene_list, copynumbers_frame, copynumber_samples) {
  print(paste("Getting copynumber values of ", length(copynumber_samples), " samples for ", length(gene_list$gene) ," genes (running time: 30 min)",sep=""))
  
  #layout of gene frame
  gene_frame <- data.frame(matrix(ncol = length(copynumber_samples)+1, nrow = length(gene_list$gene)))
  colnames(gene_frame) <- copynumber_samples
  
  gene_frame <- gene_list
  gene_frame[7:(length(copynumber_samples)+6)] <- NA
  colnames(gene_frame)[7:(length(copynumber_samples)+6)] <- copynumber_samples
    
  start_cn_region <- gene_list$chromStart
  end_cn_region <- gene_list$chromEnd
  
  print ("Located start and stop positions of all genes, getting copynumber bins for gene positions")
  
  #bin_start <- sapply(start_cn_region, function(x) {which(bins$chr == gene_list$chrom[gene_list$chromStart == x] & bins$start == round(as.numeric(x),-4) + 1)} )
  #bin_end <- sapply(end_cn_region, function(x) {which(bins$chr == gene_list$chrom[gene_list$chromEnd == x] & bins$end == round(x,-4))+1})
  
  bin_size <- bins$end[1]
  bin_start <- sapply(start_cn_region, function(x) {  as.numeric(rownames(bins[which(bins$chr == gene_list$chrom[gene_list$chromStart == x])[round(as.numeric(x)/bin_size) +1],]))} )
  bin_end <- sapply(end_cn_region, function(x) {as.numeric(rownames(bins[which(bins$chr == gene_list$chrom[gene_list$chromEnd == x])[round(as.numeric(x)/bin_size) +1],]))})
  
  bin_start <- unlist(bin_start)
  bin_end <- unlist(bin_end)
  print ("Located copynumber bins for all gene positions, getting copynumber values for all samples")
  
#fill gene frame
for (sample in 1:length(copynumber_samples)) { 
  #get NAs
  nas <- which(is.na(copynumbers_frame[bin_start,sample] == copynumbers_frame[bin_end,sample]))
  gene_frame[nas,(sample+6)] <- NA
  #if copynumber level is same at beginning and end of gene position then put in copynumber leel
  equal <- which(copynumbers_frame[bin_start,sample] == copynumbers_frame[bin_end,sample])
  gene_frame[equal,(sample+6)] <- copynumbers_frame[bin_start[equal],sample]
  #if copynumber differs from beginning of gene to end, breakpoint in gene
  not_equal <- which(copynumbers_frame[bin_start,sample] != copynumbers_frame[bin_end,sample])
  gene_frame[not_equal,(sample+6)] <- "breakpoint in gene"
}

gene_frame
}






