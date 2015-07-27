#extract sample sizes and make one clinical database based on all overlapping samples
cn_samples <- length(as.vector(copynumber_clinical_data$sampleID))
mut_samples <- length(as.vector(mut_clinical_data$sampleID))
exp_samples <- length(as.vector(expression_clinical_data$sampleID))

sizes <- c(cn_samples,mut_samples,exp_samples)
#sizes <- c(cn_samples,mut_samples)
smallest_size <- sort(sizes)

if (cn_samples == smallest_size[1]){
  cn <- copynumber_clinical_data[which(copynumber_clinical_data$sampleID %in% copynumber_clinical_data$sampleID == TRUE),]
  mut <- mut_clinical_data[which(mut_clinical_data$sampleID %in% copynumber_clinical_data$sampleID == TRUE),]
  exp <- expression_clinical_data[which(expression_clinical_data$sampleID %in% copynumber_clinical_data$sampleID == TRUE),]
  all_clinical_data <- copynumber_clinical_data[which(copynumber_clinical_data$sampleID %in% copynumber_clinical_data$sampleID == TRUE),]
}

if (mut_samples == smallest_size[1]){
  cn <- copynumber_clinical_data[which(copynumber_clinical_data$sampleID %in% mut_clinical_data$sampleID == TRUE),]
  mut <- mut_clinical_data[which(mut_clinical_data$sampleID %in% mut_clinical_data$sampleID == TRUE),]
  exp <- expression_clinical_data[which(expression_clinical_data$sampleID %in% mut_clinical_data$sampleID == TRUE),]
  #all_clinical_data <- mut_clinical_data[which(mut_clinical_data$sampleID %in% mut_clinical_data$sampleID == TRUE),]
  
  #TEMP FIND SOLUTION, CN IS SMALLER
  cn <- copynumber_clinical_data[which(copynumber_clinical_data$sampleID %in% mut$sampleID == TRUE),]
  mut <- mut_clinical_data[which(mut_clinical_data$sampleID %in% cn$sampleID == TRUE),]
  exp <- expression_clinical_data[which(expression_clinical_data$sampleID %in% cn$sampleID == TRUE),]
  
  
  
  all_clinical_data <- mut_clinical_data[which(mut_clinical_data$sampleID %in% mut$sampleID == TRUE),]
}

if (exp_samples == smallest_size[1]){
  cn <- copynumber_clinical_data[which(copynumber_clinical_data$sampleID %in% expression_clinical_data$sampleID == TRUE),]
  mut <- mut_clinical_data[which(mut_clinical_data$sampleID %in% expression_clinical_data$sampleID == TRUE),]
  exp <- expression_clinical_data[which(expression_clinical_data$sampleID %in% expression_clinical_data$sampleID == TRUE),]
  all_clinical_data <- expression_clinical_data[which(expression_clinical_data$sampleID %in% expression_clinical_data$sampleID == TRUE),]
}



#filter for genes which are not in all sets
cn_genes <- colnames(cn)[132:length(colnames(cn))]
mut_genes <- colnames(mut)[132:length(colnames(mut))]
if(length(which(input_genes%in%mut_genes != TRUE))>0){
  print (paste("gene", input_genes[which(input_genes%in%mut_genes != TRUE)], "not found"))
  input_genes <- input_genes[-which(input_genes%in%mut_genes != TRUE)]
} else {print("all genes found, extracting status")}

for (gene in 1:length(input_genes)) {
  
  #all_clinical_data <- all_clinical_data[,-which(colnames(all_clinical_data)==input_genes[gene])]
  
  #add copynumber data
  all_clinical_data[,paste(input_genes[gene],"_copynumber", sep="")] <- cn[,which(colnames(cn)==input_genes[gene])]
  #add mutation data
  all_clinical_data[,paste(input_genes[gene],"_mutation", sep="")] <- mut[,which(colnames(mut)==input_genes[gene])]
  #add expression data
  #all_clinical_data[,paste(input_genes[gene],"_expression", sep="")] <- exp[,which(colnames(exp)==input_genes[gene])]
  
  
}

#####WORKS FOR NOW
##Take largest set by hand and combine

#Combine the mutation an deletion set
all_clinical_data <- copynumber_clinical_data
for (gene in 1:length(input_genes)) {
  mutation <- as.vector(mut_clinical_data$sampleID[which(mut_clinical_data[,input_genes[gene]] == paste("a_",input_genes[gene],"_mutation", sep=""))])
  
  all_clinical_data[,paste(input_genes[gene],"_mutation",sep="")] <- paste("no_",input_genes[gene],"_mutation", sep="")
  for (j in 1:length(mutation)) {
    all_clinical_data[,paste(input_genes[gene],"_mutation",sep="")][which(copynumber_clinical_data$sampleID == mutation[j])] <- paste("a_",input_genes[gene], "_mutation", sep="")  
  }
  all_clinical_data[,paste(input_genes[gene],"_copynumber",sep="")] <- copynumber_clinical_data[,paste(input_genes[gene],sep="")]

}
#Create event and check if mutation or deletion happens

event_info <- data.frame(Gene=0,Pvalue=0, Event=0, No_Event=0)
for (gene in 1:length(input_genes)) {  
  all_clinical_data[,paste(input_genes[gene],"_event", sep="")] <- paste("no_",input_genes[gene],"_event",sep="")
  for (sample in 1:nrow(all_clinical_data)) { 
    #deletion and mutation
    if ((all_clinical_data[sample,paste(input_genes[gene],"_copynumber", sep="")] == paste("a_",input_genes[gene],"_deletion",sep="") &&
         all_clinical_data[sample,paste(input_genes[gene],"_mutation", sep="")] == paste("a_",input_genes[gene],"_mutation",sep=""))==TRUE)
      {
           all_clinical_data[sample,paste(input_genes[gene],"_event", sep="")] <- paste("a_",input_genes[gene],"_event",sep="")
      }
    #deletion and no mutation
    if ((all_clinical_data[sample,paste(input_genes[gene],"_copynumber", sep="")] == paste("a_",input_genes[gene],"_deletion",sep="") &&
           all_clinical_data[sample,paste(input_genes[gene],"_mutation", sep="")] == paste("no_",input_genes[gene],"_mutation",sep=""))==TRUE)
    {
      all_clinical_data[sample,paste(input_genes[gene],"_event", sep="")] <- paste("a_",input_genes[gene],"_event",sep="")
    }
    #no deletion and a mutation
    if ((all_clinical_data[sample,paste(input_genes[gene],"_copynumber", sep="")] == paste("no_",input_genes[gene],"_deletion",sep="") &&
           all_clinical_data[sample,paste(input_genes[gene],"_mutation", sep="")] == paste("a_",input_genes[gene],"_mutation",sep=""))==TRUE)
    {
      all_clinical_data[sample,paste(input_genes[gene],"_event", sep="")] <- paste("a_",input_genes[gene],"_event",sep="")
    }
   
  }

  #KM
  column_number <- which(colnames(all_clinical_data)==paste(input_genes[gene],"_event",sep=""))
  all_clinical_data$SurvObj <- with(all_clinical_data, Surv(as.numeric(as.vector(all_clinical_data[,which(colnames(all_clinical_data)=="_OS")])) ,as.numeric(as.vector(all_clinical_data[,which(colnames(all_clinical_data)=="_EVENT")]) == 1)))
  km.by.gene <- survfit(SurvObj ~ all_clinical_data[,column_number], data = all_clinical_data, conf.type = "log-log")
  rate <- (km.by.gene$n[1]/sum(km.by.gene$n))*100

  #more than one group to test
  if (length(km.by.gene$n) > 1) {
    pval <- survdiff(SurvObj ~ all_clinical_data[,column_number], data = all_clinical_data)
    pval.by.gene <- 1-pchisq(pval$chisq,length(pval$n)-1)
    
    #pvalue below 0.05 and deletion rate above 5%
    if (pval.by.gene <= 0.05 || rate > 5) {
      print(paste("Creating KM for",input_genes[gene]))
      event_info <- rbind(event_info, c(input_genes[gene],pval.by.gene, km.by.gene$n[1],km.by.gene$n[2]))
      
      #Make KM plot 
      png(paste("events/",input_genes[gene],"_colon.png",sep=""), width=297, height=210, units='mm', res=150)
      #plot(km.by.gene, mark.time=T, xlab='Days', ylab='Survival', col=c("dodgerblue4","gold2"),lwd=2, lty=1:2, yaxt="n", main=paste(input_genes[gene], " in TCGA Stage IV Colon adenocarcinoma's  (n = ", nrow(all_clinical_data), ")",sep=""))
      plot(km.by.gene, mark.time=T, xlab='Days', ylab='Survival', col=c("dodgerblue4","gold2"),lwd=2, lty=1:2, yaxt="n", main=paste(input_genes[gene], " in TCGA Colon adenocarcinoma's  (n = ", nrow(all_clinical_data), ")",sep=""))
      axis(1, at = c(500,1500,2500,3500,4500,5500,6500))
      axis(2, at = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), las=1)
      legend("topright", c(paste(input_genes[gene],"deletion or mutation","(n =",length(which(all_clinical_data[,column_number]==paste("a_",input_genes[gene], "_event", sep=""))),")"),paste("no", input_genes[gene], "deletion or mutation","(n =",length(which(all_clinical_data[,column_number]==paste("no_",input_genes[gene],"_event",sep=""))),")")), col=c("dodgerblue4","gold2"), lwd=2)
      legend("bottomright",inset = 0.05, paste("pvalue =",pval.by.gene),box.col=NA )
      dev.off()
      
    }
  }
    
}
event_info <- event_info[-1,]


#Compare events
comparison <- "and"
event_info <- data.frame(Gene=0,Pvalue=0, Event=0, No_Event=0)
all_clinical_data[,paste(input_genes[1],"_",input_genes[2],"_group_event", sep="")] <- paste("no_",input_genes[1],"_",input_genes[2],"_group_event",sep="")
  for (sample in 1:nrow(all_clinical_data)) { 
    #event and event
    if (comparison == "and"){
      if ((all_clinical_data[sample,paste(input_genes[1],"_event", sep="")] == paste("a_",input_genes[1],"_event",sep="") &&
             all_clinical_data[sample,paste(input_genes[2],"_event", sep="")] == paste("a_",input_genes[2],"_event",sep=""))==TRUE)
      {
        all_clinical_data[sample,paste(input_genes[1],"_",input_genes[2],"_group_event", sep="")] <- paste("a_",input_genes[1],"_",input_genes[2],"_group_event",sep="")
      }
    }
    
    if (comparison == "or"){
      #event and no event
      if ((all_clinical_data[sample,paste(input_genes[1],"_copynumber", sep="")] == paste("a_",input_genes[1],"_event",sep="") &&
             all_clinical_data[sample,paste(input_genes[2],"_mutation", sep="")] == paste("no_",input_genes[2],"_event",sep=""))==TRUE)
      {
        all_clinical_data[sample,paste(input_genes[1],"_",input_genes[2],"_group_event", sep="")] <- paste("a_",input_genes[1],"_",input_genes[2],"_group_event",sep="")
      }
      #no event and a event
      if ((all_clinical_data[sample,paste(input_genes[1],"_event", sep="")] == paste("no_",input_genes[1],"_event",sep="") &&
             all_clinical_data[sample,paste(input_genes[2],"_event", sep="")] == paste("a_",input_genes[2],"_event",sep=""))==TRUE)
      {
        all_clinical_data[sample,paste(input_genes[1],"_",input_genes[2],"_group_event", sep="")] <- paste("a_",input_genes[1],"_",input_genes[2],"_group_event",sep="")
      }
    }
  }
  
  #KM
  column_number <- which(colnames(all_clinical_data)==paste(input_genes[1],"_",input_genes[2],"_group_event", sep=""))
  all_clinical_data$SurvObj <- with(all_clinical_data, Surv(as.numeric(as.vector(all_clinical_data[,which(colnames(all_clinical_data)=="_OS")])) ,as.numeric(as.vector(all_clinical_data[,which(colnames(all_clinical_data)=="_EVENT")]) == 1)))
  km.by.gene <- survfit(SurvObj ~ all_clinical_data[,column_number], data = all_clinical_data, conf.type = "log-log")
  rate <- (km.by.gene$n[1]/sum(km.by.gene$n))*100
  
  #more than one group to test
  if (length(km.by.gene$n) > 1) {
    pval <- survdiff(SurvObj ~ all_clinical_data[,column_number], data = all_clinical_data)
    pval.by.gene <- 1-pchisq(pval$chisq,length(pval$n)-1)
    
    #pvalue below 0.05 and deletion rate above 5%
    if (pval.by.gene <= 0.05 || rate > 5) {
      print(paste("Creating KM for",paste(input_genes[1],"_",comparison,"_",input_genes[2],"_group_event", sep="")))
      event_info <- rbind(event_info, c(input_genes[gene],pval.by.gene, km.by.gene$n[1],km.by.gene$n[2]))
      
      #Make KM plot 
      png(paste("events/",input_genes[1],"_",comparison,"_",input_genes[2],"_colon.png",sep=""), width=297, height=210, units='mm', res=150)
      plot(km.by.gene, mark.time=T, xlab='Days', ylab='Survival', col=c("dodgerblue4","gold2"),lwd=2, lty=1:2, yaxt="n", main=paste(input_genes[1]," ",comparison," ",input_genes[2], " in TCGA Colon adenocarcinoma's  (n = ", nrow(all_clinical_data), ")",sep=""))
      #plot(km.by.gene, mark.time=T, xlab='Days', ylab='Survival', col=c("dodgerblue4","gold2"),lwd=2, lty=1:2, yaxt="n", main=paste(input_genes[1]," ",comparison," ",input_genes[2], " in TCGA Stage IV Colon adenocarcinoma's  (n = ", nrow(all_clinical_data), ")",sep=""))
      axis(1, at = c(500,1500,2500,3500,4500,5500,6500))
      axis(2, at = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), las=1)
      legend("topright", c(paste(input_genes[1],"deletion or mutation", comparison, input_genes[2],"deletion or mutation " ,"(n =",length(which(all_clinical_data[,column_number]==paste("a_",input_genes[1],"_",input_genes[2],"_group_event",sep=""))),")"),paste("no", input_genes[1],"and", input_genes[2], "deletion or mutation","(n =",length(which(all_clinical_data[,column_number]==paste("no_",input_genes[1],"_",input_genes[2],"_group_event",sep=""))),")")), col=c("dodgerblue4","gold2"), lwd=2)
      legend("bottomright",inset = 0.05, paste("pvalue =",pval.by.gene),box.col=NA )
      dev.off()
      
    }
  }

event_info <- event_info[-1,]



  
