#check for arid1a events
combined_clinical$arid1a_event <- "no_arid1a_event"
for (i in 1:length(combined_clinical$ARID1A)) {
  if ((combined_clinical$ARID1A[i] == "a_ARID1A_mutation" && combined_clinical$arid1a_del[i] == "arid1a_deletion")==TRUE){
    combined_clinical$arid1a_event[i]<-"arid1a_double_event"
  }
  if ((combined_clinical$ARID1A[i] == "a_ARID1A_mutation" && combined_clinical$arid1a_del[i] == "no_arid1a_deletion")==TRUE){
    combined_clinical$arid1a_event[i]<-"arid1a_event"
  }
  if ((combined_clinical$ARID1A[i] == "no_ARID1A_mutation" && combined_clinical$arid1a_del[i] == "arid1a_deletion")==TRUE){
    combined_clinical$arid1a_event[i]<-"arid1a_event"
  }
  
}


combined_clinical$SurvObj <- with(combined_clinical, Surv(as.numeric(as.vector(combined_clinical[,which(colnames(mut_clinical_data)=="_OS")])) ,as.numeric(as.vector(combined_clinical[,which(colnames(mut_clinical_data)=="_EVENT")]) == 1)))
km.by.arid1a <- survfit(SurvObj ~ arid1a_event, data = combined_clinical, conf.type = "log-log")
pval <- survdiff(SurvObj ~ arid1a_event, data = combined_clinical)
pval.by.arid1a <- 1-pchisq(pval$chisq,length(pval$n)-1)

#Make KM plot 
png(paste("plots/","arid1a_colon",".png",sep=""), width=297, height=210, units='mm', res=150)
plot(km.by.arid1a, mark.time=T, xlab='Days', ylab='Survival', col=c("dodgerblue4","gold2","red","green"),lwd=2, lty=1:2, yaxt="n", main="ARID1A in Colon")
axis(1, at = c(500,1500,2500,3500,4500,5500,6500))
axis(2, at = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), las=1)
legend("topright", c(paste("ARID1A mutation or deletion","(n =",length(which(combined_clinical$arid1a_event=="arid1a_event")),")"),paste("No ARID1A mutation or deletion","(n =",length(which(combined_clinical$arid1a_event=="no_arid1a_event")),")")), col=c("dodgerblue4","gold2"), lwd=2)
text(1900,0.9,paste("pvalue =",pval.by.arid1a))
dev.off()

#check for pten events
combined_clinical$pten_event <- "no_pten_event"
for (i in 1:length(combined_clinical$PTEN_mutation)) {
  if ((combined_clinical$PTEN_mutation[i] == "a_PTEN_mutation" && combined_clinical$PTEN_copynumber[i] == "a_PTEN_deletion")==TRUE){
    combined_clinical$pten_event[i]<-"a_pten_double_event"
  }
  if ((combined_clinical$PTEN_mutation[i] == "a_PTEN_mutation" && combined_clinical$PTEN_copynumber[i] == "no_PTEN_deletion")==TRUE){
    combined_clinical$pten_event[i]<-"a_pten_event"
  }
  if ((combined_clinical$PTEN_mutation[i] == "no_PTEN_mutation" && combined_clinical$PTEN_copynumber[i] == "a_PTEN_deletion")==TRUE){
    combined_clinical$pten_event[i]<-"a_pten_event"
  }
  
}


combined_clinical$SurvObj <- with(combined_clinical, Surv(as.numeric(as.vector(combined_clinical[,which(colnames(mut_clinical_data)=="_OS")])) ,as.numeric(as.vector(combined_clinical[,which(colnames(mut_clinical_data)=="_EVENT")]) == 1)))
km.by.pten <- survfit(SurvObj ~ pten_event, data = combined_clinical, conf.type = "log-log")
pval <- survdiff(SurvObj ~ pten_event, data = combined_clinical)
pval.by.pten <- 1-pchisq(pval$chisq,length(pval$n)-1)

#Make KM plot 
png(paste("plots/","pten_colon",".png",sep=""), width=297, height=210, units='mm', res=150)
plot(km.by.pten, mark.time=T, xlab='Days', ylab='Survival', col=c("dodgerblue4","gold2","red","green"),lwd=2, lty=1:2, yaxt="n", main="PTEN in Colon")
axis(1, at = c(500,1500,2500,3500,4500,5500,6500))
axis(2, at = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), las=1)
legend("topright", c(paste("PTEN mutation or deletion","(n =",length(which(combined_clinical$pten_event=="a_pten_event")),")"),paste("No PTEN mutation or deletion","(n =",length(which(combined_clinical$pten_event=="no_pten_event")),")")), col=c("dodgerblue4","gold2"), lwd=2)
text(1900,0.9,paste("pvalue =",pval.by.pten))
dev.off()

#check for events
combined_clinical$combined_event <- "no_combined_event"
for (i in 1:length(combined_clinical$pten_event)) {
  if ((combined_clinical$pten_event[i] == "a_pten_event" && combined_clinical$arid1a_event[i] == "a_arid1a_event")==TRUE){
    combined_clinical$combined_event[i]<-"a_double_combined_event"
  }
  if ((combined_clinical$pten_event[i] == "a_pten_event" && combined_clinical$arid1a_event[i] == "no_arid1a_event")==TRUE){
    combined_clinical$combined_event[i]<-"a_combined_event"
  }
  if ((combined_clinical$pten_event[i] == "no_pten_event" && combined_clinical$arid1a_event[i] == "a_arid1a_event")==TRUE){
    combined_clinical$combined_event[i]<-"a_combined_event"
  }
  
}


combined_clinical$SurvObj <- with(combined_clinical, Surv(as.numeric(as.vector(combined_clinical[,which(colnames(mut_clinical_data)=="_OS")])) ,as.numeric(as.vector(combined_clinical[,which(colnames(mut_clinical_data)=="_EVENT")]) == 1)))
km.by.pten <- survfit(SurvObj ~ combined_event, data = combined_clinical, conf.type = "log-log")
pval <- survdiff(SurvObj ~ combined_event, data = combined_clinical)
pval.by.pten <- 1-pchisq(pval$chisq,length(pval$n)-1)

#Make KM plot 
png(paste("plots/","combined_events_colon",".png",sep=""), width=297, height=210, units='mm', res=150)
plot(km.by.pten, mark.time=T, xlab='Days', ylab='Survival', col=c("dodgerblue4","gold2","red","green"),lwd=2, lty=1:2, yaxt="n", main="PTEN in Colon")
axis(1, at = c(500,1500,2500,3500,4500,5500,6500))
axis(2, at = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), las=1)
legend("topright", c(paste("ARID1A or PTEN (mutation or deletion)","(n =",length(which(combined_clinical$combined_event=="a_combined_event")),")"),paste("No PTEN mutation or deletion","(n =",length(which(combined_clinical$combined_event=="no_combined_event")),")")), col=c("dodgerblue4","gold2"), lwd=2)
text(1900,0.9,paste("pvalue =",pval.by.pten))
dev.off()