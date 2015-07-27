library(RColorBrewer)

###Functions

cor.mat.by.region <- function(x, m, build='GRCh37', main='Correlation Matrix', ...) {
  image(x=1:(nrow(m)+1), y=1:(nrow(m)+1), z=m, zlim=c(-1,1), col=brewer.pal(10, 'RdBu'), xlab=NA, ylab=NA, useRaster=TRUE, main=main, xaxt='n', yaxt='n', ...)
  
  chr.change <- chromosomes(x) != c(chromosomes(x)[-1], 0)
  chr.change <- which(chr.change) + 1
  abline(h=chr.change[-length(chr.change)], v=chr.change[-length(chr.change)], lty='dashed')
  
  if (nrow(m) < 100) {
    a <- add.cytobands(fData(x), genome.build=build)
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

cor.mat.by.size <- function(x, m, build='GRCh37', main='Correlation Matrix', ...) {
  chrom <- chromosomes(x)
  pos <- bpstart(x)
  pos2 <- bpend(x)
  uni.chrom <- unique(chrom)
  chrom.lengths <- CGHbase:::.getChromosomeLengths(build)[as.character(uni.chrom)]
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


