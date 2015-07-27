setMethod("frequencyPlot", signature(x="cghRegions", y="missing"),
          function (x, y, main='Frequency Plot', gaincol='blue', losscol='red', misscol=NA, build='GRCh36',... )
          {
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