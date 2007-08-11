plotProfile <- function(CGHcall.result, samples="all", export="no") {

    dataprob        <- CGHcall.result[[1]]
    calls           <- CGHcall.result[[2]]
    segment.data    <- CGHcall.result[[3]]
    nsamples        <- ncol(calls) - 3
    nclass          <- ((ncol(dataprob) - 3) / nsamples) - 1
    chrom           <- dataprob[,2]
    chrom.labels    <- unique(chrom)
    nclone          <- length(chrom)
    
    if (samples[1] == "all") samples <- 1:nsamples
    if (export != "no") postscript(export)

    for (i in samples) {
    
        cat("Plotting sample", i, "\n")
    
        start.index     <- 4+(i-1)*(nclass+1)
        genomdat        <- dataprob[,start.index]
        probsdraw       <- dataprob[,(start.index+1):(start.index+nclass)]
        lt              <- 0
        
        if (nclass==4) {
            ticksamp    <- which(probsdraw[,nclass] >= 0.5)
            lt          <- length(ticksamp)        
            probsdraw   <- cbind(probsdraw[,1:2], probsdraw[,3] + probsdraw[,4])
        }
           
        segment         <- segment.data[segment.data[,1] == i,]
        
        widths          <- segment[,4] - segment[,3] + 1
        plot.data       <- unique(probsdraw)
    
        par(mar=c(5, 4, 4, 4) + 0.2)
        
        ### Plot the probability bars
        barplot(t(plot.data), width=widths, border=F, space=0, col=c("red","white","green"), las=1, cex.axis=1, cex.lab=1, xaxt="n")
        
        lim <- par("usr")
        lim[3:4] <- c(-5, 5)
        par(usr=lim)
        dticks <- seq(-5, 5, by=1)
        axis(4, at=dticks, labels=dticks, srt=270, las=1, cex.axis=1, cex.lab=1)
        if (lt > 0) {
            axis(3,at=ticksamp, labels=FALSE, col = "blue", col.axis="black",srt=270,las=1,cex.axis=1,cex.lab=1)
        }        
        box()
        
        ### Add axis labels
        mtext("log2 ratio", side=4, line=3, srt=270)
        mtext("probability", side=2, line=3, srt=270)
        
        #### add vert lines at chromosome ends
        abline(h=0) 
        for (iii in 1:length(cumsum(table(chrom)))) {
            segments(cumsum(table(chrom))[[iii]], -5, cumsum(table(chrom))[[iii]], 5, lty=2)
        }
        
        title(names(calls)[i+3])
        
        ### Add log2ratios
        points((1:nclone)-.5,genomdat,cex=.1)
        
        ### X-axis with chromosome labels
        ax<-(cumsum(table(chrom))+c(0, cumsum(table(chrom))[-length(cumsum(table(chrom)))]))/2
        axis(side=1, at=ax, labels=chrom.labels, cex=.2, lwd=.5, las=1, cex.axis=1, cex.lab=1) # bottom axis
            
        ### Blue lines for segment means
        for (jjj in (1:nrow(segment)))
        segments(segment[jjj,3], segment[jjj,2], segment[jjj,4], segment[jjj,2], col="blue", lwd=3)        
    }
    if (export != "no") dev.off()
}
