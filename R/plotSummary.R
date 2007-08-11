plotSummary <- function(CGHcall.result, samples="all", export="no") { 

    calls       <- CGHcall.result[[2]]
    dataprob    <- CGHcall.result[[1]]
    nsamples    <- ncol(calls) - 3
    nclass      <- ((ncol(dataprob) - 3) / nsamples) - 1
    chrom       <- dataprob[,2]
    
    if (samples == "all") {
        sample <- 1:nsamples
    } else {
        sample  <- as.numeric(samples)
    }
    
    ncpp    <- nclass+1
    n.sam   <- length(sample)
    nrowdp  <- nrow(dataprob)

    nclone          <- length(chrom)
    chrom.labels    <- as.character(unique(chrom))
    
    prob.loss.ind   <- 0
    prob.none.ind   <- 0
    prob.gain.ind   <- 0
    prob.amp.ind    <- 0
    
    for (j in sample) {
        cat("Adding sample", j, "to summary plot.\n");
        jj <- sample[j]
        prob.loss.ind   <- prob.loss.ind + dataprob[,(3+ncpp*(jj-1)+2)]
        prob.none.ind   <- prob.none.ind + dataprob[,(3+ncpp*(jj-1)+3)]
        prob.gain.ind   <- prob.gain.ind + dataprob[,(3+ncpp*(jj-1)+4)]
        if (nclass==4) {
            prob.amp.ind    <- prob.amp.ind + dataprob[,(3+ncpp*(jj-1)+5)]
            prob.gain.amp   <- prob.gain.ind + prob.amp.ind
        }
    }
    
    if (nclass==4) {
        p2 <- prob.amp.ind/n.sam
        ticksamp <- which(p2 >= 0.5)
        lt <- length(ticksamp)
        probsdraw <- cbind(prob.loss.ind,prob.none.ind,prob.gain.amp)/n.sam
    }
    if (nclass==3) {
        lt <- 0 
        probsdraw <- cbind(prob.loss.ind,prob.none.ind,prob.gain.ind)/n.sam
    }
    
    if (export != "no") postscript(export)
    
    par(mar = c(5, 4, 4, 4) + 0.2)
    barplot(t(probsdraw), border=F, space=0, col=c("red","white","green"), las=1, cex.axis=1, cex.lab=1)
    
    lim <- par("usr")
    
    lim[3:4] <- c(-5,5)
    
    par(usr=lim)
    if (lt > 0) {
        axis(3,at=ticksamp, labels=FALSE, col = "blue", col.axis="black", srt=270, las=1, cex.axis=1, cex.lab=1)
    }
    box()
    mtext("Mean Probability",side=2,line=3,srt=270)
    
    ax<-(cumsum(table(chrom))+c(0,cumsum(table(chrom))[-length(cumsum(table(chrom)))]))/2
    axis(side=1, at=ax, labels=chrom.labels, cex=.2, lwd=.5, las=1, cex.axis=1, cex.lab=1)
    
    title("Summary plot")
    
    if (export != "no") dev.off()
}
