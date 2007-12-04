setMethod("initialize", "cghCall",
        function(.Object,
                assayData       = assayDataNew(copynumber=copynumber, segmented=segmented, calls=calls, probloss=probloss, probnorm=probnorm, probgain=probgain, ...),
                phenoData       = annotatedDataFrameFrom(assayData, byrow=FALSE),
                featureData     = CGHcall:::.makeEmptyFeatureData(assayData),
                experimentData  = new("MIAME"),
                annotation      = character(),
                copynumber      = new("matrix"),
                segmented       = new("matrix"),
                calls           = new("matrix"),
                probloss        = new("matrix"),
                probnorm        = new("matrix"),
                probgain        = new("matrix"),
                ... ) {
            callNextMethod(.Object,
                            assayData       = assayData,
                            phenoData       = phenoData,
                            featureData     = featureData,
                            experimentData  = experimentData,
                            annotation      = annotation)
})

setValidity("cghCall", function(object) {
    msg <- Biobase:::validMsg(NULL, Biobase:::isValidVersion(object, "cghCall"))
    msg <- Biobase:::validMsg(msg, Biobase:::assayDataValidMembers(assayData(object), c("copynumber", 
                                                                                        "segmented", 
                                                                                        "calls", 
                                                                                        "probloss", 
                                                                                        "probnorm", 
                                                                                        "probgain"
                                                                                        )))
    msg <- Biobase:::validMsg(msg, CGHcall:::.featureDataRequiredColumns(featureData(object), c("Chromosome", "Start", "End")))
    if (is.null(msg)) TRUE else msg
})

setMethod("plot", signature(x="cghCall", y="missing"),
function (x, y, ... )
{

    calls           <- calls(x)
    nsamples        <- ncol(x)
    if (!is.null(probamp(x))) nclass <- 4
    else nclass     <- 3
    chrom           <- chromosomes(x)
    chrom.labels    <- unique(chrom)
    nclone          <- length(chrom)

    for (i in 1:ncol(x)) {
    
        cat("Plotting sample", sampleNames(x)[i], "\n")
    
        genomdat        <- copynumber(x)[,i]
        probsdraw       <- cbind(probloss(x)[,i], probnorm(x)[,i], probgain(x)[,i])
        if (!is.null(probamp(x))) probsdraw <- cbind(probsdraw, probamp(x)[,i])
        lt              <- 0
        
        if (nclass==4) {
            ticksamp    <- which(probsdraw[,nclass] >= 0.5)
            lt          <- length(ticksamp)        
            probsdraw   <- cbind(probsdraw[,1:2], probsdraw[,3] + probsdraw[,4])
        }
           
        segment         <- CGHcall:::.makeSegments(cbind(chromosomes(x), segmented(x)[,i]))
        
        widths          <- segment[,3] - segment[,2] + 1
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
            axis(3,at=ticksamp, labels=FALSE, col = "blue", col.axis="black", srt=270, las=1, cex.axis=1, cex.lab=1)
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
        
        title(sampleNames(x)[i])
        
        ### Add log2ratios
        points((1:nclone)-.5,genomdat,cex=.1)
        
        ### X-axis with chromosome labels
        ax<-(cumsum(table(chrom))+c(0, cumsum(table(chrom))[-length(cumsum(table(chrom)))]))/2
        axis(side=1, at=ax, labels=chrom.labels, cex=.2, lwd=.5, las=1, cex.axis=1, cex.lab=1) # bottom axis
            
        ### Blue lines for segment means
        for (jjj in (1:nrow(segment)))
        segments(segment[jjj,2], segment[jjj,1], segment[jjj,3], segment[jjj,1], col="blue", lwd=3)        
    }
})

setMethod("plot.summary", signature(x="cghCall", y="missing"),
function (x, y, ... )
{
    
    calls           <- calls(x)
    n.sam           <- ncol(x)
    if (!is.null(probamp(x))) nclass <- 4
    else nclass     <- 3
    chrom           <- chromosomes(x)
    
    nclone          <- length(chrom)
    chrom.labels    <- as.character(unique(chrom))
    
    prob.loss.ind   <- 0
    prob.none.ind   <- 0
    prob.gain.ind   <- 0
    prob.amp.ind    <- 0
    
    for (j in 1:ncol(x)) {
        cat("Adding sample", sampleNames(x)[j], "to summary plot.\n");
        prob.loss.ind   <- prob.loss.ind + probloss(x)[,j]
        prob.none.ind   <- prob.none.ind + probnorm(x)[,j]
        prob.gain.ind   <- prob.gain.ind + probgain(x)[,j]
        if (nclass==4) {
            prob.amp.ind    <- prob.amp.ind + probamp(x)[,j]
            prob.gain.amp   <- prob.gain.ind + prob.amp.ind
        }
    }
    
    if (nclass==4) {
        p2 <- prob.amp.ind/n.sam
        ticksamp <- which(p2 >= 0.5)
        lt <- length(ticksamp)
        probsdraw <- cbind(prob.loss.ind, prob.none.ind, prob.gain.amp) / n.sam
    }
    if (nclass==3) {
        lt <- 0 
        probsdraw <- cbind(prob.loss.ind, prob.none.ind, prob.gain.ind) / n.sam
    }
    
    par(mar = c(5, 4, 4, 4) + 0.2)
    barplot(t(probsdraw), border=F, space=0, col=c("red","white","green"), las=1, cex.axis=1, cex.lab=1, xaxt="n")
    
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
    
})

setMethod("chromosomes", "cghCall", function(object) pData(featureData(object))[,"Chromosome"])
setMethod("bpstart", "cghCall", function(object) pData(featureData(object))[,"Start"])
setMethod("bpend", "cghCall", function(object) pData(featureData(object))[,"End"])

setMethod("copynumber", signature(object="cghCall"),
        function(object) Biobase:::assayDataElement(object, "copynumber"))
        
setReplaceMethod("copynumber", signature(object="cghCall", value="matrix"),
                function(object, value) assayDataElementReplace(object, "copynumber", value))   

setMethod("segmented", signature(object="cghCall"),
        function(object) Biobase:::assayDataElement(object, "segmented"))
        
setReplaceMethod("segmented", signature(object="cghCall", value="matrix"),
                function(object, value) assayDataElementReplace(object, "segmented", value))  
                
setMethod("calls", signature(object="cghCall"),
        function(object) Biobase:::assayDataElement(object, "calls"))
        
setReplaceMethod("calls", signature(object="cghCall", value="matrix"),
                function(object, value) assayDataElementReplace(object, "calls", value))                  

setMethod("probloss", signature(object="cghCall"),
        function(object) Biobase:::assayDataElement(object, "probloss"))
        
setReplaceMethod("probloss", signature(object="cghCall", value="matrix"),
                function(object, value) assayDataElementReplace(object, "probloss", value))
                
setMethod("probnorm", signature(object="cghCall"),
        function(object) Biobase:::assayDataElement(object, "probnorm"))
        
setReplaceMethod("probnorm", signature(object="cghCall", value="matrix"),
                function(object, value) assayDataElementReplace(object, "probnorm", value))
                
setMethod("probgain", signature(object="cghCall"),
        function(object) Biobase:::assayDataElement(object, "probgain"))
        
setReplaceMethod("probgain", signature(object="cghCall", value="matrix"),
                function(object, value) assayDataElementReplace(object, "probgain", value))    
                
setMethod("probamp", signature(object="cghCall"),
        function(object) Biobase:::assayDataElement(object, "probamp"))
        
setReplaceMethod("probamp", signature(object="cghCall", value="matrix"),
                function(object, value) assayDataElementReplace(object, "probamp", value))  
