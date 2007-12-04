setMethod("initialize", "cghSeg",
        function(.Object,
                assayData       = assayDataNew(copynumber=copynumber, segmented=segmented, ...),
                phenoData       = annotatedDataFrameFrom(assayData, byrow=FALSE),
                featureData     = CGHcall:::.makeEmptyFeatureData(assayData),
                experimentData  = new("MIAME"),
                annotation      = character(),
                copynumber      = new("matrix"),
                segmented       = new("matrix"),
                ... ) {
            callNextMethod(.Object,
                            assayData       = assayData,
                            phenoData       = phenoData,
                            featureData     = featureData,
                            experimentData  = experimentData,
                            annotation      = annotation)
})

setMethod("plot", signature(x="cghSeg", y="missing"),
function (x, y, ... )
{
    for (i in 1:ncol(x)) {
        cat("Plotting sample", sampleNames(x)[i], "\n")
        segment         <- CGHcall:::.makeSegments(cbind(chromosomes(x), segmented(x)[,i]))
        chrom           <- chromosomes(x)
        data            <- data.frame(chrom, bpstart(x), copynumber(x)[,i])
        colnames(data)  <- c("chromosome", "position", "ratio")
        chrom.labels    <- as.character(unique(chrom))
        plot(data[,3], pch=".", main=sampleNames(x)[i], ylab="log2ratio", xlab="chromosomes", ylim=c(-2,5), xaxt="n", xaxs="i")
        abline(h=0) 
        for (iii in 1:length(cumsum(table(chrom)))) {
            segments(cumsum(table(chrom))[[iii]],-5,cumsum(table(chrom))[[iii]],5,lty=2)
        }
        ax<-(cumsum(table(chrom))+c(0,cumsum(table(chrom))[-length(cumsum(table(chrom)))]))/2
        axis(side=1,at=ax,labels=chrom.labels,cex=.2,lwd=.5,las=1,cex.axis=1,cex.lab=1)
        for (jjj in (1:nrow(segment))) {
            segments(segment[jjj,2], segment[jjj,1], segment[jjj,3], segment[jjj,1], col="blue", lwd=3)        
        }
    }
})

setValidity("cghSeg", function(object) {
    msg <- Biobase:::validMsg(NULL, Biobase:::isValidVersion(object, "cghSeg"))
    msg <- Biobase:::validMsg(msg, Biobase:::assayDataValidMembers(assayData(object), c("copynumber", "segmented")))
    msg <- Biobase:::validMsg(msg, CGHcall:::.featureDataRequiredColumns(featureData(object), c("Chromosome", "Start", "End")))
    if (is.null(msg)) TRUE else msg
})

setMethod("chromosomes", "cghSeg", function(object) pData(featureData(object))[,"Chromosome"])
setMethod("bpstart", "cghSeg", function(object) pData(featureData(object))[,"Start"])
setMethod("bpend", "cghSeg", function(object) pData(featureData(object))[,"End"])

setMethod("copynumber", signature(object="cghSeg"),
        function(object) Biobase:::assayDataElement(object, "copynumber"))
        
setReplaceMethod("copynumber", signature(object="cghSeg", value="matrix"),
                function(object, value) assayDataElementReplace(object, "copynumber", value))    

setMethod("segmented", signature(object="cghSeg"),
        function(object) Biobase:::assayDataElement(object, "segmented"))
        
setReplaceMethod("segmented", signature(object="cghSeg", value="matrix"),
                function(object, value) assayDataElementReplace(object, "segmented", value))  
