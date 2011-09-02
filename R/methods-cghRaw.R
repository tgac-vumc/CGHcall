setMethod("initialize", "cghRaw",
        function(.Object,
                assayData       = assayDataNew(copynumber=copynumber, ...),
                phenoData       = annotatedDataFrameFrom(assayData, byrow=FALSE),
                featureData     = CGHcall:::.makeEmptyFeatureData(assayData),
                experimentData  = new("MIAME"),
                annotation      = character(),
                copynumber      = new("matrix"),
                ... ) {
            callNextMethod(.Object,
                            assayData       = assayData,
                            phenoData       = phenoData,
                            featureData     = featureData,
                            experimentData  = experimentData,
                            annotation      = annotation)
})

setMethod("plot", signature(x="cghRaw", y="missing"),
function (x, y, ... )
{
    for (i in 1:ncol(x)) {
        cat("Plotting sample", sampleNames(x)[i], "\n")
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
    }
})

setValidity("cghRaw", function(object) {
    msg <- Biobase:::validMsg(NULL, Biobase:::isValidVersion(object, "cghRaw"))
    msg <- Biobase:::validMsg(msg, Biobase:::assayDataValidMembers(assayData(object), c("copynumber")))
    msg <- Biobase:::validMsg(msg, CGHcall:::.featureDataRequiredColumns(featureData(object), c("Chromosome", "Start", "End")))
    if (is.null(msg)) TRUE else msg
})

setMethod("chromosomes", signature(object="cghRaw"), function(object) pData(featureData(object))[,"Chromosome"])
setMethod("bpstart", signature(object="cghRaw"), function(object) pData(featureData(object))[,"Start"])
setMethod("bpend", signature(object="cghRaw"), function(object) pData(featureData(object))[,"End"])

setMethod("copynumber", signature(object="cghRaw"),
        function(object) Biobase:::assayDataElement(object, "copynumber"))
        
setReplaceMethod("copynumber", signature(object="cghRaw", value="matrix"),
                function(object, value) assayDataElementReplace(object, "copynumber", value))    
