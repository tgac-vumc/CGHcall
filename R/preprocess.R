make_cghRaw <- function(input) {
    if (class(input) == "character") input  <- read.table(input, header=T, sep="\t", fill=T, quote="")
    copynumber  <- as.matrix(input[,5:ncol(input)])
    rownames(copynumber) <- input[,1]
    annotation  <- data.frame(Chromosome=input[,2], Start=input[,3], End=input[,4], row.names=input[,1])
    metadata    <- data.frame(labelDescription=c("Chromosomal position", "Basepair position start", "Basepair position end"), row.names=c("Chromosome", "Start", "End"))    
    dimLabels   <- c("featureNames", "featureColumns")
    annotation  <- new("AnnotatedDataFrame", data=annotation, dimLabels=dimLabels, varMetadata=metadata)   
    result      <- new("cghRaw", copynumber=copynumber, featureData=annotation)
}

preprocess <- function(input, maxmiss=30, nchrom=22, ...) {

    ## Version 2.0
    ## Author: Sjoerd Vosse
    ## Email: s.vosse@vumc.nl
    ## Date: 23-10-2007
    
    ## Changes since previous version
    ## - input is of class cghRaw

    ## This function preprocesses your arrayCGH data:
    ## - It throws out data with unknown positions
    ## - It shrinks your data to nchrom chromosomes
    ## - It removes rows with more than maxmiss % missing values
    ## - It imputes missing values using knnimpute
    
    ## Sort data by bppos, just to make sure
    input   <- input[order(bpstart(input)),]
    input   <- input[order(chromosomes(input)),]

    ## Delete data with unknown position or chromosome number
    input   <- input[!is.na(chromosomes(input)) & !is.na(bpstart(input))];
    input   <- input[chromosomes(input) != 0 & bpstart(input) != 0,];
    
    ## Shrink to required number of chromosomes
    input   <- input[chromosomes(input) <= nchrom,];
    
    ## Remove rows with too many missing values
    countMissing <- function(x, maxMissing) {
        if (length(x[is.na(x)]) <= maxMissing) return(TRUE);
        return(FALSE);
    }
    
    allowed     <- ncol(input) * maxmiss / 100;
    filter      <- apply(copynumber(input), 1, countMissing, allowed);
    input       <- input[filter,]
    
    ## Impute data using impute.knn from package 'impute'
    if (!exists("k")) k <- 10
    new.k       <- ceiling((1 - maxmiss / 100) * ncol(input))
    new.k       <- min(new.k, k)
    if (new.k != 10) cat("Changing impute.knn parameter k from", k, "to", new.k, "due to small sample size.\n")
    result      <- impute::impute.knn(copynumber(input), k=new.k, ...);
    copynumber(input) <- CGHcall:::.assignNames(result$data, input)
    
    input
    
}

normalize <- function(input, method="median", cellularity=1, smoothOutliers=TRUE, ...) {

    ## Version 2.0
    ## Author: Sjoerd Vosse
    ## Email: s.vosse@vumc.nl
    ## Date: 24-10-2007
    
    ## Changes since previous version
    ## - input is now cghRaw class

    ## This function normalizes the arrayCGH data using either global median, 
    ## global mean or no normalization. It can also adjust for the cellularity
    ## of your samples.
    
    ### Normalization
    values  <- c();
    if (method == "none") {
        cat("Skipping normalization ... \n");
    } else {
        if (method == "median") {
            cat("Applying median normalization ... \n");
            for (i in 1:ncol(input)) {
                values <- c(values, median(copynumber(input)[,i]));
            }    
        } else if (method == "mode") {
            cat("Applying mode normalization ... \n");        
            for (i in 1:ncol(input)) {
                density <- density(input[,i]);
                value   <- density$x[which(density$y == max(density$y))];
                values  <- c(values, value);
            }        
        }
        matrixValues    <- matrix(rep(values, nrow(input)), ncol=ncol(input), byrow=TRUE);
        copynumber(input) <- copynumber(input) - matrixValues;
    }
    
    ### Adjusting for cellularity
    adjustForCellularity <- function(matrix, cellularity) {
        cat("Adjusting for cellularity ... \n");
        result  <- c();
        adjustCellularity <- function(value, cellularity) {
            corrected   <- (2^value / cellularity - (1 - cellularity) / cellularity)
            if (corrected < 2^(-5)) {
                corrected <- 2^value;
            }
            new.value   <- log2(corrected)
            return(new.value)
        }
        for (i in 1:ncol(matrix)) {
            cat("Cellularity sample", i, ": ", cellularity[i], "\n");
            if (cellularity[i] < 1) {
                new.column  <- sapply(matrix[,i], adjustCellularity, cellularity[i]);
                result      <- cbind(result, new.column);
            } else {
                result      <- cbind(result, matrix[,i]);
            }
        }
        return(result);
    }
    
    if (smoothOutliers) {
        cat("Smoothing outliers ... \n");
        CNA.object      <- DNAcopy::smooth.CNA(DNAcopy::CNA(copynumber(input), chromosomes(input), bpstart(input), data.type="logratio"), ...);
        for (i in 1:ncol(input)) {
            copynumber(input)[,i] <- CNA.object[[i+2]];
        }
    }
    
    if (length(cellularity) < ncol(input)) cellularity <- rep(cellularity, ncol(input));
    if (length(cellularity) > ncol(input)) cellularity <- cellularity[1:ncol(input)];
    
    copynumber(input) <- CGHcall:::.assignNames(adjustForCellularity(copynumber(input), cellularity), input)
    
    input
    
}

segmentData <- function(input, method="DNAcopy", ...) {

    ## Version 2.0
    ## Author: Sjoerd Vosse
    ## Email: s.vosse@vumc.nl
    ## Date: 24-10-2007

    ## This function is a simple wrapper that applies existing segmentation
    ## algorithms (currently only DNAcopy) to your arrayCGH data.

    ## Changes since previous version
    ## - input is now cghRaw class
    ## - output is cghSeg class
    
    if (method == "DNAcopy") {
    
        CNA.object  <- DNAcopy::CNA(copynumber(input), chromosomes(input), bpstart(input), data.type="logratio");
    
        cat("Start data segmentation .. \n");
        segmented  <- DNAcopy::segment(CNA.object, ...);
        
        ## Convert DNAcopy result to CGHcall format
        numclone    <- segmented$output$num.mark
        smrat       <- segmented$output$seg
        numsmrat    <- cbind(smrat, numclone)
        repdata     <- function(row) {
            rep(row[1], row[2])
        }
        makelist    <- apply(numsmrat, 1, repdata)
        joined      <- c()
        for (j in 1:length(makelist)) {
            joined  <- c(joined, makelist[[j]])
        }
        joined      <- matrix(joined, ncol=ncol(input), byrow=FALSE)
        joined      <- CGHcall:::.assignNames(joined, input)
        result      <- CGHcall:::.segFromRaw(input, joined)
    }
    
    result
}
