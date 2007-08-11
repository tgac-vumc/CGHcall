preprocess <- function(input, type="file", maxmiss=30, nchrom=22, ...) {

    ## Version 1.0
    ## Author: Sjoerd Vosse
    ## Email: s.vosse@vumc.nl
    ## Date: 05-06-2007

    ## This function preprocesses your arrayCGH data:
    ## - It throws out data with unknown positions
    ## - It shrinks your data to nchrom chromosomes
    ## - It removes rows with more than maxmiss % missing values
    ## - It imputes missing values using knnimpute

    ## Some constants
    name        <- 1    ## Clone name column
    chromosome  <- 2    ## Chromosome column
    position    <- 3    ## Position column

    ## Read input data
    inputData   <- CGHcall:::.readInput(input, type)
    ## Delete data with unknown position or chromosome number
    inputData   <- inputData[!is.na(inputData[,chromosome]) & !is.na(inputData[,position]),];
    inputData   <- inputData[inputData[,chromosome] != 0 & inputData[,position] != 0,];
    
    ## Shrink to required number of chromosomes
    inputData   <- inputData[inputData[,chromosome] <= nchrom,];
    
    ## Make data matrix
    infoData    <- inputData[,c(name, chromosome, position)];
    matrixData  <- as.matrix(inputData[,-c(name, chromosome, position)]);    
    
    ## Remove rows with too many missing values
    countMissing <- function(x, maxMissing) {
        if (length(x[is.na(x)]) <= maxMissing) return(TRUE);
        return(FALSE);
    }
    
    allowed     <- ncol(matrixData) * maxmiss / 100;
    filter      <- apply(matrixData, 1, countMissing, allowed);
    infoData    <- infoData[filter,];
    matrixData  <- matrixData[filter,];
    
    ## Impute data using impute.knn from package 'impute'
    if (!exists("k")) k <- 10
    new.k       <- ceiling((1 - maxmiss / 100) * ncol(matrixData))
    new.k       <- min(new.k, k)
    if (new.k != 10) cat("Changing impute.knn parameter k from", k, "to", new.k, "due to small sample size.\n")
    result          <- impute::impute.knn(matrixData, k=new.k, ...);
    imputedData     <- cbind(infoData, result$data);
    
    return(imputedData);
    
}

normalize <- function(input, type="dataframe", method="median", cellularity=1, smoothOutliers=TRUE, ...) {

    ## Version 1.0
    ## Author: Sjoerd Vosse
    ## Email: s.vosse@vumc.nl
    ## Date: 05-06-2007

    ## This function normalizes the arrayCGH data using either global median, 
    ## global mean or no normalization. It can also adjust for the cellularity
    ## of your samples.

    ## Some constants
    name        <- 1;   ## Clone name column
    chromosome  <- 2;   ## Chromosome column
    position    <- 3;   ## Position column

    ## Read input data
    inputData   <- CGHcall:::.readInput(input, type)
    
    ## Make data matrix
    infoData    <- inputData[,c(name, chromosome, position)];
    matrixData  <- as.matrix(inputData[,-c(name, chromosome, position)]);
    
    ### Normalization
    values  <- c();
    if (method == "none") {
        cat("Skipping normalization ... \n");
        normalizedData <- matrixData;
    } else {
        if (method == "median") {
            cat("Applying median normalization ... \n");
            for (i in 1:ncol(matrixData)) {
                values <- c(values, median(matrixData[,i]));
            }    
        } else if (method == "mode") {
            cat("Applying mode normalization ... \n");        
            for (i in 1:ncol(matrixData)) {
                density <- density(matrixData[,i]);
                value   <- density$x[which(density$y == max(density$y))];
                values  <- c(values, value);
            }        
        }
        matrixValues    <- matrix(rep(values, nrow(matrixData)), ncol=ncol(matrixData), byrow=TRUE);
        normalizedData  <- matrixData - matrixValues;
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
        CNA.object      <- DNAcopy::smooth.CNA(DNAcopy::CNA(normalizedData, infoData[,chromosome], infoData[,position], data.type="logratio"), ...);
        smoothData      <- normalizedData;
        for (i in 1:ncol(matrixData)) {
            smoothData[,i] <- CNA.object[[i+2]];
        }
    } else {
        smoothData      <- normalizedData;
    }
    
    if (length(cellularity) < ncol(matrixData)) cellularity <- rep(cellularity, ncol(matrixData));
    if (length(cellularity) > ncol(matrixData)) cellularity <- cellularity[1:ncol(matrixData)];
    
    adjustedData <- adjustForCellularity(smoothData, cellularity);
    returnData   <- cbind(infoData, adjustedData);
    colnames(returnData) <- colnames(inputData);
    return(returnData);
    
}

segmentData <- function(input, type="dataframe", method="DNAcopy", ...) {

    ## Version 1.0
    ## Author: Sjoerd Vosse
    ## Email: s.vosse@vumc.nl
    ## Date: 05-06-2007

    ## This function is a simple wrapper that applies existing segmentation
    ## algorithms (currently only DNAcopy) to your arrayCGH data.

    ## Some constants
    name        <- 1;   ## Clone name column
    chromosome  <- 2;   ## Chromosome column
    position    <- 3;   ## Position column

    ## Read input data
    inputData   <- CGHcall:::.readInput(input, type)
    
    ## Make data matrix
    infoData    <- inputData[,c(name, chromosome, position)];
    matrixData  <- as.matrix(inputData[,-c(name, chromosome, position)]);
    
    if (method == "DNAcopy") {
    
        CNA.object  <- DNAcopy::CNA(matrixData, infoData[,chromosome], infoData[,position], data.type="logratio");
    
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
        joined      <- matrix(joined, ncol=ncol(matrixData), byrow=FALSE)
        segmented   <- data.frame(infoData, joined)
        colnames(segmented) <- colnames(inputData)
    }
    return(segmented);
}
