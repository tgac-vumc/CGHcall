preprocess <- 
function (input, maxmiss = 30, nchrom = 23, ...) 
{
    input <- input[!is.na(chromosomes(input)) & !is.na(bpstart(input))]
    input <- input[chromosomes(input) != 0 & bpstart(input) != 
        0, ]
    input <- input[chromosomes(input) <= nchrom, ]
    countMissing <- function(x, maxMissing) {
        if (length(x[is.na(x)]) <= maxMissing) 
            return(TRUE)
        return(FALSE)
    }
    nummis <- length(copynumber(input)[is.na(copynumber(input))])
    if (nummis > 0) {
        allowed <- ncol(input) * maxmiss/100
        filter <- apply(copynumber(input), 1, countMissing, allowed)
        input <- input[filter, ]
    }
    nummis <- length(copynumber(input)[is.na(copynumber(input))])
    if (nummis > 0) {
        if (!exists("k")) 
            k <- 10
        new.k <- ceiling((1 - maxmiss/100) * ncol(input))
        new.k <- min(new.k, k)
        if (new.k != 10) 
            cat("Changing impute.knn parameter k from", k, "to", 
                new.k, "due to small sample size.\n")
        if (new.k > 1) {
            nrows <- nrow(copynumber(input))
            nrpi <- floor(nrows/100)
            resultsall <- c()
            for (i in 1:(100 - 1)) {
                datforimp <- copynumber(input)[((i - 1) * nrpi + 
                  1):(nrpi * i), ]
                result <- impute::impute.knn(datforimp, k = new.k)
                resultsall <- rbind(resultsall, result$data)
            }
            datforimp <- copynumber(input)[((100 - 1) * nrpi + 
                1):nrows, ]
            result <- impute::impute.knn(datforimp, k = new.k)
            resultsall <- rbind(resultsall, result$data)
            copynumber(input) <- CGHcall:::.assignNames(resultsall, 
                input)
        }
    }
    input
}


normalize <- function(input, method="median", smoothOutliers=TRUE, ...) {

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
                density <- density(copynumber(input[,i]));
                value   <- density$x[which(density$y == max(density$y))];
                values  <- c(values, value);
            }        
        }
        matrixValues    <- matrix(rep(values, nrow(input)), ncol=ncol(input), byrow=TRUE);
        copynumber(input) <- copynumber(input) - matrixValues;
    }
    

    
    if (smoothOutliers) {
        cat("Smoothing outliers ... \n");
        CNA.object      <- DNAcopy::smooth.CNA(DNAcopy::CNA(copynumber(input), chromosomes(input), bpstart(input), data.type="logratio"), ...);
        for (i in 1:ncol(input)) {
            copynumber(input)[,i] <- CNA.object[[i+2]];
        }
    }
      
    #copynumber(input) <- CGHcall:::.assignNames(adjustForCellularity(copynumber(input), cellularity), input)
    
    input
    
}

#NEW: ALLOWS TWO SD.UNDO PARAMETERS, ONE FOR LONG AND ONE FOR SHORT SEGMENTS

segmentData <- function (input, clen=10, relSDlong=3, method = "DNAcopy", ...) 
{
    if (method == "DNAcopy") {
        CNA.object <- DNAcopy::CNA(copynumber(input), chromosomes(input), 
            bpstart(input), data.type = "logratio")
        cat("Start data segmentation .. \n")
        #segmented <- DNAcopy::segment(CNA.object, ...)
        segmented <- .segment2(CNA.object, clen, relSDlong, ...)
        numclone <- segmented$output$num.mark
        smrat <- segmented$output$seg
        numsmrat <- cbind(smrat, numclone)
        repdata <- function(row) {
            rep(row[1], row[2])
        }
        makelist <- apply(numsmrat, 1, repdata)
        joined <- unlist(makelist)
        rm(makelist)
        joined <- matrix(joined, ncol = ncol(input), byrow = FALSE)
        joined <- CGHcall:::.assignNames(joined, input)
        result <- CGHcall:::.segFromRaw(input, joined)
    }
    result
}

postsegnormalize <- function(segmentData,inter=c(-0.1,0.1)){
    #segmentData <- seg
    seg <- segmented(segmentData)
    
    values <- c()
    for (i in 1:ncol(seg)) {
                values <- c(values, median(seg[,i]));
            }    
    matrixValues    <- matrix(rep(values, nrow(seg)), ncol=ncol(seg), byrow=TRUE);
    seg <- seg - matrixValues #postseg works best when data are median normalized 
    countlevall <- apply(seg,2,function(x) {as.data.frame(table(x))})
    
    intcount <- function(int,sv){
    #int<-c(-0.5,0);sv<-segvec
        sv1 <- as.numeric(as.vector(sv[,1]))
        wh <- which(sv1<=int[2] & sv1>=int[1])
        return(sum(sv[wh,2]))
    }
    
    postsegnorm <- function(segvec,int=inter,intnr=3){
    #segvec<-countlevall[[1]];int=c(-0.30,0.1);intnr=3
        intlength <- (int[2]-int[1])/2
        gri <- intlength/intnr
        intst <- int[1]+(0:intnr)*gri
        intend <- intst+intlength
        ints <- cbind(intst,intend)
        intct <- apply(ints,1,intcount,sv=segvec)
        whmax <- which.max(intct)
        return(ints[whmax,]) 
    }
    
    postsegnorm_rec <- function(segvec,int,intnr=3){
    #segvec<-countlevall[[2]];int=c(-0.25,0.1);intnr=3
        newint <- postsegnorm(segvec,int,intnr)
        newint <- postsegnorm(segvec,newint,intnr)
        newint <- postsegnorm(segvec,newint,intnr)
        newint <- postsegnorm(segvec,newint,intnr)
        newint <- postsegnorm(segvec,newint,intnr)
        return(newint[1]+(newint[2]-newint[1])/2)
    }
    listres <- lapply(countlevall,postsegnorm_rec,int=inter)
    vecres <- c();for(i in 1:length(listres)){vecres <- c(vecres,listres[[i]])}
    
    segmented(segmentData) <- t(t(seg)-vecres)
    copynumber(segmentData) <- t(t(copynumber(segmentData)-matrixValues)-vecres)
    return(segmentData)
}
