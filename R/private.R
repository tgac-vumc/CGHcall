.segFromRaw <- function(raw, segmatrix) {
    result <- new("cghSeg", copynumber=copynumber(raw), 
                            segmented=segmatrix, 
                            phenoData=phenoData(raw), 
                            featureData=featureData(raw), 
                            annotation=annotation(raw), 
                            experimentData=experimentData(raw)
                            )
    result
}

.callFromSeg <- function(seg, assayData) {
    result <- new("cghCall", assayData=assayData,
                            phenoData=phenoData(seg), 
                            featureData=featureData(seg), 
                            annotation=annotation(seg), 
                            experimentData=experimentData(seg)
                            )
    result
}

######## COPIED FROM DNACOPY, BUT ADAPTED TO ALLOW FOR TWO UNDO-PARAMETERS

.segment2 <- function (x,clen=10, relSDlong=3, alpha = 0.01, nperm = 10000, p.method = c("hybrid", 
    "perm"), min.width = 2, kmax = 25, nmin = 200, eta = 0.05, 
    sbdry = NULL, trim = 0.025, undo.splits = c("none", "prune", 
        "sdundo"), undo.prune = 0.05, undo.SD = 3, verbose = 1) 
{
    if (!inherits(x, "CNA")) 
        stop("First arg must be a copy number array object")
    call <- match.call()
    if (min.width < 2 | min.width > 5) 
     stop("minimum segment width should be between 2 and 5")
    if (missing(sbdry)) {
        if (nperm == 10000 & alpha == 0.01 & eta == 0.05) {
            sbdry <- default.DNAcopy.bdry
        }
        else {
            max.ones <- floor(nperm * alpha) + 1
            sbdry <- DNAcopy:::getbdry(eta, nperm, max.ones)
        }
    }
    sbn <- length(sbdry)
    nsample <- ncol(x) - 2
    sampleid <- colnames(x)[-(1:2)]
    uchrom <- unique(x$chrom)
    data.type <- attr(x, "data.type")
    p.method <- match.arg(p.method)
    undo.splits <- match.arg(undo.splits)
    segres <- list()
    segres$data <- x
    allsegs <- list()
    allsegs$ID <- NULL
    allsegs$chrom <- NULL
    allsegs$loc.start <- NULL
    allsegs$loc.end <- NULL
    allsegs$num.mark <- NULL
    allsegs$seg.mean <- NULL
    for (isamp in 1:nsample) {
        if (verbose >= 1) 
            cat(paste("Analyzing:", sampleid[isamp], "\n"))
        genomdati <- x[, isamp + 2]
        ina <- which(!is.na(genomdati) & !(abs(genomdati) == 
            Inf))
        genomdati <- genomdati[ina]
        trimmed.SD <- sqrt(DNAcopy:::trimmed.variance(genomdati, trim))
        chromi <- x$chrom[ina]
        sample.lsegs <- NULL
        sample.segmeans <- NULL
        for (ic in uchrom) {
            if (verbose >= 2) 
                cat(paste("  current chromosome:", ic, "\n"))
            segci <- .changepoints2(genomdati[chromi == ic], data.type, 
                alpha, sbdry, sbn, nperm, p.method, min.width, 
                kmax, nmin, trimmed.SD, undo.splits, undo.prune, 
                undo.SD, clen, relSDlong, verbose)
            sample.lsegs <- c(sample.lsegs, segci$lseg)
            sample.segmeans <- c(sample.segmeans, segci$segmeans)
        }
        sample.nseg <- length(sample.lsegs)
        sample.segs.start <- ina[cumsum(c(1, sample.lsegs[-sample.nseg]))]
        sample.segs.end <- ina[cumsum(sample.lsegs)]
        allsegs$ID <- c(allsegs$ID, rep(isamp, sample.nseg))
        allsegs$chrom <- c(allsegs$chrom, x$chrom[sample.segs.end])
        allsegs$loc.start <- c(allsegs$loc.start, x$maploc[sample.segs.start])
        allsegs$loc.end <- c(allsegs$loc.end, x$maploc[sample.segs.end])
        allsegs$num.mark <- c(allsegs$num.mark, sample.lsegs)
        allsegs$seg.mean <- c(allsegs$seg.mean, sample.segmeans)
    }
    allsegs$ID <- sampleid[allsegs$ID]
    allsegs$seg.mean <- round(allsegs$seg.mean, 4)
    allsegs <- as.data.frame(allsegs)
    allsegs$ID <- as.character(allsegs$ID)
    segres$output <- allsegs
    segres$call <- call
    class(segres) <- "DNAcopy"
    segres
}


.changepoints2 <- function (genomdat, data.type = "logratio", alpha = 0.01, sbdry, 
    sbn, nperm = 10000, p.method = "hybrid", min.width = 2, kmax = 25, 
    nmin = 200, trimmed.SD = NULL, undo.splits = "none", undo.prune = 0.05, 
    undo.SD = 3, clen, relSDlong, verbose = 1, ngrid = 100, tol = 1e-06) 
{
    n <- length(genomdat)
    if (missing(trimmed.SD)) 
        trimmed.SD <- mad(diff(genomdat))/sqrt(2)
    seg.end <- c(0, n)
    k <- length(seg.end)
    change.loc <- NULL
    while (k > 1) {
        current.n <- seg.end[k] - seg.end[k - 1]
        if (verbose >= 3) 
            cat(".... current segment:", seg.end[k - 1] + 1, 
                "-", seg.end[k], "\n")
        if (current.n >= 2 * min.width) {
            current.genomdat <- genomdat[(seg.end[k - 1] + 1):seg.end[k]]
            current.genomdat <- current.genomdat - mean(current.genomdat)
            current.tss <- sum(current.genomdat^2)
            hybrid <- FALSE
            delta <- 0
            if ((p.method == "hybrid") & (nmin < current.n)) {
                hybrid <- TRUE
                delta <- (kmax + 1)/current.n
            }
            zzz <- .Fortran("fndcpt", n = as.integer(current.n), 
                x = as.double(current.genomdat), tss = as.double(current.tss), 
                px = double(current.n), sx = double(current.n), 
                nperm = as.integer(nperm), cpval = as.double(alpha), 
                ncpt = integer(1), icpt = integer(2), ibin = as.logical(data.type == 
                  "binary"), hybrid = as.logical(hybrid), al0 = as.integer(min.width), 
                hk = as.integer(kmax), delta = as.double(delta), 
                ngrid = as.integer(ngrid), sbn = as.integer(sbn), 
                sbdry = as.integer(sbdry), tol = as.double(tol), 
                PACKAGE = "DNAcopy")
        }
        else {
            zzz <- list()
            zzz$ncpt <- 0
        }
        if (zzz$ncpt == 0) 
            change.loc <- c(change.loc, seg.end[k])
        seg.end <- switch(1 + zzz$ncpt, seg.end[-k], c(seg.end[1:(k - 
            1)], seg.end[k - 1] + zzz$icpt[1], seg.end[k]), c(seg.end[1:(k - 
            1)], seg.end[k - 1] + zzz$icpt, seg.end[k]))
        k <- length(seg.end)
        if (verbose >= 3) 
            cat(".... segments to go:", seg.end, "\n")
    }
    seg.ends <- rev(change.loc)
    nseg <- length(seg.ends)
    lseg <- diff(c(0, seg.ends))
    if (nseg > 1) {
        if (undo.splits == "prune") {
            lseg <- changepoints.prune(genomdat, lseg, undo.prune)
        }
        if (undo.splits == "sdundo") {
            lseg <- .changepoints.sdundo2(genomdat, lseg, trimmed.SD, 
                undo.SD, clen, relSDlong)
        }
    }
    segmeans <- 0 * lseg
    ll <- uu <- 0
    for (i in 1:length(lseg)) {
        uu <- uu + lseg[i]
        segmeans[i] <- mean(genomdat[(ll + 1):uu])
        ll <- uu
    }
    list(lseg = lseg, segmeans = segmeans)
}


.changepoints.sdundo2 <- function (genomdat, lseg, trimmed.SD, change.SD = 3, clen, relSDlong) 
{
    change.SD <- trimmed.SD * change.SD
    cpt.loc <- cumsum(lseg)
    sdundo <- TRUE
    while (sdundo) {
        k <- length(cpt.loc)
        if (k > 1) {
            segments0 <- cbind(c(1, 1 + cpt.loc[-k]), cpt.loc)
                segmed <- apply(segments0, 1, function(i, x) {
                median(x[i[1]:i[2]])}, genomdat)
            ens <- apply(segments0, 1, function(i, x) {i[2]+1-i[1]
            }, genomdat)
            el <- length(ens)
            ens2 <- c(2,ens)
            ens1 <- c(ens,2)
            allens <- apply(cbind(ens1,ens2),1,function(ro) ifelse(min(ro[1],ro[2])>clen,relSDlong,1))[-c(1,(el+1))] 
            
            adsegmed <- abs(diff(segmed)*allens)
            if (min(adsegmed) < change.SD) {
                i <- which(adsegmed == min(adsegmed))
                cpt.loc <- cpt.loc[-i]
            }
            else {
                sdundo <- FALSE
            }
        }
        else {
            sdundo <- FALSE
        }
    }
    lseg.sdundo <- diff(c(0, cpt.loc))
    lseg.sdundo
}

######################################################################


.assignNames <- function(matrix, object) {
    colnames(matrix) <- sampleNames(object)
    rownames(matrix) <- featureNames(object)
    matrix
}

.makeEmptyFeatureData <- function(object) {
    dims        <- Biobase:::assayDataDims(object)
    n           <- dims[1,1]
    features    <-         
    if (is(object, "environment")) ls(object)
    else names(object)
    nms         <- rownames(object[[features[[1]]]])
    data        <- data.frame(Chromosome=numeric(n), Start=numeric(n), End=numeric(n), row.names=nms)
    dimLabels   <- c("featureNames", "featureColumns")
    metadata    <- data.frame(labelDescription=c("Chromosomal position", "Basepair position start", "Basepair position end"), row.names=c("Chromosome", "Start", "End"))
    new("AnnotatedDataFrame", data=data, dimLabels=dimLabels, varMetadata=metadata)                  
}

.featureDataRequiredColumns <- function(featureData, columns) {
    msg     <- NULL
    absent  <- columns[!(columns %in% rownames(varMetadata(featureData)))]
    if (length(absent) != 0) {
        msg <- paste(msg, paste("missing columns' ", absent ,"' in featureData" , sep = "", collapse = "\n\t"), sep="\n")
    }
    if (is.null(msg)) TRUE else msg
}



.countcl <- function(k, regionsdat) {
    regionsdat[k,2]-regionsdat[k,1]+1
}

.sumreg <- function(k, dat, regionsdat) {
    sum(dat[regionsdat[k,1]:regionsdat[k,2]])
}

.sumsqreg <- function(k, dat, regionsdat) {
    sum((dat[regionsdat[k,1]:regionsdat[k,2]])^2)
}

.varregtimescount <- function(k, counts, dat, regionsdat) {
    var(dat[regionsdat[k,1]:regionsdat[k,2]])*counts[k]
}

.profreg <- function(k) {
    profile[k]
}

.varproffun <- function(prof, vcnmat, profile) {
    vcnprof <- vcnmat[profile == prof & !is.na(vcnmat[,1]),]
    if(!is.null(dim(vcnprof))) {
        return(sum(vcnprof[,1])/sum(vcnprof[,3]))
    } else {
        return(vcnprof[1]/vcnprof[3])
    }
}

.xgivenknotrunc <- function(k, class, pm, varprofall, allsum, allsumsq, allnc, robustsig, allcell,betas=c(1,0)) {

    sumx    <- allsum[k]
    sumxsq  <- allsumsq[k]
    ncl     <- allnc[k]
    v1      <- 2*varprofall[k]
    sd1     <- sqrt(v1)
    cell <- allcell[k]

    if (class==6) {
        mu <- (log2(2)/log2(1.5))*(0.10 + exp(-pm[4]))+0.3+exp(-pm[5])
    } else if (class==5) {
        mu <- (log2(2)/log2(1.5))*(0.10 + exp(-pm[4]))  
    } else if (class==4) {
        mu <- 0.10 + exp(-pm[class])
    } else if (class==3) {
        mu <- -0.05+0.1*exp(-(pm[class])^2) 
    } else if (class==1) {
        mu <- -0.10-exp(-pm[2])-0.3 - exp(-pm[1])
    } else {
        mu <- -0.10-exp(-pm[class])
    }
    
    munorm <- -0.05+0.1*exp(-(pm[class])^2) 
    #NEW
    mu <- ((1-cell)*munorm + cell*mu)*(betas[1]+betas[2]*sd1)
    
    v2 <- if(class == 4 | class == 2) 2*((pm[5+class])^2 + 0.0001) 
    else {if(class==3) {if(robustsig=="yes") 1/4*((pm[7])^2 + 0.0001)+1/4*((pm[9])^2 + 0.0001)+2*((pm[8]^2)+0.0001)
    else {2*((pm[8])^2 + 0.0001)}}
    else {if(class==1) 2*(pm[7]^2+(pm[6])^2 + 0.0001) 
    else {if(class==5) 2*((pm[10])^2+(pm[9])^2 + 0.0001)
    else 2*((pm[11])^2 + (pm[10])^2 + (pm[9])^2 + 0.0001)}}} # for stability
    
    v2norm <- if(robustsig=="yes") {1/4*((pm[7])^2 + 0.0001)+1/4*((pm[9])^2 + 0.0001)+2*((pm[8]^2)+0.0001)}
    else {2*((pm[8])^2 + 0.0001)}
    
    #NEW
    v2 <- if(class==3) v2norm else {(1-cell)^2*v2norm + cell^2*v2}
    #v2 <- (1-cell)*v2norm + cell^2*v2
    
    A       <- (sumx*v2+mu*v1)/(ncl*v2+v1)
    B       <- v1*v2/(ncl*v2+v1)
    C       <- sumxsq/v1 + mu^2/v2 - A^2/B
    logD    <- log(B)/2 -log(v2)/2
    logres  <- - C + logD
    res     <- exp(logres)
    return(logres)
}

.posteriorp <- function(k, priorp, pm, varprofall, allsum, allsumsq, allnc,robustsig,allcell) {
   all3        <- sapply(c(1,2,3,4,5,6), .xgivenknotrunc, k=k, pm=pm, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc,robustsig=robustsig,allcell=allcell)
    maxim       <- max(all3)
    all         <- exp(all3-maxim)
    allprior    <- all*priorp[k,]
    tot         <- all%*%priorp[k,]
    return(allprior/tot)
}

.reallikk4 <- function(k, alpha, pm, varprofall, allsum, allsumsq, allnc,robustsig,allcell,betas=c(1,0)) {
#betas=c(1,0);allcell=allcell_pr*0.2;k=1;pm=bstart; varprofall =varprof_allall_pr;allsum= allsumall_pr;allsumsq= allsumsqall_pr;allnc= allncall_pr;robustsig= robustsig;
    all3    <-sapply(c(1,2,3,4,5,6), .xgivenknotrunc, k=k, pm=pm, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, 
    allnc=allnc,robustsig=robustsig,allcell=allcell,betas=betas)
    maxim   <- max(all3)
    all     <- exp(all3-maxim)
    tot     <- maxim + log(all%*%alpha[k,])
    return(tot)
}

.reallik4 <- function(nreg, alpha, pm, varprofall, allsum, allsumsq, allnc,robustsig,allcell,betas=c(1,0)) {
    return(-sum(sapply(1:nreg, .reallikk4, pm=pm, alpha=alpha, varprofall =varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc,robustsig=robustsig,allcell=allcell,betas=betas)))
}

#.reallik4b <- function(nreg, alpha, pm, varprofall, allsum, allsumsq, allnc,robustsig,allcell,betas=c(1,0)) {
#    return(-sapply(1:nreg, .reallikk4, pm=pm, alpha=alpha, varprofall =varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc,robustsig=robustsig,allcell=allcell,betas=betas))
#}

.alphafun <- function(k, profile, priorp, pm, varprofall, allsum, allsumsq, allnc=allnc, robustsig, allcell, prior) {
    prof        <- profile[k]
    regionsk    <- which(profile==prof)
    nregk       <- length(regionsk)
    if (prior=="all") {
        nregk       <- length(profile)
        regionsk    <- 1:nregk
    }
    totpost     <- rep(1,nregk)%*%t(sapply(regionsk, .posteriorp, priorp=priorp, pm=pm, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc,robustsig=robustsig,allcell=allcell))
    return(totpost/nregk)
}

.alpha0all <- function(nreg, profile, priorp, pm, varprofall, allsum, allsumsq, allnc, robustsig, allcell, prior) {
    prevprof    <- 0
    alphaall    <- c()
    alpha1      <- .alphafun(1, profile, priorp, pm, varprofall, allsum, allsumsq, allnc, robustsig, allcell, prior)
    for (i in (1:nreg)) {
        if (prior != "all") {
            curprof <- profile[i]
            if (curprof==prevprof) {
                newalpha <- oldalpha
            }
            if (curprof != prevprof) {
                newalpha <- .alphafun(i, profile, priorp, pm, varprofall, allsum, allsumsq, allnc, robustsig, allcell, prior)
            }
            oldalpha <- newalpha
            prevprof <- curprof
        }
        if (prior == "all") {
            newalpha <- alpha1
        }
        alphaall <- rbind(alphaall,newalpha)
    }
    return(alphaall)
}


.minusEloglikreg <- function(k, posteriorprev, alphaprev, pm, varprofall, allsum, allsumsq, allnc, robustsig, allcell) {
    all3 <- sapply(c(1,2,3,4,5,6), .xgivenknotrunc, k=k, pm=pm, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc, robustsig=robustsig, allcell=allcell)
    return(-all3%*%posteriorprev[,k] - (log(alphaprev[k,])%*%posteriorprev[,k]))
} 
  
.totallik <- function(nreg, pm, posteriorprev, alphaprev, varprofall, allsum, allsumsq, allnc, robustsig,allcell,ncpus=1) {
    if(ncpus==1) {
    return(sum(sapply(1:nreg, .minusEloglikreg, pm=pm, posteriorprev=posteriorprev, alphaprev=alphaprev, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc,robustsig=robustsig,allcell=allcell)))
    } else {
    sfExport("pm")
    return(sum(sfSapply(1:nreg, .minusEloglikreg, pm=pm, posteriorprev=posteriorprev, alphaprev=alphaprev, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc,robustsig=robustsig,allcell=allcell)))
    } 
}

#.mwp <- function(prof) {
#    meanpost <- cbind(allmean,posteriorfin2)
#    return(meanpost[profile==prof,])
#}

.MakeData <- function(smratall,chrnum) {  #updated 16/07/10
    #smratall <- datmat[,-(1:nc)];chrnum<-chr;
    nc <- ncol(smratall)
    allregions  <- vector("list", nc)
    ls      <- length(chrnum)
    chrnumfw <- c(100,chrnum[-ls])
    for (j in 1:(nc)) {
        smrat   <- smratall[,j] #changed 19/06/2009, much faster
        smratshfw <- c(100,smrat[-ls])
        smratshbw <- c(smrat[-1],100)
        mult <- (smrat-smratshfw)*(smratshbw-smratshfw) + (chrnumfw-chrnum)
        wm <- which(mult!=0)
        smwh <- smrat[wm]
        wmend <- c(wm[-1]-1,ls)
        regions <- cbind(wm,wmend,smwh)
        allregions[[j]] <- regions
#  breaks  <- c(0)
#        for(i in 2:(ls-1)) {
#            breaks <- c(breaks, ifelse(((smrat[i] != smrat[i-1]) & (smrat[i-1] != smrat[i+1])) | chrnum[i] != chrnum[i-1], 1, 0))
#        } 
#        breaks  <- c(breaks,0)
#        ind     <- 1:ls
#        difind  <- cbind(ind, breaks)
#        br      <- difind[breaks==1,2]
#        if (length(br) < 1) breakpos <- c()
#        if (length(br) ==1) breakpos <- difind[breaks == 1,][1]
#        if (length(br) > 1) breakpos <- difind[breaks == 1,][,1]
#        regions <- cbind(append(breakpos, 1, 0), append(breakpos-1, ls, length(breakpos)))
    rm(smrat,smratshfw,smratshbw,mult,smwh,wm);gc()
    }
    return(allregions)
}
