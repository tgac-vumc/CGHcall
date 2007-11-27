CGHcall <- function(inputSegmented, prior="auto", nclass=3, organism="human") {

    ## Version 2.0
    ## Author: Mark van de Wiel
    ## Maintainer: Sjoerd Vosse
    ## Email: s.vosse@vumc.nl
    ## Date: 24-10-2007
    
    ## Changes since previous version
    ## - input is of class cgh

    timeStarted <- proc.time()

    ## Read input data
    normalizedData  <- copynumber(inputSegmented)
    segmentedData   <- segmented(inputSegmented)
    combined        <- data.frame(featureNames(inputSegmented), chromosomes(inputSegmented), bpstart(inputSegmented), normalizedData, segmentedData)
    
    datareg     <- CGHcall:::.MakeData(combined)
    posit       <- datareg[[1]][,3]
    chr         <- datareg[[1]][,2]
    naam        <- datareg[[1]][,1]
    dataprob    <- data.frame(naam, chr, posit)
    nc          <- (ncol(datareg[[1]])-3)/2
    
    ## Determine method for prior probabilities
    if (prior == 'auto') {
        if (nc >  20) prior <- "not all";
        if (nc <= 20) prior <- "all";
    }
    
    ### Convert to chromosome arm if neccessary
    if (prior != "all" && organism == "human") {
        temp            <- datareg[[1]];
        datareg[[1]]    <- CGHcall:::.convertChromosomeToArm(temp);
    }

    cat("EM algorithm started ... \n");
    
    dat         <- c()
    datsm       <- c()
    regions     <- c()
    regionsdat  <- c()
    nclones     <- length(chr)
    profile     <- c()
    
    for (j in (1:nc)) {
        dat         <- c(dat, datareg[[1]][,(3+j)])
        datsm       <- c(datsm, datareg[[1]][,(3+nc+j)])
        regions1    <- datareg[[2]][[j]]
        nreg1       <- nrow(regions1)
        profile     <- c(profile, rep(j,nreg1))
        regions     <- rbind(regions, regions1)
        toadd       <- (j-1)*c(nclones, nclones)
        regionsdat  <- rbind(regionsdat, regions1+toadd)
    }

    takechr     <- function(reg, chrom) {
        chrom[reg[1]]
    }
    
    chrarmreg   <- as.vector(apply(as.matrix(regions), 1, takechr, chrom=chr))

    outthr          <- 2
    toolarge        <- which(abs(datsm)>=2)
    datold          <- dat
    dat[toolarge]   <- 0

    nreg        <- nrow(regionsdat)
    allnc       <- sapply(1:nreg, CGHcall:::.countcl, regionsdat=regionsdat)
    allsum      <- sapply(1:nreg, CGHcall:::.sumreg, dat=dat, regionsdat=regionsdat)
    allsumsq    <- sapply(1:nreg, CGHcall:::.sumsqreg, dat=dat, regionsdat=regionsdat)
    varregall   <- sapply(1:nreg, CGHcall:::.varregtimescount, counts=allnc, dat=dat, regionsdat=regionsdat)
    lev         <- 1:nc
    varprofnc   <- cbind(varregall, profile, allnc)
    varprof     <- sapply(lev, CGHcall:::.varproffun, vcnmat=varprofnc, profile=profile)

    selk <- function(k, varprof) {
        profk <- profile[k]
        return(max(0.001,varprof[profk]))
    }
    varprofall <- sapply(1:nreg, selk, varprof=varprof)

    allmean     <- allsum/allnc
    allmeanl    <- allmean[allmean < -0.15 & allmean > -0.7]
    allmeang    <- allmean[allmean>0.15 & allmean<0.7]
    allmean0    <- allmean[abs(allmean) <= 0.15]
    sdlst       <- if(length(allmeanl) <= 1) 0.01 else mad(allmeanl)
    meanlst     <- if(length(allmeanl) <= 0) -0.3 else mean(allmeanl)
    sdgst       <- if(length(allmeang) <= 1) 0.01 else mad(allmeang)
    meangst     <- if(length(allmeang) <= 0) 0.3 else mean(allmeang)
    sd0st       <- if(length(allmean0) <= 1) 0.01 else mad(allmean0)

    profchrom   <- chrarmreg
    bstart      <- c(-log(0.2), -log(-(meanlst+0.1)), sqrt(-log(0.5)), -log(meangst-0.1), -log(0.2), 0.0001, sdlst, sd0st, sdgst, 0.0001)
    mus         <- c(-0.10-exp(-bstart[2])-0.3 - exp(-bstart[1]), -0.10-exp(-bstart[2]), -0.05+0.1*exp(-(bstart[3])^2), 0.1 + exp(-bstart[4]), 2*(0.10 + exp(-bstart[4]))+0.3+exp(-bstart[5]))
    priorp      <- matrix(rep(c(0.01, 0.09, 0.8, 0.08, 0.01, 0.01), nreg), ncol=6, byrow=TRUE)
    alpha0      <- CGHcall:::.alpha0all(nreg, profchrom, priorp, bstart, varprofall, allsum, allsumsq, allnc, prior)

    maxiter     <- 10
    stop        <- 0
    iter        <- 1
    thrpara     <- 0.01
    
    while (stop == 0 & iter <= maxiter) {
        posterior0  <- sapply(1:nreg, CGHcall:::.posteriorp, priorp=alpha0, pm=bstart, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc)
        likprev     <- CGHcall:::.totallik(bstart, nreg=nreg, posteriorprev=posterior0, alphaprev=alpha0, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc)
        nlmres      <- nlm(CGHcall:::.totallik, bstart, nreg=nreg, posteriorprev=posterior0, alphaprev=alpha0, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc)
        bprev       <- bstart
        bstart      <- nlmres$est
        alpha0      <- CGHcall:::.alpha0all(nreg, profchrom, alpha0, bstart, varprofall, allsum, allsumsq, allnc, prior)
        rl          <- CGHcall:::.reallik4(nreg, alpha0, bstart, varprofall, allsum, allsumsq, allnc)
        llmin       <- nlmres$min
        musprev     <- c(-0.10-exp(-bprev[2])-0.3 - exp(-bprev[1]), -0.10-exp(-bprev[2]), -0.05+0.1*exp(-(bprev[3])^2), 0.1 + exp(-bprev[4]), 2*(0.10 + exp(-bprev[4]))+0.3+exp(-bprev[5]))
        mus         <- c(-0.10-exp(-bstart[2])-0.3 - exp(-bstart[1]), -0.10-exp(-bstart[2]), -0.05+0.1*exp(-(bstart[3])^2), 0.1 + exp(-bstart[4]),2*(0.10 + exp(-bstart[4]))+0.3+exp(-bstart[5]))
        param       <- c(mus, bstart[6], bstart[7], bstart[8], bstart[9], bstart[10])
        paramprev   <- c(musprev, bprev[6], bprev[7], bprev[8], bprev[9], bprev[10])
        cat("Calling iteration", iter, ":\n")
        print(c(j,rl,mus,bstart[6],bstart[7],bstart[8],bstart[9],bstart[10]))
        if (max(abs(param[2:4]-paramprev[2:4])) <= thrpara) stop <- 1
        iter <- iter+1
    }
    
    cat("EM algorithm done ...\n")

    best            <- bstart
    allsumold       <- sapply(1:nreg, CGHcall:::.sumreg, dat=datold, regionsdat=regionsdat)
    allsumsqold     <- sapply(1:nreg, CGHcall:::.sumsqreg, dat=datold, regionsdat=regionsdat)
    posteriorfin    <- t(sapply(1:nreg, CGHcall:::.posteriorp, priorp=alpha0, pm=best, varprofall=varprofall, allsum=allsumold, allsumsq=allsumsqold, allnc=allnc))
    posteriorfin1   <- cbind(allnc, allmean, posteriorfin)
    amps            <- posteriorfin1[posteriorfin1[,7]+posteriorfin1[,8]>=0.5,]

    amporgain <- function(row) {
        if (row[6] >= 0.5) return(c(row[1]+row[2], row[3], row[4]+row[5], row[6]))
        else return(c(row[1]+row[2], row[3], row[4]+row[5], row[6]))
    }

    if (nclass == 3) posteriorfin1 <- cbind(posteriorfin[,1]+posteriorfin[,2], posteriorfin[,3], posteriorfin[,4]+posteriorfin[,5]+posteriorfin[,6], posteriorfin[,6])
    if (nclass == 4) posteriorfin1 <- t(apply(posteriorfin, 1, amporgain))

    posteriorfin2   <- cbind(profile, posteriorfin1)
    regionsprof     <- cbind(profile, regions)
    dataprob        <- data.frame(naam, chr, posit)

    for (k in (1:nc)) {
        post        <- (posteriorfin2[profile==k,])[,-1]
        regionsk    <- (regionsprof[profile==k,])[,-1]
        nregk       <- nrow(post)
        probs       <- c()
        for (i in (1:nregk)) {
            regl    <- regionsk[i,2]-regionsk[i,1]+1
            togeth  <- post[i,(1:nclass)]
            probs   <- c(probs,rep(togeth,regl))
        }
        allprobs    <- matrix(probs, ncol=nclass, byrow=TRUE)
        datk        <- datareg[[1]][,(3+k)]
        dataprob    <- data.frame(dataprob, datk, allprobs)
    }
    
    ### Add column names
    info        <- c("Name", "Chromosome", "Position")
    samplenames <- sampleNames(inputSegmented)
    col.names   <- rep(samplenames, each=nclass+1)
    
    paste.stuff <- c("log2ratio", "loss", "normal", "gain")
    if (nclass == 4) paste.stuff <- c(paste.stuff, "amplification")
    
    col.names   <- paste(paste.stuff, col.names, sep="_")
    col.names   <- c(info, col.names)
    
    colnames(dataprob) <- col.names
    
    ### Define calls from probabilities
    ncpp        <- nclass+1
    nclone      <- nrow(dataprob)
    ncolscl     <- ncol(normalizedData)
    classify.res        <- array(NA,c(nclone,ncolscl))
    for (i in 1:nc) {
        genomdat        <- dataprob[,(3+ncpp*(i-1)+1)]
        prob.loss.ind   <- dataprob[,(3+ncpp*(i-1)+2)]
        prob.none.ind   <- dataprob[,(3+ncpp*(i-1)+3)]
        prob.gain.ind   <- dataprob[,(3+ncpp*(i-1)+4)]
        if (nclass==4) {
            prob.amp.ind    <- dataprob[,(3+ncpp*(i-1)+5)]
            ticksamp        <- which(prob.amp.ind >= 0.5)
            lt              <- length(ticksamp)
        }        
        
        if (nclass==3) {
            classify.res[(1:nclone),i] <- (2*(as.numeric(prob.gain.ind>prob.loss.ind))-1)*(as.numeric(apply(cbind(prob.loss.ind,prob.gain.ind),1,max)>prob.none.ind))
        }
        
        if (nclass==4) {
            prob.gain.amp <- prob.gain.ind + prob.amp.ind
            classify.res[(1:nclone),i] <- (as.numeric(prob.amp.ind>0.5)+1)*(2*(as.numeric(prob.gain.amp>prob.loss.ind))-1)*(as.numeric(apply(cbind(prob.loss.ind,prob.gain.ind),1,max)>prob.none.ind))
        }              
    }
    
    ### Make result list
    
    splitProbs <- function(probs, nsample=nc, classes=nclass) {
        probs   <- probs[,4:ncol(probs)]
        count   <- classes+1
        probloss <- c()
        probnorm <- c()
        probgain <- c()
        probamp <- c()
        for (i in 0:(nsample-1)) {
            probloss <- cbind(probloss, probs[,i*count+2])
            probnorm <- cbind(probnorm, probs[,i*count+3])
            probgain <- cbind(probgain, probs[,i*count+4])
            if (nclass==4) probamp <- cbind(probamp, probs[,i*count+5])
        }
        if (nclass == 3) result  <- list(loss=probloss, normal=probnorm, gain=probgain)
        else if (nclass == 4) result  <- list(loss=probloss, normal=probnorm, gain=probgain, amp=probamp)
        result
    }
    
    allmeanold  <- allsumold/allnc  #include amplifications
    profmeannc  <- cbind(profile, allmeanold, regions)
    
    print(head(dataprob))
    
    probs       <- splitProbs(dataprob)
    
    print(probs$amp[1:5,])
    
    if (nclass == 3) assayData <- assayDataNew( copynumber  = copynumber(inputSegmented),
                                                segmented   = segmented(inputSegmented),
                                                calls       = CGHcall:::.assignNames(classify.res, inputSegmented), 
                                                probloss    = CGHcall:::.assignNames(probs$loss, inputSegmented), 
                                                probnorm    = CGHcall:::.assignNames(probs$normal, inputSegmented), 
                                                probgain    = CGHcall:::.assignNames(probs$gain, inputSegmented)
                                              )
    else if (nclass == 4) assayData <- assayDataNew(copynumber  = copynumber(inputSegmented),
                                                    segmented   = segmented(inputSegmented),
                                                    calls       = CGHcall:::.assignNames(classify.res, inputSegmented), 
                                                    probloss    = CGHcall:::.assignNames(probs$loss, inputSegmented), 
                                                    probnorm    = CGHcall:::.assignNames(probs$normal, inputSegmented), 
                                                    probgain    = CGHcall:::.assignNames(probs$gain, inputSegmented),
                                                    probamp     = CGHcall:::.assignNames(probs$amp, inputSegmented)
                                                   )
                                                            
    result  <- CGHcall:::.callFromSeg(inputSegmented, assayData)
    
    cat("FINISHED!\n")
    timeFinished <- round((proc.time() - timeStarted)[1] / 60)
    cat("Total time:", timeFinished, "minutes\n")
    result
}
