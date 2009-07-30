CGHcall <- function(inputSegmented, prior="auto", nclass=3, organism="human", robustsig="yes",digits = 3,nsegfit=3000,maxnumseg=100,minlsforfit=0.5) {
    #load("seg.Rdata");inputSegmented <- seg[,1:5]; prior="auto"; nclass=4; organism="human"; robustsig="yes";digits=3;nsegfit=100;maxnumseg=100;minlsforfit=0.5
    print("changed")
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
    bw <- round(20*minlsforfit*nrow(normalizedData)/44000) 
    
    datareg     <- CGHcall:::.MakeData(combined)            
    #datareg <- .MakeData(combined)
    #save(datareg,file="datareg.Rdata")
    #load("datareg.Rdata")
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
       # datareg[[1]]    <- .convertChromosomeToArm(temp);
        chr <- datareg[[1]][,2] #BUG;repaired 22/06/09
    }

    cat("EM algorithm started ... \n");
    
    dat         <- c()
    datsm       <- c()
    datall         <- c()
    datsmall       <- c()
    regions     <- c()
    regionsall  <- c()
    regionsdat  <- c()
    regionsdatall  <- c()
    nclones     <- length(chr)
     profile     <- c()
     profileall <- c()
    countreg    <- 0
    newregions <- list()
#    nclonesj <- 0
    thr <- 2
   
    #here is the selection of regions used for fitting. First select regions long enough and with smoothed signal <= thr. Then, maximize #regions per prof to maxseg
    #nc<-2
    for (j in (1:nc)) {
        datall        <- c(datall, datareg[[1]][,(3+j)])
        datsmall       <- c(datsmall, datareg[[1]][,(3+nc+j)])
        regions1    <- datareg[[2]][[j]]
        nreg1       <- nrow(regions1)
        profileall     <- c(profileall, rep(j,nreg1))
        regionsall     <- rbind(regionsall, regions1)
        toaddall       <- (j-1)*c(nclones, nclones)
        regionsdatall  <- rbind(regionsdatall, regions1[,-3]+toaddall)

        segl        <- regions1[,2]-regions1[,1]+1
        whseg       <- which(segl>=bw)
        regions2    <- regions1[whseg,] #selects regions long enough to participate in fitting
        regions2    <- regions2[which(abs(regions2[,3]) <= thr),] #selects regions with smoothed signal small enough
        nreg2       <- nrow(regions2)
        countreg <- countreg + min(maxnumseg,nreg2)
        newregions[[j]] <- regions2
        }
        fracforfit <- min(1,nsegfit/countreg)
        
        for (j in (1:nc)) {
        regions2 <- newregions[[j]]
        nreg2       <- nrow(regions2)
        ordsig      <- order(regions2[,3])
        regions2b   <- data.frame(regions2,genord = 1:nreg2)
        regions2c   <- regions2b[ordsig,]
        if(nreg2 >= maxnumseg | fracforfit<1){
            takeseg     <- floor(seq(1,nreg2,length.out=round(min(nreg2,maxnumseg)*fracforfit)))
            regions2c <- regions2c[takeseg,]
            }
        regions2d <- regions2c[order(regions2c[,4]),]
        
        nreg2       <- nrow(regions2d)
        profile    <- c(profile, rep(j,nreg2))
        regions     <- rbind(regions, regions2d)
        #takerow <- c();for(k in 1:nrow(regions2d)) {takerow <- c(takerow,regions2d[k,1]:regions2d[k,2])}
#        dat <- c(dat,datareg[[1]][takerow,(3+j)])
#        datsm <- c(datsm,datareg[[1]][takerow,(3+nc+j)])
        toadd       <- (j-1)*c(nclones, nclones)
        regionsdat  <- rbind(regionsdat, regions2d[,-(3:4)]+toadd) #delete smrat value and index
        #        nclonesj <- sum((regions2d[,2]+1)-regions2d[,1]) 
    }

    takechr     <- function(reg, chrom) {
        chrom[reg[1]]
    }
    
    chrarmreg   <- as.vector(apply(as.matrix(regions), 1, takechr, chrom=chr))
    
    nreg        <- nrow(regionsdat)
    nregall     <- nrow(regionsdatall)
    print(paste("Total number of segments present in the data:",nregall))
    print(paste("Number of segments used for fitting the model:",nreg))
    allnc       <- sapply(1:nreg, CGHcall:::.countcl, regionsdat=regionsdat)
    allsum      <- sapply(1:nreg, CGHcall:::.sumreg, dat=datall, regionsdat=regionsdat)
    allsumsq    <- sapply(1:nreg, CGHcall:::.sumsqreg, dat=datall, regionsdat=regionsdat)
    varregall   <- sapply(1:nreg, CGHcall:::.varregtimescount, counts=allnc, dat=datall, regionsdat=regionsdat)
    lev         <- 1:nc
    varprofnc   <- cbind(varregall, profile, allnc)
    varprof     <- sapply(lev, CGHcall:::.varproffun, vcnmat=varprofnc, profile=profile)

    selk <- function(k, varprof,prof) {  #changed 22/06/09
        profk <- prof[k]
        return(max(0.001,varprof[profk]))
    }
    varprofall <- sapply(1:nreg, selk, varprof=varprof,prof=profile)

    allmean     <- allsum/allnc
    allmeanl    <- allmean[allmean < -0.15 & allmean > -0.7]
    allmeang    <- allmean[allmean>0.15 & allmean<0.7]
    allmean0    <- allmean[abs(allmean) <= 0.15]
    sdlst       <- if(length(allmeanl) <= 1) 0.01 else mad(allmeanl)
    meanlst     <- if(length(allmeanl) <= 0) -0.3 else mean(allmeanl)
    sdgst       <- if(length(allmeang) <= 1) 0.01 else mad(allmeang)
    meangst     <- if(length(allmeang) <= 0) 0.3 else mean(allmeang)
    sd0st       <- if(length(allmean0) <= 1) 0.01 else mad(allmean0)
    #sd0st       <- if(length(allmean0) <= 1) 0.01 else sd(allmean0)

    profchrom   <- chrarmreg
    if(robustsig=="no") {
        bstart      <- c(-log(0.2), -log(-(meanlst+0.1)), sqrt(-log(0.5)), -log(meangst-0.1), -log(0.2), 0.0001, sdlst, sd0st, sdgst, 0.0001, 0.0001)
    } else {
        bstart      <- c(-log(0.2), -log(-(meanlst+0.1)), sqrt(-log(0.5)), -log(meangst-0.1), -log(0.2), 0.0001, sdlst, sd0st+max(0,sqrt(sdgst^2/8+sdlst^2/8)-sd0st), sdgst, 0.0001, 0.0001) #robust option 
    }
    mus         <- c(-0.10-exp(-bstart[2])-0.3 - exp(-bstart[1]), -0.10-exp(-bstart[2]), -0.05+0.1*exp(-(bstart[3])^2), 0.1 + exp(-bstart[4]), 2*(0.10 + exp(-bstart[4]))+0.3+exp(-bstart[5]))
    priorp      <- matrix(rep(c(0.01, 0.09, 0.8, 0.08, 0.01, 0.01), nreg), ncol=6, byrow=TRUE)
    alpha0      <- CGHcall:::.alpha0all(nreg, profchrom, priorp, bstart, varprofall, allsum, allsumsq, allnc,robustsig,prior)

    maxiter     <- 10
    stop        <- 0
    iter        <- 1
    thrpara     <- 0.015 #changed 15/6/09
    
    while (stop == 0 & iter <= maxiter) {
        cat("Calling iteration", iter, ":\n")
        posterior0  <- sapply(1:nreg, CGHcall:::.posteriorp, priorp=alpha0, pm=bstart, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc,robustsig=robustsig)
        likprev     <- CGHcall:::.totallik(bstart, nreg=nreg, posteriorprev=posterior0, alphaprev=alpha0, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc,robustsig=robustsig)
        nlmres      <- nlm(CGHcall:::.totallik, bstart, nreg=nreg, posteriorprev=posterior0, alphaprev=alpha0, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc,robustsig=robustsig)
        bprev       <- bstart
        bstart      <- nlmres$est
        alpha0      <- CGHcall:::.alpha0all(nreg, profchrom, alpha0, bstart, varprofall, allsum, allsumsq, allnc, robustsig, prior)
        rl          <- CGHcall:::.reallik4(nreg, alpha0, bstart, varprofall, allsum, allsumsq, allnc, robustsig)
        llmin       <- nlmres$min
        musprev     <- c(-0.10-exp(-bprev[2])-0.3 - exp(-bprev[1]), -0.10-exp(-bprev[2]), -0.05+0.1*exp(-(bprev[3])^2), 0.1 + exp(-bprev[4]),(log2(2)/log2(1.5))*(0.1 + exp(-bprev[4])), (log2(2)/log2(1.5))*(0.10 + exp(-bprev[4]))+0.3+exp(-bprev[5]))
        mus         <- c(-0.10-exp(-bstart[2])-0.3 - exp(-bstart[1]), -0.10-exp(-bstart[2]), -0.05+0.1*exp(-(bstart[3])^2), 0.1 + exp(-bstart[4]),(log2(2)/log2(1.5))*(0.1 + exp(-bstart[4])),(log2(2)/log2(1.5))*(0.10 + exp(-bstart[4]))+0.3+exp(-bstart[5]))
        param       <- c(mus, bstart[6], bstart[7], bstart[8], bstart[9], bstart[10],bstart[11])
        paramprev   <- c(musprev, bprev[6], bprev[7], bprev[8], bprev[9], bprev[10],bprev[11])
     
        if(robustsig=="no") {
            printmat <- matrix(c(j,rl,mus,sqrt(bstart[7]^2+(bstart[6])^2 + 0.0001),bstart[7],bstart[8],bstart[9],sqrt((bstart[10])^2+(bstart[9])^2 + 0.0001),sqrt((bstart[11])^2 + (bstart[10])^2+ (bstart[9])^2+ 0.0001)),nrow=1)
            colnames(printmat) <- c("j","rl","mudl","musl","mun","mug","mudg","mua","sddl","sdsl","sdn","sdg","sddg","sda")
            print(printmat)
        } else {
            printmat<-matrix(c(j,rl,mus,sqrt(bstart[7]^2+(bstart[6])^2 + 0.0001),bstart[7],sqrt(1/8*((bstart[7])^2 + 0.0001)+1/8*((bstart[9])^2 + 0.0001)+((bstart[8]^2)+0.0001)),bstart[9],sqrt((bstart[10])^2+(bstart[9])^2 + 0.0001),sqrt((bstart[11])^2 + (bstart[10])^2 + (bstart[9])^2 + 0.0001)),nrow=1)  #robust option
            colnames(printmat) <- c("j","rl","mudl","musl","mun","mug","mudg","mua","sddl","sdsl","sdn","sdg","sddg","sda")
            print(printmat)
        }
        if (max(abs(param[c(2:4,8)]-paramprev[c(2:4,8)])) <= thrpara) stop <- 1 #changed 15/6/09; also sdn should converge
        iter <- iter+1
    }
    
    cat("EM algorithm done ...\n")
    best            <- bstart
    
    #now start computing posteriors for ALL data
    cat("Computing posterior probabilities for all segments ...\n")
    #nregall        <- nrow(regionsdatall)
    allncall       <- sapply(1:nregall, CGHcall:::.countcl, regionsdat=regionsdatall)
    allsumall      <- sapply(1:nregall, CGHcall:::.sumreg, dat=datall, regionsdat=regionsdatall)
    allmeanall     <- allsumall/allncall  
    allsumsqall    <- sapply(1:nregall, CGHcall:::.sumsqreg, dat=datall, regionsdat=regionsdatall)
    varregallall   <- sapply(1:nregall, CGHcall:::.varregtimescount, counts=allncall, dat=datall, regionsdat=regionsdatall)
    lev             <- 1:nc
    varprofncall   <- cbind(varregallall, profileall, allncall)
    varprof_all     <- sapply(lev, CGHcall:::.varproffun, vcnmat=varprofncall, profile=profileall)
    varprof_allall <- sapply(1:nregall, selk, varprof=varprof_all,prof=profileall)

    
    if(prior!="all"){
    chrarm_alpha0 <- cbind(chrarmreg,alpha0)
    chrarm_alpha0uni <- unique(chrarm_alpha0) 
    chrarmregall   <- as.vector(apply(as.matrix(regionsall), 1, takechr, chrom=chr))
    alpha0_all      <- t(sapply(chrarmregall,function(x){chrarm_alpha0uni[chrarm_alpha0uni[,1]==x,-1]}))
    } else alpha0_all <- matrix(rep(alpha0[1,],nregall),byrow=T,nrow=nregall)
    
    posteriorfin    <- t(sapply(1:nregall, CGHcall:::.posteriorp, priorp=alpha0_all, pm=best, varprofall=varprof_allall, allsum=allsumall, allsumsq=allsumsqall, allnc=allncall, robustsig=robustsig))
    posteriorfin1   <- cbind(allncall, allmeanall, posteriorfin)
    amps            <- posteriorfin1[posteriorfin1[,7]+posteriorfin1[,8]>=0.5,]

    amporgain <- function(row) {
        if (row[6] >= 0.5) return(c(row[1]+row[2], row[3], row[4]+row[5], row[6]))
        else return(c(row[1]+row[2], row[3], row[4]+row[5], row[6]))
    }

    #if(length(toolargeseg) >= 0) {for(i in toolargeseg){
#    if(allmeanold[i] < 0) {posteriorfin[i,] <- c(1,0,0,0,0,0)} else {posteriorfin[i,] <- c(0,0,0,0,0,1)}}}
#    
    if (nclass == 3) posteriorfin1 <- cbind(posteriorfin[,1]+posteriorfin[,2], posteriorfin[,3], posteriorfin[,4]+posteriorfin[,5]+posteriorfin[,6], posteriorfin[,6])
    if (nclass == 4) posteriorfin1 <- t(apply(posteriorfin, 1, amporgain))

    posteriorfin2   <- cbind(profileall, posteriorfin1)
    regionsprof     <- cbind(profileall, regionsall)
    dataprob        <- data.frame(naam, chr, posit)

    for (k in (1:nc)) {
        post        <- (posteriorfin2[profileall==k,])[,-1]
        regionsk    <- (regionsprof[profileall==k,])[,-1]
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
    
    #allmeanold  <- allsumold/allnc  #include amplifications
#    profmeannc  <- cbind(profile, allmeanold, regions)
#    
    probs       <- splitProbs(dataprob)
    
    if (nclass == 3) assayData <- assayDataNew( copynumber  = round(copynumber(inputSegmented),digits),
                                                segmented   = round(segmented(inputSegmented),digits),
                                                calls       = CGHcall:::.assignNames(classify.res, inputSegmented), 
                                                probloss    = CGHcall:::.assignNames(round(probs$loss,digits), inputSegmented), 
                                                probnorm    = CGHcall:::.assignNames(round(probs$normal,digits), inputSegmented), 
                                                probgain    = CGHcall:::.assignNames(round(probs$gain,digits), inputSegmented)
                                                ) else if (nclass == 4) assayData <- assayDataNew(copynumber  = round(copynumber(inputSegmented),digits),
                                                    segmented   = round(segmented(inputSegmented),digits),
                                                    calls       = CGHcall:::.assignNames(classify.res, inputSegmented), 
                                                    probloss    = CGHcall:::.assignNames(round(probs$loss,digits), inputSegmented), 
                                                    probnorm    = CGHcall:::.assignNames(round(probs$normal,digits), inputSegmented), 
                                                    probgain    = CGHcall:::.assignNames(round(probs$gain,digits), inputSegmented),
                                                    probamp     = CGHcall:::.assignNames(round(probs$amp,digits), inputSegmented)
                                                    )
                                                            
    result  <- CGHcall:::.callFromSeg(inputSegmented, assayData)                                                            
    
    cat("FINISHED!\n")
    timeFinished <- round((proc.time() - timeStarted)[1] / 60)
    cat("Total time:", timeFinished, "minutes\n")
    result
}
