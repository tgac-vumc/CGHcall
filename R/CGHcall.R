CGHcall <- function(inputSegmented, prior="auto", nclass=4, organism="human", robustsig="yes",nsegfit=3000,maxnumseg=100,minlsforfit=0.5) {
  #  library(CGHcall)
#    setwd("C:\\VUData\\Maria")
#    #rm(list=ls())
#    memory.limit(size = 4000)
#    load("seg.Rdata")
#    inputSegmented <- seg; prior="auto"; nclass=4; organism="human"; robustsig="yes";nsegfit=500;maxnumseg=100;minlsforfit=0.5
    timeStarted <- proc.time()
    print("changed")
    gc() #memory usage in Mb
    
    ## Version 2.0
    ## Author: Mark van de Wiel
    ## Maintainer: Sjoerd Vosse
    ## Email: s.vosse@vumc.nl
    ## Date: 24-10-2007
    
    ## Changes since previous version
    ## - input is of class cgh
   
    timeStarted <- proc.time()

    ## Read input data
    #normalizedData  <- copynumber(inputSegmented)
    #combined        <- data.frame(featureNames(inputSegmented), chromosomes(inputSegmented), bpstart(inputSegmented), copynumber(inputSegmented),
    #segmented(inputSegmented))
    bw <- round(20*minlsforfit*nrow(copynumber(inputSegmented))/44000) 
    ncolscl     <- ncol(copynumber(inputSegmented))
  #  datareg     <- CGHcall:::.MakeData(data.frame(featureNames(inputSegmented), chromosomes(inputSegmented), bpstart(inputSegmented), copynumber(inputSegmented), segmented(inputSegmented)))  
    whna <- !is.na(chromosomes(inputSegmented)) & !is.na(featureNames(inputSegmented))
    datmat <- cbind(copynumber(inputSegmented), segmented(inputSegmented))[whna,]
    posit       <- bpstart(inputSegmented)[whna]
    chr         <- chromosomes(inputSegmented)[whna]
    naam        <- featureNames(inputSegmented)[whna]
    rm(whna,inputSegmented);gc()
    nc    <- ncol(datmat)/2 
    nclone <-nrow(datmat)
    datareg  <- .MakeData(datmat[,-(1:nc)],chr)  
   
    #rm(combined); 
    
    gc() #added 24/11/2009
              
    #datareg <- .MakeData(combined)
    #save(datareg,file="datareg.Rdata")
    #load("datareg.Rdata")
     dataprob    <- data.frame(naam, chr, posit)
    
 
    
    ## Determine method for prior probabilities
    if (prior == 'auto') {
        if (nc >  20) prior <- "not all";
        if (nc <= 20) prior <- "all";
    }
    
    ### Convert to chromosome arm if neccessary
    if (prior != "all" && organism == "human") {
       # temp            <- datareg[[1]];
        dataprob    <- CGHcall:::.convertChromosomeToArm(dataprob); 
       # datareg[[1]]    <- .convertChromosomeToArm(temp);
        chr <- dataprob[,2] #BUG;repaired 22/06/09
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
   
    #here is the selection of regions used for fitting. First select regions long enough and with smoothed signal <= thr. Then, maximize 
    #regions per prof to maxseg
    #nc<-2
    gc() 
    datall <- as.vector(datmat[,1:nc]) #changed 16/7/10
    datsmall      <- as.vector(datmat[,-(1:nc)])
    rm(datmat);gc();
  
    for (j in (1:nc)) {
        regions1    <- datareg[[j]]
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
    gc()

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
    
    gc()

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
    
    gc()
    
    while (stop == 0 & iter <= maxiter) {
        print(gc())
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
   # pm <- as.vector(printmat)
#    plot(function(x) dnorm(x,mean=pm[4],sd=pm[11]),-1,1,col="red")
#   plot(function(x) dnorm(x,mean=pm[5],sd=pm[12]),-1,1,col="black",add=TRUE)
#   plot(function(x) dnorm(x,mean=pm[6],sd=pm[13]),-1,1,col="green",add=TRUE)
#    cat("EM algorithm done ...\n")
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
    gc()
    
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
   
    if (nclass == 3) posteriorfin1 <- cbind(posteriorfin[,1]+posteriorfin[,2], posteriorfin[,3], posteriorfin[,4]+posteriorfin[,5]+posteriorfin[,6], posteriorfin[,6])
    if (nclass == 4) posteriorfin1 <- t(apply(posteriorfin, 1, amporgain))

    posteriorfin2   <- cbind(profileall, posteriorfin1)
    rm(posteriorfin,posteriorfin1);gc() #added24/11/2009
    regionsprof     <- cbind(profileall, regionsall)
    listcall <- list(posteriorfin2,nclone,nc,nclass,regionsprof)
    #save(listcall,file="listcall.Rdata")
    gc()
    timeFinished <- round((proc.time() - timeStarted)[1] / 60)
    cat("Total time:", timeFinished, "minutes\n")
    return(listcall)
    }
    

ExpandCGHcall <- function(listcall,inputSegmented, digits=3,divide=4, memeff = FALSE, fileoutpre="Callobj_"){
  timeStarted <- proc.time()
  posteriorfin2 <- listcall[[1]];nclone<-listcall[[2]];nctot <- listcall[[3]];nclass <- listcall[[4]];regionsprof<-listcall[[5]]
  #digits=3;divide=4; memeff = FALSE; fileoutpre="Callobj_"
  copynumber(inputSegmented)<-round(copynumber(inputSegmented),digits)
  segmented(inputSegmented)<-round(segmented(inputSegmented),digits)
  if (divide > nctot) divide <- nctot
  nperturn <- floor(nctot/divide)
  
  for(part in 1:divide){
    print(part)
    if (part < divide) whprof <- ((part-1)*nperturn+1):(part*nperturn) else whprof <- ((part-1)*nperturn+1):nctot
    IS <- inputSegmented[,whprof]
    nc <- length(whprof)
    dataprob<-array(0,c(nclone,nclass*nc))
    print(gc())
     for (k in 1:nc) {
        post        <- (posteriorfin2[posteriorfin2[,1]==whprof[k],])[,-1]
        regionsk    <- (regionsprof[regionsprof[,1]==whprof[k],])[,-1]
        nregk       <- nrow(post)
        probs       <- c()
        for (i in (1:nregk)) {
            regl    <- regionsk[i,2]-regionsk[i,1]+1
            togeth  <- post[i,(1:nclass)]
            probs   <- c(probs,rep(togeth,regl))
        }
        probs <- round(probs,digits)
        allprobs    <- matrix(probs, ncol=nclass, byrow=TRUE)
      #  datk        <- datareg[[1]][,k]
        dataprob[,(nclass*(k-1)+1):(nclass*k)] <- allprobs
        rm(probs,post,regionsk,togeth);print(gc())
    }
    print(gc())
   
    #ncolscl     <- ncol(normalizedData) #25/11/09 moved upward for efficiency reasons
    classify.res        <- array(NA,c(nclone,nc))
    for (i in 1:nc) {
   #     genomdat        <- dataprob[,(ncpp*(i-1)+1)] #removed 16/7/10
        prob.loss.ind   <- dataprob[,(nclass*(i-1)+1)]
        prob.none.ind   <- dataprob[,(nclass*(i-1)+2)]
        prob.gain.ind   <- dataprob[,(nclass*(i-1)+3)]
        if (nclass==4) {
            prob.amp.ind    <- dataprob[,(nclass*(i-1)+4)]
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
    print(gc())
    ncolprob <- ncol(dataprob)
    neworder<-as.vector(sapply(1:nclass,function(x)seq(x,ncolprob,by=nclass)))
    dataprob   <- dataprob[,neworder]
    calls       <- CGHcall:::.assignNames(classify.res, IS)
    rm(classify.res);print(gc());
    probloss <- CGHcall:::.assignNames(dataprob[,1:nc,drop=FALSE], IS)
    dataprob <- dataprob[,-(1:nc),drop=FALSE]
    print(gc())
    probnorm <- CGHcall:::.assignNames(dataprob[,1:nc,drop=FALSE], IS)
    dataprob <- dataprob[,-(1:nc),drop=FALSE]
    print(gc())
    probgain <- CGHcall:::.assignNames(dataprob[,1:nc,drop=FALSE], IS)
    dataprob <- dataprob[,-(1:nc),drop=FALSE]
    print(gc())
    if(nclass==4) probamp <- CGHcall:::.assignNames(dataprob[,1:nc,drop=FALSE], IS)
    rm(dataprob);print(gc())
    if (nclass == 4) {assayData <-assayDataNew(copynumber=copynumber(IS),segmented=segmented(IS),calls=calls,probloss=probloss,probnorm=probnorm,probgain=probgain,probamp=probamp)} else {
    assayData <-assayDataNew(copynumber=copynumber(IS),segmented=segmented(IS),calls=calls,probloss=probloss,probnorm=probnorm,probgain=probgain) }
    rm(probloss,probnorm,probgain); if(nclass==4) rm(probamp);
    print(gc())
   
    result0  <- CGHcall:::.callFromSeg(IS, assayData)  
    if(!memeff){  
        if (part==1) {result <- result0} else {
                result <- combine(result,result0)  
                }
        } else {
        fileout <- paste(fileoutpre,"profiles",whprof[1],"to",whprof[length(whprof)],".Rdata",sep="")  
        save(result0,file=fileout)
    }
    rm(assayData,result0);print(gc())
    }
    
    cat("FINISHED!\n")
    timeFinished <- round((proc.time() - timeStarted)[1] / 60)
    cat("Total time:", timeFinished, "minutes\n")
    if (memeff) print("Results printed to separate files") else return(result)
}
