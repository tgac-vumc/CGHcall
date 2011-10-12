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

.getCentromere <- function() {
    ### Centromere data from http://genome.ucsc.edu
    ### Database March 2006
    ### Date: 11 september 2006
    
    centromere       <- matrix(NA, 23, 2);
    centromere[1,1]  <- 121236957;    
    centromere[1,2]  <- 123476957;
    centromere[2,1]  <- 91689898;
    centromere[2,2]  <- 94689898;
    centromere[3,1]  <- 90587544;
    centromere[3,2]  <- 93487544;
    centromere[4,1]  <- 49354874;
    centromere[4,2]  <- 52354874;
    centromere[5,1]  <- 46441398; 
    centromere[5,2]  <- 49441398;
    centromere[6,1]  <- 58938125; 
    centromere[6,2]  <- 61938125;
    centromere[7,1]  <- 58058273; 
    centromere[7,2]  <- 61058273;
    centromere[8,1]  <- 43958052; 
    centromere[8,2]  <- 46958052;
    centromere[9,1]  <- 47107499; 
    centromere[9,2]  <- 50107499;
    centromere[10,1] <- 39244941; 
    centromere[10,2] <- 41624941;
    centromere[11,1] <- 51450781; 
    centromere[11,2] <- 54450781;
    centromere[12,1] <- 34747961; 
    centromere[12,2] <- 36142961;
    centromere[13,1] <- 16000000; 
    centromere[13,2] <- 17868000;
    centromere[14,1] <- 15070000; 
    centromere[14,2] <- 18070000;
    centromere[15,1] <- 15260000; 
    centromere[15,2] <- 18260000;
    centromere[16,1] <- 35143302; 
    centromere[16,2] <- 36943302;
    centromere[17,1] <- 22187133;
    centromere[17,2] <- 22287133;
    centromere[18,1] <- 15400898; 
    centromere[18,2] <- 16764896;
    centromere[19,1] <- 26923622; 
    centromere[19,2] <- 29923622;
    centromere[20,1] <- 26267569; 
    centromere[20,2] <- 28033230;
    centromere[21,1] <- 10260000; 
    centromere[21,2] <- 13260000;
    centromere[22,1] <- 11330000; 
    centromere[22,2] <- 14330000;
    centromere[23,1] <- 58598737; 
    centromere[23,2] <- 61598737;
    centromere       <- apply(centromere, 1, mean);
    return(centromere);
}

#.convertChromosomeToArm <- function(dataframe) {
#    cat("Dividing chromosomes into arms:\n\n");
#    centromere  <- CGHcall:::.getCentromere();
#    arm         <- 1;
#    output      <- c();
#    previous    <- 1;
#    a           <- 0;
#    prevpos     <- dataframe[1,3];
#    cat("New chromosome:\t\t", previous, "\t\tArm:\t", arm, "\n");
#    temp <- dataframe[,2:3];
#    for (i in 1:length(temp[[1]])) {
#        if (temp[i,1] != previous) {
#            arm <- arm + 1;
#            cat("New chromosome:\t\t", temp[i,1], "\t\tArm:\t", arm, "\n");
#        }
#        else {
#            if (prevpos < centromere[temp[i,1]] && temp[i,2] >= centromere[temp[i,1]]) {
#                arm <- arm + 1;            
#                cat("Centromere found:\t", centromere[temp[i,1]], "\tArm:\t", arm, "\n");
#            }
#        }
#        output   <- c(output, arm);
#        previous <- temp[i,1];
#        prevpos  <- temp[i,2];
#    }
#    dataframe[,2] <- output;
#    return(dataframe);    
#}


.convertChromosomeToArm <- function(dataframe) { #changed 22/06/2009; 
    cat("Dividing chromosomes into arms:\n\n");
    centromere  <- CGHcall:::.getCentromere();
    chr <- dataframe[,2]
    bp <- dataframe[,3]
    chrlev <- unique(chr)
    a<-1
    chrarms <- c()
    for(i in chrlev){
    print(i)
    chri <- which(chr==i)
    bpi <- bp[chri]
    wbpi <- length(which(bpi<=centromere[i]))
    wbpil <- length(which(bpi>centromere[i]))
    if(wbpi>0) {chrarms <- c(chrarms,rep(a,wbpi));a<-a+1}
    if(wbpil>0) {chrarms <- c(chrarms,rep(a,wbpil));a<-a+1}  
    }  
    dataframe[,2] <- chrarms
    return(dataframe)
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

.xgivenknotrunc <- function(k, class, pm, varprofall, allsum, allsumsq, allnc, robustsig) {

    sumx    <- allsum[k]
    sumxsq  <- allsumsq[k]
    ncl     <- allnc[k]
    v1      <- 2*varprofall[k]

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
    
    v2 <- if(class == 4 | class == 2) 2*((pm[5+class])^2 + 0.0001) 
    else {if(class==3) {if(robustsig=="yes") 1/4*((pm[7])^2 + 0.0001)+1/4*((pm[9])^2 + 0.0001)+2*((pm[8]^2)+0.0001)
    else {2*((pm[8])^2 + 0.0001)}}
    else {if(class==1) 2*(pm[7]^2+(pm[6])^2 + 0.0001) 
    else {if(class==5) 2*((pm[10])^2+(pm[9])^2 + 0.0001)
    else 2*((pm[11])^2 + (pm[10])^2 + (pm[9])^2 + 0.0001)}}} # for stability


    A       <- (sumx*v2+mu*v1)/(ncl*v2+v1)
    B       <- v1*v2/(ncl*v2+v1)
    C       <- sumxsq/v1 + mu^2/v2 - A^2/B
    logD    <- log(B)/2 -log(v2)/2
    logres  <- - C + logD
    res     <- exp(logres)
    return(logres)
}

.posteriorp <- function(k, priorp, pm, varprofall, allsum, allsumsq, allnc,robustsig) {
   all3        <- sapply(c(1,2,3,4,5,6), CGHcall:::.xgivenknotrunc, k=k, pm=pm, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc,robustsig=robustsig)
    maxim       <- max(all3)
    all         <- exp(all3-maxim)
    allprior    <- all*priorp[k,]
    tot         <- all%*%priorp[k,]
    return(allprior/tot)
}

.reallikk4 <- function(k, alpha, pm, varprofall, allsum, allsumsq, allnc,robustsig) {
    all3    <-sapply(c(1,2,3,4,5,6), CGHcall:::.xgivenknotrunc, k=k, pm=pm, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc,robustsig=robustsig)
    maxim   <- max(all3)
    all     <- exp(all3-maxim)
    tot     <- maxim + log(all%*%alpha[k,])
    return(tot)
}

.reallik4 <- function(nreg, alpha, pm, varprofall, allsum, allsumsq, allnc,robustsig) {
    return(sum(sapply(1:nreg, CGHcall:::.reallikk4, pm=pm, alpha=alpha, varprofall =varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc,robustsig=robustsig)))
}

.alphafun <- function(k, profile, priorp, pm, varprofall, allsum, allsumsq, allnc=allnc, robustsig, prior) {
    prof        <- profile[k]
    regionsk    <- which(profile==prof)
    nregk       <- length(regionsk)
    if (prior=="all") {
        nregk       <- length(profile)
        regionsk    <- 1:nregk
    }
    totpost     <- rep(1,nregk)%*%t(sapply(regionsk, CGHcall:::.posteriorp, priorp=priorp, pm=pm, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc=allnc,robustsig=robustsig))
    return(totpost/nregk)
}

.alpha0all <- function(nreg, profile, priorp, pm, varprofall, allsum, allsumsq, allnc, robustsig, prior) {
    prevprof    <- 0
    alphaall    <- c()
    alpha1      <- CGHcall:::.alphafun(1, profile, priorp, pm, varprofall, allsum, allsumsq, allnc, robustsig, prior)
    for (i in (1:nreg)) {
        if (prior != "all") {
            curprof <- profile[i]
            if (curprof==prevprof) {
                newalpha <- oldalpha
            }
            if (curprof != prevprof) {
                newalpha <- CGHcall:::.alphafun(i, profile, priorp, pm, varprofall, allsum, allsumsq, allnc, robustsig, prior)
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


.minusEloglikreg <- function(k, posteriorprev, alphaprev, pm, varprofall, allsum, allsumsq, allnc, robustsig) {
    all3 <- sapply(c(1,2,3,4,5,6), CGHcall:::.xgivenknotrunc, k=k, pm=pm, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc,robustsig=robustsig)
    return(-all3%*%posteriorprev[,k] - (log(alphaprev[k,])%*%posteriorprev[,k]))
} 
  
.totallik <- function(nreg, pm, posteriorprev, alphaprev, varprofall, allsum, allsumsq, allnc, robustsig) {
    return(sum(sapply(1:nreg, CGHcall:::.minusEloglikreg, pm=pm, posteriorprev=posteriorprev, alphaprev=alphaprev, varprofall=varprofall, allsum=allsum, allsumsq=allsumsq, allnc,robustsig=robustsig)))
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
