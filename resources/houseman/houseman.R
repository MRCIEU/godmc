options(stringsAsFactors=FALSE)
library(parallel)
library(limma)
library(nlme)

#source('~/src/rstuff/rmodels.R')
source('wbcInference-V112.R')


ilogit = function(M){
    res = exp(M) / (1 + exp(M))
    rownames(res) = rownames(M)
    colnames(res) = colnames(M)
    res
}
logit = function(b){
    res = log(b) - log(1 - b)
    rownames(res) = rownames(b)
    colnames(res) = colnames(b)
    res
}

penFitOne = function(y, Zmat){
  adj = y
  is.obs = !is.na(y)
  Z = Zmat[is.obs,,drop=FALSE]
  y = y[is.obs]
  id = rep(1,length(y))
  lmod = try(lme(y~1, random=list(id=pdIdent(~Z-1))), silent=TRUE)
  if(!inherits(lmod,"try-error")){
    adj[is.obs] = resid(lmod) + lmod$coef$fixed[1]
  }# else {
  #    message(lmod)
  #}
  list(modelFit=lmod, adjusted=adj)
}

penFitAll = function(Ymat, Zmat){
  nFeature = dim(Ymat)[1]

  mu = sigma = tau = rep(NA, nFeature)
  beta = matrix(NA, nFeature, dim(Zmat)[2])
  adjusted = matrix(NA, nFeature, dim(Ymat)[2])
  nbad = 0
  for(i in 1:nFeature){
    pf = penFitOne(Ymat[i,], Zmat)

    if(inherits(pf$modelFit,"try-error")){
      nbad = nbad + 1
      mu[i] = mean(Ymat[i,],na.rm=TRUE)
      sigma[i] = var(Ymat[i,],na.rm=TRUE)
      tau[i] = 0
      beta[i,] = 0
    }
    else{
      mu[i] = pf$modelFit$coef$fixed[1]
      sigma[i] = pf$modelFit$sigma^2
      tau[i] = getVarCov(pf$modelFit)[1,1]
      beta[i,] = pf$modelFit$coef$random$id
    }
    adjusted[i,] = pf$adjusted
  }
  message(paste("total with error in model fit:", nbad))
  list(mu=mu, beta=beta, tau=tau, sigma=sigma, adjusted=adjusted)
}

adjust.beta = function(B, top_n=500, mc.cores=24, cell.coefs="houseman-dmrs.txt", est.only=FALSE){
    if (! all(B > 0, na.rm=TRUE)){
        message("performing inverse logit to get values from 0 to 1")
        B = ilogit(B)
        print(dim(B))
        do.logit = TRUE
    } else { do.logit = FALSE }
    stopifnot(all((B > 0) & (B < 1), na.rm=TRUE))

    # after adjusting, set values < 0 or > 1 to the smallest observed value
    epsilon.min = min(B)
    epsilon.max = min(1 - B)
    
    dmr.coefs = read.delim(cell.coefs, row.names=1)
    # take shared probes. may be differences if beta is from 450k
    dmr.coefs = dmr.coefs[rownames(dmr.coefs) %in% rownames(B),]
    dmr.coefs = as.matrix(dmr.coefs[1:min(top_n, nrow(dmr.coefs)),])
    stopifnot(nrow(dmr.coefs) > 0)

    # reporter matrix
    Lwbc = diag(8)[-(1:2),]
    Lwbc[1:2,2] = 1 ; Lwbc[,1] = 1
    rownames(Lwbc) = colnames(dmr.coefs)[-(1:2)]
    colnames(Lwbc) = colnames(dmr.coefs)
    dmrs = rownames(dmr.coefs)

    message(paste0("using ", length(dmrs), " dmrs"))

    omega.mix = projectWBC(B[dmrs,],
        dmr.coefs[dmrs,],
        contrastWBC=Lwbc,
        nonnegative=TRUE,
        lessThanOne=FALSE)


    omega.mix = (1 / apply(omega.mix, 1, sum)) * omega.mix
    if(est.only){ return(omega.mix) }
    #write.table(omega.mix, sep="\t", row.names=T, quote=F, file="tmp/omega.mix.txt")
    #message("WROTE")

    message("adjusting beta (this will take a while)...")
    tmpList = lapply(1:mc.cores, function(i){ seq(from=i, to=nrow(B), by=mc.cores) })
    tmpAdj = mclapply(tmpList, function(ix){ penFitAll(B[ix,], omega.mix) }, mc.cores=mc.cores)

    adjBeta = matrix(NA, nrow(B), ncol(B))
    for (i in 1:length(tmpList)){
        adjBeta[tmpList[[i]],] = tmpAdj[[i]]$adjusted
    }
    nCpGs = nrow(B) * ncol(B)
    message(paste("% values <= 0:", sum(adjBeta <= 0, na.rm=TRUE) / nCpGs * 100))
    message(paste("% values >= 1:", sum(adjBeta >= 1, na.rm=TRUE) / nCpGs * 100))
    message("these will be changed to 0+epsilon, 1-epsilon respectively")

    adjBeta[(adjBeta <= 0)] = epsilon.min
    adjBeta[(adjBeta >= 1)] = epsilon.max
    rownames(adjBeta) = rownames(B)
    colnames(adjBeta) = colnames(B)

    if(do.logit){ return(logit(adjBeta)) };
    return(adjBeta)
}