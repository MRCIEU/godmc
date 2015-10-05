############################################
# WBC INFERENCE VERSION 1.12
# January 8, 2013
#
# Author: E. Andres Houseman
#
############################################
#
#  Cite:  Houseman EA, Accomando WP, et al.,
#   BMC Bioinformatics (2012)
#
############################################

############################################
#
# Permission to use, copy, modify, and distribute this software and its
# documentation for any purpose other than its incorporation into a
# commercial product is hereby granted without fee, provided that the
# above copyright notice appear in all copies and that both that
# copyright notice and this permission notice appear in supporting
# documentation, and that the name of Brown University not be used in
# advertising or publicity pertaining to distribution of the software
# without specific, written prior permission.

# BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
# PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
# ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#
############################################


options(stringsAsFactors=FALSE)
library(parallel)
library(limma)
library(nlme)



ilogit <- function(M)
{
  res = exp(M) / (1 + exp(M))
  rownames(res) = rownames(M)
  colnames(res) = colnames(M)
  res
}

logit <- function(b)
{
  res = log(b) - log(1 - b)
  rownames(res) = rownames(b)
  colnames(res) = colnames(b)
  res
}

penFitOne <- function(y, Zmat)
{
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

penFitAll <- function(Ymat, Zmat)
{
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


rntransform <- function(x)
{
  out <- rank(x) - 0.5
  out[is.na(x)] <- NA
  mP <- 0.5/max(out, na.rm = T)
  out <- out/(max(out, na.rm = T) + 0.5)
  out <- scale(qnorm(out))
  out
}

inverse.rank.transform <- function(B, mc.cores=mc.cores)
{
  tmpList = lapply(1:mc.cores, function(i){ seq(from=i, to=nrow(B), by=mc.cores) })

  message("Inverse rank transforming data, may take a few minutes...")
  tmpAdj = mclapply(tmpList, function(ix){ apply(B[ix,], 1, rntransform) }, mc.cores=mc.cores)

  message("Reducing results...")
  adjBeta = matrix(NA, nrow(B), ncol(B))
  for (i in 1:length(tmpList)){
      adjBeta[tmpList[[i]],] = t(tmpAdj[[i]])
  }
  rownames(adjBeta) <- rownames(B)
  colnames(adjBeta) <- colnames(B)
  return(adjBeta)
}



estimate.cellcounts <- function(B, top_n=500, cell.coefs)
{
  if (! all(B > 0, na.rm=TRUE)){
      message("performing inverse logit to get values from 0 to 1")
      B = ilogit(B)
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
  return(omega.mix)
}



adjust.beta <- function(B, omega.mix, mc.cores=24)
{
  if (! all(B > 0, na.rm=TRUE)){
      message("performing inverse logit to get values from 0 to 1")
      B = ilogit(B)
      do.logit = TRUE
  } else { do.logit = FALSE }
  stopifnot(all((B > 0) & (B < 1), na.rm=TRUE))

  # after adjusting, set values < 0 or > 1 to the smallest observed value
  epsilon.min = min(B)
  epsilon.max = min(1 - B)

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

  if(do.logit){ 
    adjBeta <- logit(adjBeta)
  }
  return(adjBeta)
}


projectWBC <- function(Y, coefWBC, contrastWBC=NULL, nonnegative=TRUE, lessThanOne=FALSE){ 

  if(is.null(contrastWBC)) Xmat = coefWBC
  else Xmat = coefWBC %*% t(contrastWBC) 

  nCol = dim(Xmat)[2]
  nSubj = dim(Y)[2]

  mixCoef = matrix(0, nSubj, nCol)
  rownames(mixCoef) = colnames(Y)
  colnames(mixCoef) = colnames(Xmat)

  if(nonnegative){
    library(quadprog)

    if(lessThanOne){
      Amat = cbind(rep(-1,nCol), diag(nCol))
      b0vec = c(-1,rep(0,nCol))
    }
    else{
      Amat = diag(nCol)
      b0vec = rep(0,nCol)
    }

    for(i in 1:nSubj){
      obs = which(!is.na(Y[,i])) 
      Dmat = t(Xmat[obs,])%*%Xmat[obs,]
      mixCoef[i,] = solve.QP(Dmat, t(Xmat[obs,])%*%Y[obs,i], Amat, b0vec)$sol
    }
  }
  else{
    for(i in 1:nSubj){
      obs = which(!is.na(Y[,i])) 
      Dmat = t(Xmat[obs,])%*%Xmat[obs,]
      mixCoef[i,] = solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
    }
  }

  return(mixCoef)
}

############################################
############################################

validationWBC <- function(Y, pheno, modelFix, modelBatch=NULL, L.forFstat = NULL){
  N = dim(pheno)[1]
  pheno$y = rep(0, N)
  xTest <- model.matrix(modelFix, pheno)
  sizeModel = dim(xTest)[2]
  
  M = dim(Y)[1]

  if(is.null(L.forFstat)){
     L.forFstat = diag(sizeModel)[-1,]  #All non-intercept coefficients
     colnames(L.forFstat) = colnames(xTest) 
     rownames(L.forFstat) = colnames(xTest)[-1] 
  }

  # Initialize various containers
  sigmaResid = sigmaIcept = nObserved = nClusters = Fstat = rep(NA, M)
  coefEsts = matrix(NA, M, sizeModel)
  coefVcovs =list()

  for(j in 1:M){ # For each CpG

    #Remove missing methylation values
    ii = !is.na(Y[j,])
    nObserved[j] = sum(ii)
    pheno$y = Y[j,]

    if(j%%round(M/10)==0) cat(j,"\n") # Report progress

    try({ # Try to fit a mixed model to adjust for plate
       if(!is.null(modelBatch)){
         fit = try(lme(modelFix, random=modelBatch, data=pheno[ii,]))
         OLS = inherits(fit,"try-error") # If LME can't be fit, just use OLS
       }
       else OLS = TRUE

       if(OLS){
          fit = lm(modelFix, data=pheno[ii,])
          fitCoef = fit$coef
          sigmaResid[j] = summary(fit)$sigma
          sigmaIcept[j] = 0
          nClusters[j] = 0
        }
        else{ 
          fitCoef = fit$coef$fixed
          sigmaResid[j] = fit$sigma
          sigmaIcept[j] = sqrt(getVarCov(fit)[1])
          nClusters[j] = length(fit$coef$random[[1]])
        }
        coefEsts[j,] = fitCoef
        coefVcovs[[j]] = vcov(fit)

        useCoef = L.forFstat %*% fitCoef
        useV = L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
        Fstat[j] = (t(useCoef) %*% solve(useV, useCoef))/sizeModel
    })
  }
  # Name the rows so that they can be easily matched to the target data set
  rownames(coefEsts) = rownames(Y)
  colnames(coefEsts) = names(fitCoef)
  
  degFree = nObserved - nClusters - sizeModel + 1

  # Get P values corresponding to F statistics
  Pval = 1-pf(Fstat, sizeModel, degFree)

  out = list(coefEsts=coefEsts, coefVcovs=coefVcovs, modelFix=modelFix, modelBatch=modelBatch,
     sigmaIcept=sigmaIcept, sigmaResid=sigmaResid, L.forFstat=L.forFstat, Pval=Pval, orderFstat=order(-Fstat),
     Fstat=Fstat, nClusters=nClusters, nObserved=nObserved, degFree=degFree)

  class(out) = "validationWBC"
  out
}

############################################

print.validationWBC <- function(x, digits=4){
  cat("WBC validation object\n\n",sep="")
  cat(dim(x$coefEsts)[1],"CpGs assessed\n\n")
  cat("Fixed effects model:\n")
  print(x$modelFix)
  if(is.null(x$modelBatch)){
    cat("\n(No batch adjustment, lm used)")
  }
  else{
    cat("\nRandom effects model used in lme for batch:\n")
    print(x$modelBatch)
  }
  cat("\nAverage sample size: ", round(mean(x$nObserved),2))
  cat("\nAverage number of batches: ", round(mean(x$nClusters),2))
  cat("\nAverage degrees-of-freedom: ", round(mean(x$degFree),2))
  cat("\n")
}

############################################
############################################

inferWBCbyLme <- function(Y, pheno, modelFix, modelRand, coefWBC, 
   contrastWBC=NULL, detailLevel=1, silentErrors=TRUE){

  M = dim(coefWBC)[1]
  n = dim(Y)[2]
  if(dim(Y)[1] != M) stop("Y does not match coefWBC in dimension\n")

  if(detailLevel == -1){
    lmeVcomp = list(
      mu = matrix(NA,M,n),
      e = matrix(NA,M,n),
      a = list()
     )
  }

  sigmaResid = sigmaIcept = nObserved = nClusters = rep(NA, M)
  coefEsts = list()
  for(j in 1:M){
    ii = !is.na(Y[j,])
    nObserved[j] = sum(ii)
    pheno$y = Y[j,]

    try({
      fit = try(lme(modelFix, random=modelRand, data=pheno[ii,]), 
        silent=silentErrors)

      if(inherits(fit,"try-error")){
        fit = lm(modelFix, data=pheno[ii,])
        fitCoef = fit$coef
        sigmaResid[j] = summary(fit)$sigma
        sigmaIcept[j] = 0
        nClusters[j] = 0
        if(detailLevel == -1){
           lmeVcomp$mu[j,ii] = predict(fit)
           lmeVcomp$e[j,ii] = residuals(fit)
           lmeVcomp$a[[j]] = list()
        }
      }
      else{
        fitCoef = fit$coef$fixed
        sigmaResid[j] = fit$sigma
        sigmaIcept[j] = sqrt(getVarCov(fit)[1])
        nClusters[j] = length(fit$coef$random[[1]])
        if(detailLevel == -1){
           lmeVcomp$mu[j,ii] = predict(fit,level=0)
           lmeVcomp$e[j,ii] = residuals(fit,level=1)
           lmeVcomp$a[[j]] = fit$coef$random
        }
      }
      coefEsts[[j]] = fitCoef
    })
  }
  coefEsts = t(matrix(unlist(coefEsts), length(fitCoef), M))
  rownames(coefEsts) = rownames(coefWBC)
  colnames(coefEsts) = names(fitCoef)

  if(is.null(contrastWBC)) Z = cbind(1,coefWBC)
  else Z = cbind(1, coefWBC %*% t(contrastWBC))
  colnames(Z)[1] = "<Intercept>"

  ZtZ = t(Z) %*% Z
  ZtZqr = qr(ZtZ)
  G = solve(ZtZqr, t(Z) %*% coefEsts)

  out = list(method="lme", GammaMatrix=G, Beta1Matrix=coefEsts, 
    sigma=cbind(resid=sigmaResid, intercept=sigmaIcept),
    N=cbind(obs=nObserved,clusters=nClusters))

  pheno$y[] <- 0 # EAH 12-08-05 (Prevents NAs in the last CpG from killing the next loop
  if(detailLevel>=1){
    out$var.unscaled=solve(ZtZqr)

    d = dim(Z)[2]
    out$predicted = Z %*% G
    residual2 = coefEsts - out$predicted
    out$Sigma = (1/(M-d)) * t(residual2)%*%residual2

    X = model.matrix(modelFix, pheno)
    Xbar = apply(X,2,mean)
    Xctr = (X - matrix(1,n,1) %*% Xbar)
    residual2z = Xctr %*% t(residual2)
    out$ssu = apply(residual2z*residual2z, 2, sum)

    residual1 = Y - coefEsts %*% t(X) 
    out$sse = apply(residual1*residual1, 1, sum, na.rm=TRUE)

    residual0 = Y - apply(Y,1,mean,na.rm=TRUE) %*% matrix(1,1,n)
    out$ss0 = apply(residual0*residual0, 1, sum, na.rm=TRUE)

    if(detailLevel>=1){
       out$residuals = list(stage1=residual1, stage2=residual2, inner=t(residual2z))
    }
  }
  if(detailLevel==-1){
    trueLme = which(sapply(lmeVcomp$a, length)>0)
    ch = sort(unique(unlist(lapply(lmeVcomp$a[trueLme], function(u)rownames(u[[1]])))))
    nch = length(ch)
    a = matrix(0, M, nch)
    colnames(a) = ch
    for(j in trueLme){
      aa = lmeVcomp$a[[j]][[1]]
      a[j,rownames(aa)] = aa
    }
    lmeVcomp$ch = names(lmeVcomp$a[[trueLme[1]]])
    lmeVcomp$a = a
    out$lmeVcomp = lmeVcomp
  }

  class(out) = "inferWBC"
  out
}

############################################

summary.inferWBC <- function(x, xboot=NULL){

  out = list(Coef=x$Gamma)
  class(out) = "inferWBCsummary"

  if(is.null(x$var.unscaled)){
     return(out)
  }

  sUnscaled = sqrt(diag(x$var.unscaled))
  sigma = sqrt(diag(x$Sigma))
  se = outer(sUnscaled,sigma,"*")
  rownames(se) = rownames(x$Gamma)
  out$StdErrNaive=se

  class(out) = "inferWBCsummary"

  if(!is.null(x$ss0)){
    out$RsquareTotal = 1-sum(x$ssu+x$sse)/sum(x$ss0)
    out$RsquareStage1 = 1-sum(x$ssu)/sum(x$ss0-x$sse)
  }

  if(is.null(xboot)) return(out)

  sboot = summary(xboot)

  out$Bias.Single = x$Gamma - sboot$Mean.Single
  out$Bias.Double = x$Gamma - sboot$Mean.Double 
  out$SD.Single = sboot$SD.Single
  out$SD.Double = sboot$SD.Double
  out$numBoots = sboot$numBoots

  out
}

############################################

print.inferWBC <- function(x, digits=4){
  cat("WBC inference object (method = ",x$method,")\n\n",sep="")
  print(summary(x), digits=digits)
}

############################################

print.inferWBCsummary <- function(x, digits=4, displayFactor=100){
  vars = setdiff(colnames(x$Coef),"(Intercept)")
  if(length(x)==1) {
    cat("(no inference data)\n")
    print(x$Coef[,vars])
    return()
  }

  tag = flag = "?"
  tab = NULL
  for(v in vars){
    cat(v,"\n")
    tab = cbind(Est=x$Coef[,v], StdErr0=x$StdErrNaive[,v])
    if(is.null(x$SD.Single)){
       flag = "StdErr0"
       tag = "naive"
    }
    else{
       tab = cbind(tab, StdErr1=x$SD.Single[,v])
       if(is.null(x$SD.Double)){
          flag = "StdErr1"
          tag = "single bootstrap"
       }
       else{
         tab = cbind(tab, StdErr2=x$SD.Double[,v])
          flag = "StdErr2"
          tag = "double bootstrap"
       }
    }
    tab = tab*displayFactor
    tab = cbind(tab, Zscore = as.vector(x$Coef[,v]/tab[,flag])*displayFactor)
    tab = cbind(tab, Pvalue = 2*pnorm(-abs(tab[,"Zscore"])))
    print(tab, digits=digits)
    cat("\n")
  }

  cat("Inference based on ", tag, " standard errors (", flag, ")\n",sep="")
  if(dim(tab)[2]>4) cat("\tfrom",x$numBoots,"bootstrap iterations\n\n")
  else cat("\n")

  cat("Proportion of total variation explained by WBC:", 
       format(x$RsquareTotal,digits=digits),"\n")
  cat("Proportion of stage 1 model explained by WBC:", 
       format(x$RsquareStage1,digits=digits),"\n")

}

############################################

bootInferWBCbyLme <- function(Y, pheno, modelFix, modelRand, 
  coefWBC, contrastWBC=NULL, strata=NULL, R=10,
  vcovWBC=NULL, # Supply this as a list to compute double-bootstrap
  degFree=NULL,  # Supply this as a vector to use t-distribution for 2-boot
  saveBetaCoefficients = FALSE,
  silentErrors = TRUE
){

  M = dim(coefWBC)[1]
  d = dim(coefWBC)[2]

  if(is.null(contrastWBC)) {
     dBeta0 = d
     contrastWBC = diag(d)
     rownames(contrastWBC) = colnames(coefWBC)
  }
  else {
     dBeta0 = dim(contrastWBC)[1]
  }

  n = dim(pheno)[1]
  if(is.null(strata)) stratList = list(1:n)
  else{
    stratList = split(1:n, strata)
  }
  nStrata = length(stratList)
  stratListN = sapply(stratList, length)

  doDoubleBoot = FALSE
  if(!is.null(vcovWBC)){
     doDoubleBoot = TRUE
     oneR = matrix(1,R,1)
     Beta0 = array(NA, dim=c(M, dBeta0, R))
     for(j in 1:M){
       vchol = chol(vcovWBC[[j]])
       if(is.null(degFree)) noise = rnorm(R*d)
       else noise = rt(R*d, degFree[j])

       bootWBC = t(oneR %*% coefWBC[j,] + matrix(noise, R, d) %*% vchol)
       Beta0[j,,] = contrastWBC %*% bootWBC
     }
  }

  fit0 = inferWBCbyLme(Y, pheno, modelFix, modelRand, 
   coefWBC, contrastWBC, detailLevel=-1, silentErrors=silentErrors)
  nchip = dim(fit0$lmeVcomp$a)[2]
  chips = as.character(pheno[[fit0$lmeVcomp$ch]])

  GammaList = Beta1List = list()
  for(r in 1:R){
    bootsL = list()
    for(s in 1:nStrata){
       bootsL[[s]] = sample(stratList[[s]], stratListN[s], replace=TRUE)
    }
    boots = unlist(bootsL)

    Ystar = fit0$lmeVcomp$mu + fit0$lmeVcomp$e[,boots] 
    astar = fit0$lmeVcomp$a[,sample(1:nchip, nchip, replace=TRUE)]

    colnames(astar) = colnames(fit0$lmeVcomp$a)

    Ystar = Ystar + astar[,chips]

    GammaObject = inferWBCbyLme(Ystar, pheno, modelFix, modelRand, 
       coefWBC, contrastWBC, detailLevel=0, silentErrors=silentErrors)

    if(r %% 10==0) cat(r,"\n")
    GammaList[[r]] = GammaObject$GammaMatrix
    Beta1List[[r]] = GammaObject$Beta1Matrix
  }

  Z1 = cbind(1, coefWBC %*% t(contrastWBC))

  out = list()
  #out = list(numBoots=R, tDistribution=!is.null(degFree))

  out$Gamma1 = array(NA, dim=c(dBeta0+1, dim(Beta1List[[1]])[2], R))
  for(r in 1:R) {
    out$Gamma1[,,r] = solve(t(Z1)%*%Z1, t(Z1) %*% Beta1List[[r]])
  }

  dimnames(out$Gamma1) = c(dimnames(GammaList[[1]]), list(1:R))

  if(doDoubleBoot){
    out$Gamma2 = array(NA, dim=c(dBeta0+1, dim(Beta1List[[1]])[2], R))
    for(r in 1:R) {
      Z2 = cbind(1, Beta0[,,r])
      out$Gamma2[,,r] = solve(t(Z2)%*%Z2, t(Z2) %*% Beta1List[[r]])
    }
    dimnames(out$Gamma2) = c(dimnames(GammaList[[1]]), list(1:R))
  }

  if(saveBetaCoefficients) {
    if(doDoubleBoot) out$Beta0 = Beta0
    out$Beta1 = array(NA, dim=c(dim(Beta1List[[1]]),R))
    for(r in 1:R) {
      out$Beta1[,,r] = Beta1List[[r]]
    }
    dimnames(out$Beta1) = c(dimnames(Beta1List[[1]]), list(1:R))
  }

  class(out) = "inferWBCBoot"
  out
}

############################################

summary.inferWBCBoot <- function(x){
  out = list(numBoots = dim(x$Gamma1)[3])
  out$Mean.Single = apply(x$Gamma1,1:2,mean)
  out$SD.Single = apply(x$Gamma1,1:2,sd)
  if(!is.null(x$Gamma2)){
    out$Mean.Double = apply(x$Gamma2,1:2,mean)
    out$SD.Double = apply(x$Gamma2,1:2,sd)
  }
  out
}

############################################

print.inferWBCBoot <- function(x, digits=NULL){
  print(summary(x), digits=digits)
}

############################################
############################################

inferWBCbyLm <- function(Y, pheno, modelFix, coefWBC, 
   contrastWBC=NULL, detailLevel=1, silentErrors=TRUE){

  M = dim(coefWBC)[1]
  n = dim(Y)[2]
  if(dim(Y)[1] != M) stop("Y does not match coefWBC in dimension\n")

  if(detailLevel == -1){
    lmVcomp = list(
      mu = matrix(NA,M,n),
      e = matrix(NA,M,n)
     )
  }

  sigmaResid = sigmaIcept = nObserved = nClusters = rep(NA, M)
  coefEsts = list()
  for(j in 1:M){
    ii = !is.na(Y[j,])
    nObserved[j] = sum(ii)
    pheno$y = Y[j,]

     fit = lm(modelFix, data=pheno[ii,])
     fitCoef = fit$coef
     sigmaResid[j] = summary(fit)$sigma
     sigmaIcept[j] = 0
     nClusters[j] = 0
     if(detailLevel == -1){
        lmVcomp$mu[j,ii] = predict(fit)
        lmVcomp$e[j,ii] = residuals(fit)
     }
     coefEsts[[j]] = fitCoef
  }
  coefEsts = t(matrix(unlist(coefEsts), length(fitCoef), M))
  rownames(coefEsts) = rownames(coefWBC)
  colnames(coefEsts) = names(fitCoef)

  if(is.null(contrastWBC)) Z = cbind(1,coefWBC)
  else Z = cbind(1, coefWBC %*% t(contrastWBC))
  colnames(Z)[1] = "<Intercept>"

  ZtZ = t(Z) %*% Z
  ZtZqr = qr(ZtZ)
  G = solve(ZtZqr, t(Z) %*% coefEsts)

  out = list(method="lm", GammaMatrix=G, Beta1Matrix=coefEsts, 
    sigma=sigmaResid,N=nObserved)

  pheno$y[] <- 0 # EAH 12-08-05 (Prevents NAs in the last CpG from killing the next loop
  if(detailLevel>=1){
    out$var.unscaled=solve(ZtZqr)

    d = dim(Z)[2]
    out$predicted = Z %*% G
    residual2 = coefEsts - out$predicted
    out$Sigma = (1/(M-d)) * t(residual2)%*%residual2

    X = model.matrix(modelFix, pheno)
    Xbar = apply(X,2,mean)
    Xctr = (X - matrix(1,n,1) %*% Xbar)
    residual2z = Xctr %*% t(residual2)
    out$ssu = apply(residual2z*residual2z, 2, sum)

    residual1 = Y - coefEsts %*% t(X) 
    out$sse = apply(residual1*residual1, 1, sum, na.rm=TRUE)

    residual0 = Y - apply(Y,1,mean,na.rm=TRUE) %*% matrix(1,1,n)
    out$ss0 = apply(residual0*residual0, 1, sum, na.rm=TRUE)

    if(detailLevel>=1){
       out$residuals = list(stage1=residual1, stage2=residual2, inner=t(residual2z))
    }
  }
  if(detailLevel==-1){
    out$lmVcomp = lmVcomp
  }

  class(out) = "inferWBC"
  out
}

############################################

bootInferWBCbyLm <- function(Y, pheno, modelFix, 
  coefWBC, contrastWBC=NULL, strata=NULL, R=10,
  vcovWBC=NULL, # Supply this as a list to compute double-bootstrap
  degFree=NULL,  # Supply this as a vector to use t-distribution for 2-boot
  saveBetaCoefficients = FALSE,
  silentErrors = TRUE
){

  M = dim(coefWBC)[1]
  d = dim(coefWBC)[2]

  if(is.null(contrastWBC)) {
     dBeta0 = d
     contrastWBC = diag(d)
     rownames(contrastWBC) = colnames(coefWBC)
  }
  else {
     dBeta0 = dim(contrastWBC)[1]
  }

  n = dim(pheno)[1]
  if(is.null(strata)) stratList = list(1:n)
  else{
    stratList = split(1:n, strata)
  }
  nStrata = length(stratList)
  stratListN = sapply(stratList, length)

  doDoubleBoot = FALSE
  if(!is.null(vcovWBC)){
     doDoubleBoot = TRUE
     oneR = matrix(1,R,1)
     Beta0 = array(NA, dim=c(M, dBeta0, R))
     for(j in 1:M){
       vchol = chol(vcovWBC[[j]])
       if(is.null(degFree)) noise = rnorm(R*d)
       else noise = rt(R*d, degFree[j])

       bootWBC = t(oneR %*% coefWBC[j,] + matrix(noise, R, d) %*% vchol)
       Beta0[j,,] = contrastWBC %*% bootWBC
     }
  }

  fit0 = inferWBCbyLm(Y, pheno, modelFix,
   coefWBC, contrastWBC, detailLevel=-1, silentErrors=silentErrors)

  GammaList = Beta1List = list()
  for(r in 1:R){
    bootsL = list()
    for(s in 1:nStrata){
       bootsL[[s]] = sample(stratList[[s]], stratListN[s], replace=TRUE)
    }
    boots = unlist(bootsL)
    Ystar = fit0$lmVcomp$mu + fit0$lmVcomp$e[,boots] 

    GammaObject = inferWBCbyLm(Ystar, pheno, modelFix,
       coefWBC, contrastWBC, detailLevel=0, silentErrors=silentErrors)

    if(r %% 10==0) cat(r,"\n")
    GammaList[[r]] = GammaObject$GammaMatrix
    Beta1List[[r]] = GammaObject$Beta1Matrix
  }

  Z1 = cbind(1, coefWBC %*% t(contrastWBC))

  out = list()

  out$Gamma1 = array(NA, dim=c(dBeta0+1, dim(Beta1List[[1]])[2], R))
  for(r in 1:R) {
    out$Gamma1[,,r] = solve(t(Z1)%*%Z1, t(Z1) %*% Beta1List[[r]])
  }

  dimnames(out$Gamma1) = c(dimnames(GammaList[[1]]), list(1:R))

  if(doDoubleBoot){
    out$Gamma2 = array(NA, dim=c(dBeta0+1, dim(Beta1List[[1]])[2], R))
    for(r in 1:R) {
      Z2 = cbind(1, Beta0[,,r])
      out$Gamma2[,,r] = solve(t(Z2)%*%Z2, t(Z2) %*% Beta1List[[r]])
    }
    dimnames(out$Gamma2) = c(dimnames(GammaList[[1]]), list(1:R))
  }

  if(saveBetaCoefficients) {
    if(doDoubleBoot) out$Beta0 = Beta0
    out$Beta1 = array(NA, dim=c(dim(Beta1List[[1]]),R))
    for(r in 1:R) {
      out$Beta1[,,r] = Beta1List[[r]]
    }
    dimnames(out$Beta1) = c(dimnames(Beta1List[[1]]), list(1:R))
  }

  class(out) = "inferWBCBoot"
  out
}

############################################
############################################

initializeBootSampleValidation <- function(validationObj, selection, Lwbc = NULL){
  if(is.null(Lwbc)) Lwbc = diag(dim(validationObj$coefEsts)[2])
  coefEstSel = validationObj$coefEsts[selection,] %*% t(Lwbc)

  coefVcovSel <- lapply(validationObj$coefVcovs[selection], function(v) Lwbc %*% v %*% t(Lwbc))
  LcholSel <- lapply(coefVcovSel, function(v) chol(v))
  degFreeSel <- validationObj$degFree[selection]

  list(coefEst = coefEstSel, Lchol=LcholSel, degfree=degFreeSel, d=dim(coefEstSel)[2], J0=dim(coefEstSel)[1])
}


bootSampleValidation <- function(inits, tDistribution=TRUE){
  if(tDistribution) rError = function(dendf){rt(inits$d,dendf)}
  else rError=function(dendf)rnorm(inits$d)
  Ysamp <- inits$coefEst
  for(j in 1:inits$J0){
    Ysamp[j,] <- Ysamp[j,] + t(inits$Lchol[[j]]) %*% rError(inits$degfree[j])
  }
  Ysamp
}

