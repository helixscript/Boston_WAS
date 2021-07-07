getEtaLambda <- function(x,om=NULL){
  if (is.null(om)) om  <- diag(nrow(x$eta))
  rho <- as.vector(t(x$eta[,1:x$etaM])%*%rowSums(om))
  lam <- x$lambda[1:x$lambdaM]
  pobs <- ppois(0, rho%o%lam, lower.tail = FALSE) 
  etalam <- table(x$dataToEta,x$dataToLambda)/pobs
  etacells <- etalam %*% lam
  list(eta=cbind(t(x$eta[,1:x$etaM,drop=FALSE]),x$etaN[1:x$etaM],cells=etacells),
       lambda=cbind(lam,x$lambdaN[1:x$lambdaM],
		    lambdacells=colSums(etalam)*lam))
}

gini <- function(lamobj){
  lamord <- order(lamobj[,1])
  lamobj[] <- lamobj[lamord,]
  prclone <- prop.table(lamobj[,2])
  prcells <- prop.table(lamobj[,3])
  wd <- prclone+c(prclone[-1L],0)
  area <-
    sum(wd*cumsum(prcells))
  area
}

minqgini <- function(x,k=1:4,om=NULL){
  el <- getEtaLambda(x,om)
  y <- el$eta  
  c(
    unlist(lapply(k,function(kelt){
      mins <- combn(ncol(y)-2,kelt,function(xi) apply(y[,xi,drop=FALSE],1,min))
      minkelt <- cbind(mins,y[,-1:0+ncol(y),drop=FALSE])
      nc <- ncol(minkelt)
      ncells <- prop.table(minkelt[,nc])
      ncells %*%  minkelt[,1:(nc-2)] * kelt
    } )),
    gini(el$lambda))
}

minqkeep <- function(pass,log.posterior,i,nkeep){
  om <- eval.parent(quote(om))
  if (i<nkeep)
    c(minqgini(pass,1:4,om),logLik=log.posterior)
  else
    keep.default(pass,log.posterior,i,nkeep)
}

library(parallel)
library(cellTypeCompositions)
load("metadata.RData")
load("gscans.RData")
gseqs <- mclapply(gscans,update,keep=minqkeep,nkeep=201L,mc.cores=length(gscans))
save(gseqs,file="gseqs.RData")
if (!interactive()) q("no")
