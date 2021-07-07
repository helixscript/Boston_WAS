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
load("gscans.RData")
load("metadata.RData")

tiltom <- function(omelt,tilt=0.10){
  onecol <- apply(omelt,1,max) * (1+tilt) > 1.0
  sapply(1:4,function(x){
    up <- rep(1,4)
    if (onecol[x])
      up[-x] <- 1/(1+tilt)
    else
      up[x] <- 1+tilt
    down <- ifelse((1:4) == x,1/(1+tilt),1)
    list(up=omelt%*%diag(up),down=omelt%*%diag(down))
  }
  )
}

om3 <- sapply(om2,tiltom,simplify="array")
## 2 by 4 by 18 list of 4 by 4 matrices

gscan.calls <-
  mapply(
	 function(dir,col,scan) {
	   gcall <- attr(gscans[[scan]],"call")
	   gcall[["om"]] <-
	     parse(text=paste("om3[[",
			      dir,",",
			      col,",",
			      scan,"]]"))[[1L]]
	   gcall[["nkeep"]] <- 21L
	   gcall[["keep"]] <- as.symbol("minqkeep")
	   gcall
	 },
    slice.index(om3,1),slice.index(om3,2),slice.index(om3,3),
    SIMPLIFY=FALSE)


gseqvary <- mclapply(gscan.calls, eval, mc.cores=length(gscans))
save(gseqvary,file="gseqvary.RData")
if (!interactive()) q("no")
