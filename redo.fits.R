library(cellTypeCompositions)
library(parallel)
getEtaLambda <- function(x,om=NULL){
  if (is.null(om)) om  <- diag(nrow(x$eta))
  rho <- as.vector(t(x$eta[,1:x$etaM, drop=FALSE])%*%rowSums(om))
  lam <- x$lambda[1:x$lambdaM]
  lamN <- x$lambdaN[1:x$lambdaM]
  pobs <- ppois(0, rho%o%lam, lower.tail = FALSE) 
  etalam <- unname(table(x$dataToEta,x$dataToLambda))/pobs
  etacells <- etalam %*% lam
  etaclones <- rowSums(etalam)
  mat <- cbind(X=t(x$eta[,1:x$etaM,drop=FALSE]),
	       x$etaN[1:x$etaM],
	       etaclones,
	       etacells)
  colnames(mat) <- c(paste0("X",1:4),"obs.clones","clones","cells")
  list(eta=mat,
       lambda=cbind(lam,lamN,
		    lambdacells=colSums(etalam)*lam,
		    ExClones=colSums(etalam))
       )
}

## this needs to be corrected - at least for eta

gini <- function(lamobj,index=1,cells=3,clones=4){
  lamord <- order(lamobj[,index])
  lamobj[] <- lamobj[lamord,]
  prclone <- prop.table(lamobj[,clones])
  prcells <- prop.table(lamobj[,cells])
  wd <- prclone+c(prclone[-1L],0)
  area <-
    sum(wd*cumsum(prcells))
  area
}


minqeta <- function(eta,k=1:4){
  res <-
    lapply(k,
	   function(kelt) combn(ncol(eta)-3,kelt,
				function(xi) kelt*apply(eta[,xi,drop=FALSE],1,min)))
  cbind(matrix(unlist(res),nrow=nrow(eta)),
	eta[,-2:0+ncol(eta),drop=FALSE]) 
}



frac100 <- function(el) {
  tmp2 <- el$lambda
  tmp2 <- rbind(0,tmp2[order(-tmp2[,1]),])
  approx(cumsum(tmp2[,4]),cumsum(tmp2[,3]),100)$y/sum(tmp2[,3])
}


## this needs to work for nrow(el$eta)==1
minqgini <- function(x,k=1:4,om=NULL){
  el <- getEtaLambda(x,om)
  y <- el$eta  
  c(
    unlist(lapply(k,function(kelt){
      mins <- combn(ncol(y)-3,kelt,function(xi) apply(y[,xi,drop=FALSE],1,min))
      minkelt <- cbind(matrix(mins,nrow=nrow(y)),
		       y[,-2:0+ncol(y),drop=FALSE])
      nc <- ncol(minkelt)
      ncells <- prop.table(minkelt[,nc])
      ncells %*%  minkelt[,1:(nc-3)] * kelt
    } )),
    gini=gini(el$lambda),
    clones=sum(el$lambda[,"ExClones"]),
    frac100 = frac100(el)
  )
}

minqkeep <- function(pass,log.posterior,i,nkeep){
  om <- eval.parent(quote(om))
  if (i<nkeep)
    unname(c(minqgini(pass,1:4,om),logLik=log.posterior))
  else
    keep.default(pass,log.posterior,i,nkeep)
}


## note that this creates an environment that may use space unnecessarily

minqkeep.el <-
  function(eq4sum=0.0,lambdaSum=0.0) {
    function(pass,log.posterior,i,nkeep){
      om <- eval.parent(quote(om))
      if (i<nkeep){
	res <- unname(c(minqgini(pass, 1:4, om),logLik=log.posterior))
	eq4 <- apply(with(pass, eta[, 1:etaM]), 2, min)
	eq4sum <<- eq4sum + eq4[ 1L + pass[["dataToEta"]] ]
	lambdaSum <<- lambdaSum + with(pass, lambda[1L + dataToLambda])
	res
      } else {
	res <- keep.default(pass,log.posterior,i,nkeep)
	res[["eq4.average"]] <- eq4sum / (nkeep-1)
	res[["lambda.average"]] <- lambdaSum / (nkeep - 1)
	eq4sum <<- 0.0
	lambdaSum <<- 0.0
	res
      }
    }
  }



minqminimal <- function(pass,...){
  om <- eval.parent(quote(om))
  minqgini(pass,1:4,om)
}
load("metadata.RData")
load("rscans.RData")

redo.fit.fun <- function(obj,test=FALSE){
  orig.call <- attr(obj,"call")
  orig.om <- orig.call[["om"]]
  orig.fit <- minqgini(tail(obj,1L)[[1L]],om=eval(orig.om))
  redo.expr <- orig.call
  redo.expr$eta <-   redo.expr$etaN <-   redo.expr$etaM <-
    redo.expr$dataToEta <- redo.expr$alphaEta <-
      redo.expr$lambda <-   redo.expr$lambdaN <-   redo.expr$lambdaM <-
	redo.expr$dataToLambda <- redo.expr$alphaLambda <- NULL
  if (test){
    redo.expr[["nburn"]] <- 0L
    redo.expr[["nkeep"]] <- 5L
  }
  redo.expr[["wtab"]] <- as.symbol("wdata")
  redo.expr[["keep"]] <- as.symbol("minqminimal")
  redo.expr[["om"]] <- as.symbol("om.cur")
  ## simulate new table of counts
  wdata <- uniTab(rtab(obj,om=eval(orig.om), tol=0.0001))
  ## recalc om2 like using wdata
  om.index <- orig.om[[3]]
  ## psi2 was used in the original, so back it out
  om0 <- eval(orig.om)%*%diag(1/psi2[[om.index]])
  ## update psi2
  psi.cur <- colSums(wdata$tab*wdata$n)/
    as.numeric(dna.wts[[om.index]]$VCN)/
    pre.counts[[om.index]]
  psi.cur <- psi.cur/max(psi.cur)
  ## and now:
  om.cur <- om0%*%diag(psi.cur)
  redo.res <- do.call(cbind,eval(redo.expr))
  new.fit <- rowMeans(redo.res)
  lower <- apply(redo.res,1,quantile,probs=0.025)
  upper <- apply(redo.res,1,quantile,probs=1-0.025)
  list(orig=orig.fit,
       redo=new.fit,
       lower=lower,upper=upper)
}

redo.fits <-
  mclapply(rescans,
	   redo.fit.fun,mc.cores=min(detectCores(),50L),
	   mc.set.seed=TRUE,mc.preschedule=FALSE)

save(redo.fits,file="redo.fits.RData")
if (!interactive()) q("no")
