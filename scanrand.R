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
randscan <- function(wtab,om,...){
  mc <- match.call()
  wt.order <- sample(length(wtab$data.index))
  wtab$data.index <- wtab$data.index[wt.order]
  res <- gibbsScan(wtab,om,...)
  res$dataToEta[wt.order] <-  res$dataToEta
  res$dataToLambda[wt.order] <-  res$dataToLambda
  rescall <- attr(res,"call")
  mc <- mc[-1]
  rescall[names(mc)] <- mc
  attr(res,"call") <- rescall
  res
}
library(parallel)
library(cellTypeCompositions)
load("metadata.RData")
rscan.calls <- list()
for (i in seq_along(wttabs))
  rscan.calls[[ i ]] <-
    bquote(randscan(wttabs[[.(j)]],om2[[.(j)]],lambdaShape=1.0,lambdaRate=0.0001),list(j=i))
rscans <- mclapply(rscan.calls,function(x) replicate(10L, eval(x), simplify=FALSE),
		   mc.set.seed=TRUE,mc.cores=min(detectCores(),18L))
rescans <- mclapply(unlist(rscans,recursive=FALSE),
		    update, keep=minqkeep,nburn=300L,nkeep=101L,nthin=20L,
		    mc.preschedule=FALSE,
		    mc.set.seed=TRUE,mc.cores=min(detectCores(),40L))
save(rscans,rescans,file="rscans.RData")
if (!interactive()) q("no")
