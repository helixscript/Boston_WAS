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
library(parallel)
library(cellTypeCompositions)
load("rscans.RData")
load("metadata.RData")

tiltom <- function(omelt,plus=0.01){
  prop.om <- prop.table(omelt,1)
  prop.om <- prop.om - diag(4)*4*plus + plus
  new.om <- prop.om * rowSums(omelt)
  new.om
}

om3 <- lapply(om2,tiltom)

rescan.calls <-
  lapply(rescans,
	 function(x) {
	   gcall <- attr(x,"call")
	   gcall[["om"]][[2]] <- as.symbol("om3")
	   gcall[4:13] <- NULL
	   gcall
	 })
gseq3 <- mclapply(rescan.calls, eval, mc.cores=min(50L,detectCores()))
gscan.eqs <- sapply(gseq3,function(x) unlist(x[-101]))
dim(gscan.eqs) <- c(19,100,10,18)
g3.means <- apply(gscan.eqs,c(1,4),mean)
g3.upper <- apply(gscan.eqs,c(1,4),quantile,probs=0.975)
g3.lower <- apply(gscan.eqs,c(1,4),quantile,probs=0.025)
save(g3.means,g3.upper,g3.lower,file="gseq3.RData")
if (!interactive()) q("no")
