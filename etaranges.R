maxRangeQuantile <-
   function(wtab,om,etaBreaks=c(0.025,0.975),
 	   qbreaks=c(0.01,0.02,0.05),nkeep=100,...){
   tmptab <- uniTab(wtab$tab)
   tmpscan <- gibbsScan(tmptab,om)
   tmpscan2 <- tmpscan
   tmpscan2[[1]]$eta <- with(tmpscan2[[1]],eta[,dataToEta])
   tmpscan2[[1]]$etaN <- with(tmpscan2[[1]],rep(1L,length(dataToEta)))
   tmpscan2[[1]]$dataToEta <- with(tmpscan2[[1]],1:length(dataToEta))
 
   tmpscan2[[1]]$lambda <- with(tmpscan2[[1]],lambda[dataToLambda])
   tmpscan2[[1]]$lambdaN <- with(tmpscan2[[1]],rep(1L,length(dataToLambda)))
   tmpscan2[[1]]$dataToLambda <- with(tmpscan2[[1]],1:length(dataToLambda))
 
   tmpscan2[[1]]$lambdaM <- with(tmpscan2[[1]],length(dataToLambda))
   tmpscan2[[1]]$etaM <- with(tmpscan2[[1]],length(dataToLambda))
 
   tmp2 <- update(tmpscan2,ijvals=nrow(tmptab$tab),nkeep=nkeep,...)
   tmpeta <- sapply(tmp2,function(x) as.vector(x$eta))
   etasd <- apply(tmpeta,1,sd)
   etaRange <- apply(tmpeta,1,function(x) diff(quantile(x,etaBreaks)))
   maxsd <- apply(matrix(etasd,nr=4),2,max)
   maxRange <- apply(matrix(etaRange,nr=4),2,max)
   N.est <- vegan::estimateR(with(wtab,rep(rowSums(tab),n)))[1:2]
   list(range =  quantile(rep(maxRange,wtab$n),qbreaks * N.est[2] / N.est[1]),
        sd = 
 	 quantile(rep(maxsd,wtab$n), qbreaks * N.est[2] / N.est[1]))
   }
library(parallel)
library(cellTypeCompositions)
load("metadata.RData")

rscan.calls <- list()
for (i in seq_along(wttabs))
  rscan.calls[[ i ]] <-
    bquote(maxRangeQuantile(wttabs[[.(j)]],om2[[.(j)]],nkeep=1000,nthin=5),list(j=i))

etaRanges <-
  mclapply(rscan.calls, eval, 
	   mc.set.seed=TRUE,mc.cores=min(detectCores(),18L))

save(etaRanges,file="etaRanges.RData")
if (!interactive()) q("no")
