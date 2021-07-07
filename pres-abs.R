load("pWAS.list.RData")


pwas.time <- function(x){
  dot.for.dash <- function(x) sub("[-_]",".",x)
  d.nought <- function(x) sub("d0.*","D0",x,ignore.case=TRUE)
  d42.to.wk6 <- function(x) sub("^d4[25]$","wk6",x)
  pwas.unit <- function(x) c(WK = "WEEKS", M="MONTHS", Y="YEARS", D="DAYS")[
			     sub("[[:digit:]].*$","",x)]
  pwas.timeval <- function(x) sub("[[:alpha:]]+","",x)
  x <- as.character(x)
  tm <- toupper( d42.to.wk6( d.nought( dot.for.dash( x ))))
  res <- paste(pwas.timeval( tm ), pwas.unit( tm ))
  res[ res == "18 MONTHS" ]  <- "1.5 YEARS"
  res[ res == "30 MONTHS" ]  <- "2.5 YEARS"
  res[ res == "42 MONTHS" ]  <- "3.5 YEARS"
  res[ res == "32 DAYS" ]  <- "1 MONTHS"
  res[ res %in% c("45 DAYS", "42 DAYS") ]  <- "1.5 MONTHS"
  res[ res == "6 WEEKS" ]  <- "1.5 MONTHS"
  res
}

tabs <-
    lapply(pWAS.list[1:4],
           function(pt){
               time <- pwas.time(pt$d$timePoint)
               crossprod(
                   table(pt$d$posid,
                         factor(time,unique(time)))>0
               )})


