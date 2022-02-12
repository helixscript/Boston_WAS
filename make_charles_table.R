#charles tables



#2 YEARS	3 YEARS 	4 YEARS 	5 YEARS
data.frame(
'2 YEARS'=c(100,100,100,100),
'3 YEARS'=c(99,100,100,100),
'4 YEARS'=c(93,97,100,100),
'5 YEARS'=c(92,96,100,100))

was02_names <- c('2 YEARS','3 YEARS','4 YEARS','5 YEARS')
m02 <- matrix(c(100,100,100,100,99,100,100,100,93,97,100,100,92,96,100,100),
       ncol = 4,byrow=TRUE,dimnames=list(was02_names,was02_names))

was03_names <- c('9 MONTHS','1 YEARS','2 YEARS')
m03 <- matrix(c(100,99,96,92,100,90,81,90,100),
       ncol = 3,byrow=TRUE,dimnames=list(was03_names,was03_names))


was04_names <- c('6 MONTHS','1 YEARS','1.5 YEARS','2 YEARS','3 YEARS')
m04 <- matrix(c(100,33,26,33,36,34,100,45,49,52,35,55,100,99,98,25,42,90,100,99,19,31,81,89,100),
       ncol = 5,byrow=TRUE,dimnames=list(was04_names,was04_names))

was05_names <- c('1 MONTHS','3 MONTHS','1 YEARS','2 YEARS','3 YEARS')
m05 <- matrix(c(100,2,0,0,0,0,100,42,55,36,0,35,100,100,100,0,25,97,100,99,0,15,89,96,100),
       ncol = 5,byrow=TRUE,dimnames=list(was05_names,was05_names))


as.data.frame(m05) %>%
  rownames_to_column(var = "y") %>% 
  pivot_longer(!y,values_to='val', names_to='x')


ch.list <-
  list('Pt 1'=list(table=m02),
     'Pt 2'=list(table=m03),
     'Pt 4'=list(table=m04),
     'Pt 3'=list(table=m05))

long_tables <- imap(ch.list,function(tab,name){
  factor_order <- c("0 DAYS","1 MONTHS","1.5 MONTHS","3 MONTHS","6 MONTHS","9 MONTHS",
                    "1 YEARS","1.5 YEARS","2 YEARS","2.5 YEARS","3 YEARS",
                    "3.5 YEARS","4 YEARS","4.5 YEARS","5 YEARS","6 YEARS")
  as.data.frame(tab$table) %>%
    rownames_to_column(var = "y") %>% 
    pivot_longer(!y,values_to='val', names_to='x') %>%
    mutate(x=fct_relevel(x,factor_order[factor_order %in% x])) %>%
    mutate(y=fct_relevel(as.factor(y),factor_order[factor_order %in% y])) %>%
    mutate(sample_name=name) %>%
    mutate(x_s=paste(x,sample_name),y_s=paste(y,sample_name))
})
ch.final.table <- reduce(long_tables,rbind)
saveRDS(ch.final.table,file='charles_table.rds')
