library(tidyverse)

pWAS.list <- readRDS('pWAS.list.rds')

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

# tabs <-
#   lapply(pWAS.list,
#          function(pt){
#            output_list = list()
#            time <- pwas.time(pt$d$timePoint)
#            table <- table(pt$d$posid,factor(time,unique(time)))>0
#            table_sum <- crossprod(table)
#            counts <- data.frame(table) %>% summarise_all(~sum(.)) %>% simplify()
#            names(counts) <-  colnames(table_sum)
#            output_list$table  <- table_sum
#            output_list$counts <- counts
#            browser()
#            return(output_list)
#            })

tabs <-
  imap(pWAS.list,
         function(pt,x){
           output_list = list()
           time <- pwas.time(pt$d$timePoint)
           table <- table(pt$d$posid,factor(time,unique(time)))>0
           table_sum <- crossprod(table)
           counts <- data.frame(table) %>% summarise_all(~sum(.)) %>% simplify()
           names(counts) <-  colnames(table_sum)
           output_list$table  <- table_sum
           output_list$counts <- counts
           output_list$name <- x
 #          browser()
           return(output_list)
         })

# > head(table(pt$d$posid,factor(time,unique(time))))
# 
# 0 DAYS 1 MONTHS 1.5 MONTHS 3 MONTHS 6 MONTHS 1 YEARS 1.5 YEARS 2 YEARS 2.5 YEARS 3 YEARS 3.5 YEARS 4 YEARS 4.5 YEARS 5 YEARS
# chr1_GL383518v1_alt-145344      0        0          0        0        0       1         0       0         0       0         0       0         0       0
# chr1_GL383518v1_alt-145817      0        0          0        0        1       0         0       0         0       0         0       0         0       0
# chr1_GL383518v1_alt-146914      0        0          0        0        0       0         0       1         0       2         2       0         0       0
# chr1_GL383518v1_alt-147164      0        0          1        0        0       0         0       0         0       0         0       0         0       0
# chr1_GL383518v1_alt-82897       0        0          0        0        0       0         0       0         0       1         2       0         0       0
# chr1_GL383518v1_alt+133749      0        0          0        0        0       0         0       0         1       0         0       0         0       0




#t1 <- tabs$pWAS02$table
#t1_c <- tabs$pWAS02$counts
#rownames(t1)
#as.data.frame(t1) %>% rownames_to_column(var = "rowname") %>% mutate_if(is_numeric,~./sum(.))

# factor_order <- c("0 DAYS","1 MONTHS","1.5 MONTHS","3 MONTHS","6 MONTHS",
#                   "1 YEARS","1.5 YEARS","2 YEARS","2.5 YEARS","3 YEARS",
#                   "3.5 YEARS","4 YEARS","4.5 YEARS","6 YEARS")
# 
# 
# t1_df <- as.data.frame(t1) %>%
#   rownames_to_column(var = "y") %>% 
# #  mutate(across(where(is.numeric),~./max(.))) %>%
#   pivot_longer(!y,values_to='val', names_to='x') %>%
#   mutate(x=fct_relevel(x,factor_order[factor_order %in% x])) %>%
#   mutate(y=fct_relevel(as.factor(y),factor_order[factor_order %in% y])) %>%
#   mutate(y_count=t1_c[as.character(y)],x_count=t1_c[as.character(x)]) %>%
#   mutate(jaccard=(val/(x_count+y_count-val)))

#  fct_relevel(as.factor(t1_df$y),factor_order[factor_order %in% t1_df$y])

long_tables <- imap(tabs,function(tab,name){
  #browser()
  factor_order <- c("0 DAYS","1 MONTHS","1.5 MONTHS","3 MONTHS","6 MONTHS","9 MONTHS",
                    "1 YEARS","1.5 YEARS","2 YEARS","2.5 YEARS","3 YEARS",
                    "3.5 YEARS","4 YEARS","4.5 YEARS","5 YEARS","6 YEARS")
  as.data.frame(tab$table) %>%
    rownames_to_column(var = "y") %>% 
    #  mutate(across(where(is.numeric),~./max(.))) %>%
    pivot_longer(!y,values_to='val', names_to='x') %>%
    mutate(x=fct_relevel(x,factor_order[factor_order %in% x])) %>%
    mutate(y=fct_relevel(as.factor(y),factor_order[factor_order %in% y])) %>%
    mutate(y_count=tab$counts[as.character(y)],x_count=tab$counts[as.character(x)]) %>%
    mutate(jaccard=(val/(x_count+y_count-val))) %>%
    mutate(sample_name=name) %>%
    mutate(x_s=paste(x,sample_name),y_s=paste(y,sample_name))
})

factor_order <- c("0 DAYS","1 MONTHS","1.5 MONTHS","3 MONTHS","6 MONTHS","9 MONTHS",
                  "1 YEARS","1.5 YEARS","2 YEARS","2.5 YEARS","3 YEARS",
                  "3.5 YEARS","4 YEARS","4.5 YEARS","5 YEARS","6 YEARS")
long_names <- names(long_tables)
SS <- levels(as.factor(reduce(long_tables,rbind)$x_s))
is_valid <- function(x,y,ss) paste(x,y) %in% ss
valid_paste <- function(x,y,ss) ifelse(is_valid(x,y,ss),paste(x,y),'Na')

factor_order_plus <- as.vector(outer(factor_order, long_names, FUN = "paste"))
factor_order_plus_bis <- as.vector(outer(factor_order, long_names, FUN =valid_paste,ss=SS))

inv_paste <- function(x,y) paste(y,x)
factor_order_per_date <- as.vector(outer(long_names,factor_order, FUN = inv_paste))

get_first <- function(x,y) return(x)
get_second <- function(x,y) return(y)
date_order_per_date <- as.vector(outer(long_names,factor_order, FUN = get_second))
names_order_per_date <- as.vector(outer(long_names,factor_order, FUN = get_first))
per_date_sample_color <- sapply(as.integer(as.factor(names_order_per_date)),function(x) {return(hcl.colors(5, "viridis")[x])} )
sample_to_color_df <- data.frame(sample=names(long_tables),color=hcl.colors(5, "viridis"))

names_order_per_sample <- as.vector(outer(factor_order, long_names, FUN = get_second))
date_order_per_sample <- as.vector(outer(factor_order, long_names, FUN = get_first))
per_sample_sample_color <-sapply(as.integer(as.factor(names_order_per_sample)),function(x) {return(hcl.colors(5, "viridis")[x])} )

pp <- reduce(long_tables,rbind) %>% 
  mutate(x_s=fct_relevel(x_s,factor_order_per_date[factor_order_per_date %in% x_s])) %>%
  mutate(y_s=fct_relevel(y_s,factor_order_per_date[factor_order_per_date %in% y_s])) %>%
  ggplot( aes(x_s,y_s, fill= jaccard,color=sample_name)) + 
  geom_tile(size=0) +
  scale_fill_binned(breaks = c(0.001, 0.2, 0.4,0.6),low="white", high="blue",show.limits = TRUE) +
  #scale_fill_gradient(low="white", high="black",na.value = 'salmon') +
  #  geom_text(aes(label=round(jaccard,digits = 2)),size=1, color='red') +
  theme_classic() +
  scale_x_discrete(labels=date_order_per_date) +
  scale_y_discrete(labels=date_order_per_date) +
  #  theme(panel.background=element_rect(fill="grey95", colour="grey95")) +
  theme(axis.text.x = element_text(angle = 90,color = per_date_sample_color),
        axis.text.y = element_text(color = per_date_sample_color),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  #  xlab("Sample") +
  #  labs(fill='Dilution \nat threshold',
  #       title = "Quiescent") +
  #  My_Theme +
  coord_equal()


png(height = 10, width = 10,units = 'in', res=300, file = 'test3.png')
pp + scale_colour_manual(name = 'the colour', values =c('pWAS02'='black','pWAS03'='red','pWAS04'='black','pWAS05 '='red','pWAS06'='black'),
                         labels =sample_to_color_df$sample ,guide = 'legend') +
  guides(colour=guide_legend("Subject", override.aes=list(shape=22,size=0, fill=sample_to_color_df$color)))

dev.off()

pp2 <- reduce(long_tables,rbind) %>% 
  mutate(x_s=fct_relevel(x_s,factor_order_plus[factor_order_plus %in% x_s])) %>%
  mutate(y_s=fct_relevel(y_s,factor_order_plus[factor_order_plus %in% y_s])) %>%
  ggplot( aes(x_s,y_s, fill= jaccard,color=sample_name)) + 
  geom_tile(size=0) +
  scale_fill_binned(breaks = c(0.001, 0.2, 0.4,0.6),low="white", high="blue",show.limits = TRUE) +
  #scale_fill_gradient(low="white", high="black",na.value = 'salmon') +
  #  geom_text(aes(label=round(jaccard,digits = 2)),size=1, color='red') +
  theme_classic() +
  scale_x_discrete(labels=date_order_per_sample) +
  scale_y_discrete(labels=date_order_per_sample) +
  #  theme(panel.background=element_rect(fill="grey95", colour="grey95")) +
  theme(axis.text.x = element_text(angle = 90,color = per_sample_sample_color),
        axis.text.y = element_text(color = per_sample_sample_color),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  #  xlab("Sample") +
  #  labs(fill='Dilution \nat threshold',
  #       title = "Quiescent") +
  #  My_Theme +
  coord_equal()

png(height = 10, width = 10,units = 'in', res=300, file = 'test2.png')
pp2 + scale_colour_manual(name = 'the colour', values =c('pWAS02'='black','pWAS03'='red','pWAS04'='black','pWAS05 '='red','pWAS06'='black'),
                        labels =sample_to_color_df$sample ,guide = 'legend') +
  guides(colour=guide_legend("Subject", override.aes=list(shape=22,size=0, fill=sample_to_color_df$color)))
dev.off()

