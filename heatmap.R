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

original_names <- c("pWAS02","pWAS03","pWAS04","pWAS05","pWAS06")
new_names <- c('Pt 1','Pt 2','Pt 4','Pt 3','Pt 5')
names(new_names) <- original_names
names(pWAS.list) <- new_names[names(pWAS.list)]
pWAS.list <- pWAS.list[c('Pt 1','Pt 2','Pt 3','Pt 4','Pt 5')]


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
#                            0 DAYS   1 MONTHS 1.5 MONTHS 3 MONTHS 6 MONTHS 1 YEARS   1.5 YEARS  2 YEARS 2.5 YEARS 3 YEARS 3.5 YEARS  4 YEARS  4.5 YEARS   5 YEARS
# chr1_GL383518v1_alt-145344      0        0          0         0        0       1         0       0         0       0         0       0         0         0
# chr1_GL383518v1_alt-145817      0        0          0         0        1       0         0       0         0       0         0       0         0         0
# chr1_GL383518v1_alt-146914      0        0          0         0        0       0         0       1         0       2         2       0         0         0
# chr1_GL383518v1_alt-147164      0        0          1         0        0       0         0       0         0       0         0       0         0         0
# chr1_GL383518v1_alt-82897       0        0          0         0        0       0         0       0         0       1         2       0         0         0
# chr1_GL383518v1_alt+133749      0        0          0         0        0       0         0       0         1       0         0       0         0         0


long_tables <- imap(tabs,function(tab,name){
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

#factor_order_plus <- as.vector(outer(factor_order, long_names, FUN = "paste"))
factor_order_plus_tmp <- as.vector(outer(factor_order, long_names, FUN =valid_paste,ss=SS))
factor_order_plus <- factor_order_plus_tmp[factor_order_plus_tmp != 'Na']

inv_valid_paste <- function(x,y,ss) ifelse(is_valid(y,x,ss),paste(y,x),'Na')
factor_order_per_date_tmp <- as.vector(outer(long_names,factor_order, FUN = inv_valid_paste, ss=SS))
factor_order_per_date <- factor_order_per_date_tmp[factor_order_per_date_tmp != 'Na']

get_first_valid <- function(x,y,ss) ifelse(is_valid(y,x,ss),x,'Na')
get_second_valid <- function(x,y,ss) ifelse(is_valid(y,x,ss),y,'Na')

date_order_per_date_tmp <- as.vector(outer(long_names,factor_order, FUN = get_second_valid,ss=SS))
date_order_per_date <- date_order_per_date_tmp[date_order_per_date_tmp != 'Na']


names_order_per_date_tmp <- as.vector(outer(long_names,factor_order, FUN = get_first_valid,ss=SS))
names_order_per_date <- names_order_per_date_tmp[names_order_per_date_tmp != 'Na']
# old viridis text
#per_date_sample_color <- sapply(as.integer(as.factor(names_order_per_date)),function(x) {return(hcl.colors(5, "viridis")[x])} )
#sample_to_color_df <- data.frame(sample=names(long_tables),color=hcl.colors(5, "viridis"))
#div_text_color <- c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e')
div_text_color <- c('blue','red','green','orange','purple')
per_date_sample_color <- sapply(as.integer(as.factor(names_order_per_date)),function(x) {return(div_text_color[x])} )
sample_to_color_df <- data.frame(sample=names(long_tables),color=div_text_color)
#========

get_first_valid_inv <- function(x,y,ss) ifelse(is_valid(x,y,ss),x,'Na')
get_second_valid_inv <- function(x,y,ss) ifelse(is_valid(x,y,ss),y,'Na')
names_order_per_sample_tmp <- as.vector(outer(factor_order, long_names, FUN = get_second_valid_inv,ss=SS))
names_order_per_sample <- names_order_per_sample_tmp[names_order_per_sample_tmp!= 'Na']

date_order_per_sample_tmp <- as.vector(outer(factor_order, long_names, FUN = get_first_valid_inv,ss=SS))
date_order_per_sample <- date_order_per_sample_tmp[date_order_per_sample_tmp!= 'Na']
#per_sample_sample_color <-sapply(as.integer(as.factor(names_order_per_sample)),function(x) {return(hcl.colors(5, "viridis")[x])} )
per_sample_sample_color <-sapply(as.integer(as.factor(names_order_per_sample)),function(x) {return(div_text_color[x])} )

pp <- reduce(long_tables,rbind) %>% 
  mutate(x_s=fct_relevel(x_s,factor_order_per_date)) %>%
  mutate(y_s=fct_relevel(y_s,factor_order_per_date)) %>%
  ggplot( aes(x_s,y_s, fill= jaccard,color=sample_name)) + 
  geom_tile(size=0) +
  scale_fill_stepsn(breaks = c(0.001,0.1, 0.25, 0.5),colors=c('#ffeda0','#feb24c','#f03b20','#000000'),values = c(0.05,0.25,0.5,0.8),show.limits = TRUE,na.value="grey90") +
  #scale_fill_stepsn(breaks = c(0.001, 0.25, 0.5),colors=c('#ffeda0','#feb24c','#f03b20','#000000'),values = c(0,0.25,0.5,0.8),show.limits = TRUE) +
  #scale_fill_binned(breaks = c(0.001, 0.25, 0.5),colors=c('#ffeda0','#feb24c','#f03b20','#000000'),values = c(0,0.25,0.5,0.8),show.limits = TRUE) +
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


png(height = 10, width = 10,units = 'in', res=300, file = 'per_date_jaccard.png')
pp + scale_colour_manual(name = 'the colour', values =c('pWAS02'='black','pWAS03'='red','pWAS04'='black','pWAS05 '='red','pWAS06'='black'),
                         labels =sample_to_color_df$sample ,guide = 'legend') +
  guides(colour=guide_legend("Subject", override.aes=list(shape=22,size=0, fill=sample_to_color_df$color)))

dev.off()

pp2 <- reduce(long_tables,rbind) %>% 
  mutate(x_s=fct_relevel(x_s,factor_order_plus)) %>%
  mutate(y_s=fct_relevel(y_s,factor_order_plus)) %>%
  ggplot( aes(x_s,y_s, fill= jaccard,color=sample_name)) + 
  geom_tile(size=0) +
  scale_fill_stepsn(breaks = c(0.001,0.1, 0.25, 0.5),colors=c('#ffeda0','#feb24c','#f03b20','#000000'),values = c(0.05,0.25,0.5,0.8),show.limits = TRUE,na.value="grey90") +
  #scale_fill_binned(breaks = c(0.001, 0.2, 0.4,0.6),low="white", high="blue",show.limits = TRUE) +
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

png(height = 10, width = 10,units = 'in', res=300, file = 'per_subject_jaccard.png')
pp2 + scale_colour_manual(name = 'the colour', values =c('pWAS02'='black','pWAS03'='red','pWAS04'='black','pWAS05 '='red','pWAS06'='black'),
                        labels =sample_to_color_df$sample ,guide = 'legend') +
  guides(colour=guide_legend("Subject", override.aes=list(shape=22,size=0, fill=sample_to_color_df$color)))
dev.off()

# charles figures ----------------- 

ch.final.table <- readRDS('charles_table.rds')

# per sample charles -------------------------

factor_order_ch_tmp <- as.vector(outer(factor_order, long_names, FUN =valid_paste,ss=ch.final.table$y_s))
factor_order_ch <- factor_order_plus_tmp[factor_order_ch_tmp != 'Na']

names_order_ch_sample_tmp <- as.vector(outer(factor_order, long_names, FUN = get_second_valid_inv,ss=ch.final.table$y_s))
names_order_ch_sample <- names_order_ch_sample_tmp[names_order_ch_sample_tmp!= 'Na']

date_order_ch_sample_tmp <- as.vector(outer(factor_order, long_names, FUN = get_first_valid_inv,ss=ch.final.table$y_s))
date_order_ch_sample <- date_order_ch_sample_tmp[date_order_ch_sample_tmp!= 'Na']
#per_sample_sample_color <-sapply(as.integer(as.factor(names_order_per_sample)),function(x) {return(hcl.colors(5, "viridis")[x])} )
ch_sample_color <-sapply(as.integer(as.factor(names_order_ch_sample)),function(x) {return(div_text_color[x])} )


pp3 <-  ch.final.table %>% #mutate(val=val/100) %>%
  mutate(x_s=fct_relevel(x_s,factor_order_ch)) %>%
  mutate(y_s=fct_relevel(y_s,factor_order_ch)) %>%
  ggplot( aes(x_s,y_s, fill=val,color=sample_name)) + 
  geom_tile(size=0) +
  #scale_fill_gradientn(breaks = c(0.001,0.1, 0.25, 0.5),colors=c('#ffeda0','#feb24c','#f03b20','#000000'),values = c(0,0.5,0.95,1),na.value="grey90") +
  scale_fill_viridis_c(direction = -1,name="Shared Top100\nClones") +
  theme_classic() +
  scale_x_discrete(labels=date_order_ch_sample) +
  scale_y_discrete(labels=date_order_ch_sample) +
  #  theme(panel.background=element_rect(fill="grey95", colour="grey95")) +
  theme(axis.text.x = element_text(angle = 90,color = ch_sample_color, size = 10),
        axis.text.y = element_text(color = ch_sample_color, size = 10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  coord_equal()

png(height = 4, width = 5,units = 'in', res=300, file = 'per_sample_ch.png')
pp3 + scale_colour_manual(name = 'the colour', values =c('pWAS02'='black','pWAS03'='red','pWAS04'='black','pWAS05 '='red'),
                        labels =sample_to_color_df$sample[1:4] ,guide = 'legend') +
  #  guides(fill=guide_legend("Shared Top100\nClones")) +
  guides(colour=guide_legend("Subject", override.aes=list(shape=22,size=0, fill=sample_to_color_df$color[1:4])))
dev.off()

# per date charles ---------------
factor_order_ch_date_tmp <- as.vector(outer(long_names,factor_order, FUN = inv_valid_paste, ss=ch.final.table$y_s))
factor_order_ch_date <- factor_order_ch_date_tmp[factor_order_ch_date_tmp != 'Na']

date_order_ch_date_tmp <- as.vector(outer(long_names,factor_order, FUN = get_second_valid,ss=ch.final.table$y_s))
date_order_ch_date <- date_order_ch_date_tmp[date_order_ch_date_tmp != 'Na']

names_order_ch_date_tmp <- as.vector(outer(long_names,factor_order, FUN = get_first_valid,ss=ch.final.table$y_s))
names_order_ch_date <- names_order_ch_date_tmp[names_order_ch_date_tmp != 'Na']

ch_date_color <-sapply(as.integer(as.factor(names_order_ch_date)),function(x) {return(div_text_color[x])} )

pp4 <-  ch.final.table %>% #mutate(val=val/100) %>%
  mutate(x_s=fct_relevel(x_s,factor_order_ch_date)) %>%
  mutate(y_s=fct_relevel(y_s,factor_order_ch_date)) %>%
  ggplot( aes(x_s,y_s, fill=val,color=sample_name)) + 
  geom_tile(size=0) +
  #scale_fill_gradientn(breaks = c(0.001,0.1, 0.25, 0.5),colors=c('#ffeda0','#feb24c','#f03b20','#000000'),values = c(0,0.5,0.95,1),na.value="grey90") +
  scale_fill_viridis_c(direction = -1,name="Shared Top100\nClones") +
  theme_classic() +
  scale_x_discrete(labels=date_order_ch_date) +
  scale_y_discrete(labels=date_order_ch_date) +
  #  theme(panel.background=element_rect(fill="grey95", colour="grey95")) +
  theme(axis.text.x = element_text(angle = 90,color = ch_date_color,size = 10),
        axis.text.y = element_text(color = ch_date_color,size = 10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  coord_equal()

png(height = 4, width = 5,units = 'in', res=300, file = 'per_date_ch.png')
pp4 + 
  scale_colour_manual(name = 'the colour', values =c('pWAS02'='black','pWAS03'='red','pWAS04'='black','pWAS05 '='red'),
                          labels =sample_to_color_df$sample[1:4] ,guide = 'legend') +
#  guides(fill=guide_legend("Shared Top100\nClones")) +
  guides(colour=guide_legend("Subject", override.aes=list(shape=22,size=0, fill=sample_to_color_df$color[1:4])))
dev.off()

