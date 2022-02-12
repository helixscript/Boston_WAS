library(dplyr)
library(ggplot2)
library(RMySQL)
library(gt23)
library(GenomicRanges)
library(vegan)
library(untb)
library(egg)
library(gridtext)

invisible(sapply(dbListConnections(MySQL()), dbDisconnect))
dbConn  <- dbConnect(MySQL(), group='specimen_management')
samples <- dbGetQuery(dbConn, 'select * from gtsp where Trial="WAS_Boston"')

if(! file.exists('intSites.rds')){
  intSites <- gt23::getDBgenomicFragments(samples$SpecimenAccNum, 'specimen_management', 'intsites_miseq') %>%
              gt23::stdIntSiteFragments(CPUs = 2) %>%
              gt23::collapseReplicatesCalcAbunds() %>%
              gt23::annotateIntSites(CPUs = 2)
  
  saveRDS(intSites, 'intSites.rds')
} else {
  intSites <- readRDS('intSites.rds')
}

intSites$patient   <- sub('^p', '', intSites$patient)
intSites$timePoint <- sub('D32', 'M1',   intSites$timePoint)
intSites$timePoint <- sub('D42', 'M1.5', intSites$timePoint)
intSites$timePoint <- sub('D45', 'M1.5', intSites$timePoint)
intSites$timePoint <- sub('M18', 'Y1.5', intSites$timePoint)
intSites$timePoint <- sub('M30', 'Y2.5', intSites$timePoint)
intSites$timePoint <- sub('M42', 'Y3.5', intSites$timePoint)

# Tally the numbers of time points per patient / cell type
# to determine which cell types can be readily plotted.
# Whole blood / PBMC, T cells, Neutrophils, NK cells, B cells
d <- group_by(data.frame(intSites), patient, cellType) %>%
     summarise(nTimePoints = n_distinct(timePoint)) %>%
     ungroup() %>%
     arrange(patient, desc(nTimePoints)) 


original_names <- c("WAS00002","WAS00003","WAS00004","WAS00005","WAS00006")
new_names <- c('Pt 1','Pt 2','Pt 4','Pt 3','Pt 5')
names(new_names) <- original_names
#names(pWAS.list) <- new_names[names(pWAS.list)]


d <- group_by(data.frame(intSites), patient, cellType, timePoint) %>%
       mutate(patient=new_names[patient]) %>%
       summarise(Chao1   = round(estimateR(estAbund, index='chao')[2], 0),
                 Shannon = diversity(estAbund),
                 Simpson = 1-diversity(estAbund, index = "simpson"),
                 Simpson_nr = simpson(estAbund, with.replacement=FALSE),
                 Simpson_inv = diversity(estAbund, index = "invsimpson"),
                 nSites  = n_distinct(posid)) %>%
     ungroup() %>%
     filter(cellType %in% c("Whole blood", "T cells", "Neutrophils", "NK cells", "PBMC", "B cells"))



d$timePoint <- factor(d$timePoint, levels = unique(sort(unique(d$timePoint))))
d$patient   <- factor(d$patient, levels = unique(sort(unique(d$patient))))


colors_np <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_distinct(d$patient))
colors <- c('blue','red','green','orange','purple')
ShannonPatientPlot <- function(x){ 
  ggplot(x, aes(timePoint, Shannon, color = patient, group = patient)) +
  theme_bw()+
  scale_color_manual(name = 'Patient', values = colors) +
  geom_point(size = 3) +
  geom_line() +
  scale_y_continuous(breaks=seq(0, 10, by = 0.5), limits=c(4,10))+
  labs(x = 'Time point', y = 'Shannon index') +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
}

shannonPatientPlotWholeBlood <- ShannonPatientPlot(filter(d, cellType == 'Whole blood', nSites >= 100))
shannonPatientPlotPBMC <- ShannonPatientPlot(filter(d, cellType == 'PBMC', nSites >= 100))

ggsave(shannonPatientPlotWholeBlood, file = 'shannonPatientPlotWholeBlood.pdf', width = 10, units = 'in', useDingbats = FALSE)
ggsave(shannonPatientPlotPBMC,       file = 'shannonPatientPlotPBMC.pdf', width = 10, units = 'in', useDingbats = FALSE)


cellTypeShannonPlot <-
  ggplot(filter(d, nSites >= 100, cellType != 'Whole blood'), aes(timePoint, Shannon, color = cellType, group = cellType)) +
    theme_bw()+
    scale_color_manual(name = 'Cell type', values = colors_np) +
    geom_point(size = 3) +
    geom_line() +
    scale_y_continuous(breaks=seq(0, 10, by = 1.0))+
    ylim(c(3,10)) +
    labs(x = 'Time point', y = 'Shannon index') +
    facet_grid(patient~.) +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(cellTypeShannonPlot, file = 'cellTypeShannonPlot.pdf', width = 10, height = 7, units = 'in', useDingbats = FALSE)


# aqui empiezan los plots de simson index por patient---------------
SimpsonPatientPlot <- function(x){ 
  ggplot(x, aes(timePoint, Simpson_nr, color = patient, group = patient)) +
    theme_bw()+
    scale_color_manual(name = 'Patient', values = colors) +
    geom_point(size = 3) +
    geom_line() +
    scale_y_continuous(breaks=seq(0, 0.8, by = 0.1), limits=c(0,0.8))+
     labs(x = 'Time point', y = 'Simpson index') +
     theme(axis.text = element_text(size = 14),
           axis.title = element_text(size = 16),
           panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))
}

SimpsonPatientPlotWholeBlood <- SimpsonPatientPlot(filter(d, cellType == 'Whole blood', nSites >= 100))
SimpsonPatientPlotPBMC <- SimpsonPatientPlot(filter(d, cellType == 'PBMC', nSites >= 100))

ggsave(SimpsonPatientPlotWholeBlood, file = 'SimpsonPatientPlotWholeBlood.png', width = 10, units = 'in')
ggsave(SimpsonPatientPlotPBMC,       file = 'SimpsonPatientPlotPBMC.png', width = 10, units = 'in')

# aqui es por typo de celula -----------------
cellTypeSimpsonPlot <-  ggplot(filter(d, nSites >= 100, cellType != 'Whole blood'), aes(timePoint, Simpson_nr, color = cellType, group = cellType)) +
  theme_bw()+
  scale_color_manual(name = 'Cell type', values = colors_np) +
  geom_point(size = 3) +
  geom_line() +
  scale_y_continuous(breaks=seq(0, 0.8, by = 0.2), limits = c(0,0.8))+
#  ylim(c(3,10)) +
  labs(x = 'Time point', y = 'Simpson index') +
  facet_grid(patient~.) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        strip.background = element_blank(),
        strip.text.y = element_blank()
        )
cellTypeSimpsonPlot

dd <- data.frame(timePoint=c('M1','M1','M1','M1','M1'),
                 Simpson_nr = c(0.7,0.7,0.7,0.7,0.7),
                 cellType=c('PBMC','PBMC','PBMC','PBMC','PBMC'),
                 patient = c('Pt 1','Pt 2','Pt 3','Pt 4','Pt 5')
)

grob_size <- 20
grid_padd <- unit(c(3, 2, 1, 2), "pt")

grid_labels <- list(
  textbox_grob('Pt 1',width=NULL,height = NULL,padding = grid_padd, gp = gpar(fontsize = grob_size,col=colors[1]), box_gp = gpar(col = "black", fill = "grey92")),
  textbox_grob('Pt 2',width=NULL,height = NULL,padding = grid_padd, gp = gpar(fontsize = grob_size,col=colors[2]), box_gp = gpar(col = "black", fill = "grey92")),
  textbox_grob('Pt 3',width=NULL,height = NULL,padding = grid_padd, gp = gpar(fontsize = grob_size,col=colors[3]), box_gp = gpar(col = "black", fill = "grey92")),
  textbox_grob('Pt 4',width=NULL,height = NULL,padding = grid_padd, gp = gpar(fontsize = grob_size,col=colors[4]), box_gp = gpar(col = "black", fill = "grey92")),
  textbox_grob('Pt 5',width=NULL,height = NULL,padding = grid_padd, gp = gpar(fontsize = grob_size,col=colors[5]), box_gp = gpar(col = "black", fill = "grey92"))
)
#dd$grob <- grob
dd$grid_l <- grid_labels 

cellTypeSimpsonPlot_ll <- cellTypeSimpsonPlot + geom_custom(data = dd, aes(data = grid_l), grob_fun = identity) 


ggsave(cellTypeSimpsonPlot_ll, file = 'cellTypeSimpsonPlot.pdf', width = 10, height = 7, units = 'in')

# aqui acaban los plots ------------


# Add nearest feature flags.
d <- data.frame(intSites) %>%
  mutate(patient=new_names[patient]) %>%
     mutate(labeledNearestFeature = paste0(nearestFeature, ' ')) %>% 
     mutate(labeledNearestFeature = ifelse(inFeature, paste0(labeledNearestFeature, '*'), labeledNearestFeature)) 

if('nearestOncoFeatureDist' %in% names(d))
  d <- mutate(d, labeledNearestFeature = ifelse(abs(nearestOncoFeatureDist) <= 50000, paste0(labeledNearestFeature, '~'), labeledNearestFeature))

if('nearestlymphomaFeatureDist' %in% names(d))
  d <- mutate(d, labeledNearestFeature = ifelse(abs(nearestlymphomaFeatureDist) <= 50000, paste0(labeledNearestFeature, '!'), labeledNearestFeature)) 

d$posidLabel <- paste0(d$posid, '\n', d$labeledNearestFeature)


numClones <- 10

# Create data frame needed to generate relative abundance plots.
d <- filter(d, cellType == 'PBMC')
abundantClones <- bind_rows(lapply(split(d, paste(d$patient, d$cellType)), function(x){
  
  x$totalCells <- sum(x$estAbund)
  
  # Add one more clone to the patient with the most sites to balance out the figure legend.
  if (x$patient == 'Pt 1') numClones <- numClones+1 
  
  # Adjust the number of clones to return based on the number of sites per cell type.
  if(nrow(x) < numClones) numClones <- nrow(x)
  
  # Sort nearest genes by abundance.
  x <- x[order(x$estAbund, decreasing = TRUE),]
  
  # Select clones to report.
  topClones <-  unique(x$posidLabel)[1:numClones]
  
  # For each time point, create a data frame for relative abundance plots
  bind_rows(lapply(split(x, x$timePoint), function(x2){
    
    lowAbundData <- dplyr::mutate(x2, posidLabel = 'LowAbund',
                                  totalCells = sum(estAbund),
                                  relAbund   = 20) %>%
      dplyr::slice(1) %>% 
      dplyr::select(patient, cellType, timePoint, posidLabel, totalCells, timePointDays, relAbund)
    
    x3 <- subset(x2, posidLabel %in% topClones)
    if(nrow(x3) == 0) return(lowAbundData)
    x3$totalCells <- sum(x2$estAbund)
    
    lowAbundData$relAbund <- 20 - sum(x3$relAbund)
    bind_rows(lowAbundData,  dplyr::select(x3, patient, cellType, timePoint, posidLabel, totalCells, timePointDays, relAbund))
  }))
}))


x <- abundantClones
x[grepl('LOC101928595,GDPD3\\s\\*~', x$posidLabel),]$posidLabel <- 'chr16+30105436\nGDPD3 *~'


# Create named color vector for unique clones.

o <- subset(x, posidLabel != 'LowAbund')
o <- o[order(o$relAbund, decreasing = TRUE),]

cloneColorsVector <- setNames(c('#eeeeee', grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_distinct(o$posidLabel))),  
                              c('LowAbund', unique(o$posidLabel)))


x$posidLabel <- factor(x$posidLabel, levels = c('LowAbund', unique(o$posidLabel)))
x <- x[order(x$timePointDays),]
x$timePoint  <- factor(x$timePoint, levels = (unique(sort(x$timePoint)))) 

#sub('^.*\\n','',names(cloneColorsVector))

p <- ggplot(x) +
    theme_bw() +
    scale_x_discrete(drop=FALSE) + 
    geom_bar(aes(x=timePoint, y=relAbund/100, fill=posidLabel), stat='identity', color = 'black', size = 0.20) + 
    scale_fill_manual(name = 'Clones', values = cloneColorsVector) +
    #scale_shape_manual(values = c(16, 17, 15), drop = FALSE) +
    labs(x = 'Timepoint', y = 'Relative Sonic Abundance') +
    guides(fill=guide_legend(title.position = "top", ncol=1, keyheight=0.35, default.unit="inch")) +
    scale_y_continuous(labels = scales::percent) + 
    #annotate('text', x=1:length(totalCellLabel), y=1.04, label=totalCellLabel, size=2.7, angle=45, hjust=0.5) +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    guides(fill=guide_legend(ncol=2)) +
    facet_grid(patient~.) 

p

g_label <- sub('^.*\\n','',names(cloneColorsVector))
g_label <- sub('LOC101928298','GDPD3 *~',g_label)
               
ddd <- data.frame(timePoint=c(1.5,1.5,1.5,1.5,1.5),
                  relAbund = c(0.16,0.16,0.16,0.16,0.16),
                  posidLabel=c("LowAbund","LowAbund","LowAbund","LowAbund","LowAbund"),
                 patient = c('Pt 1','Pt 2','Pt 3','Pt 4','Pt 5')
)

ddd$grid_l <- grid_labels 

p_new <- p + scale_fill_manual(
  values = cloneColorsVector,
  limits = names(cloneColorsVector),
  labels = g_label
) + theme(         
  strip.background = element_blank(),
  strip.text.y = element_blank()
) + geom_custom(data = ddd, aes(x=timePoint,y=relAbund,data = grid_l), grob_fun = identity) 



ggsave(p, file = 'patientRelAbund_old.pdf', height = 10, width = 15, units = 'in', useDingbats = FALSE)
ggsave(p_new, file = 'patientRelAbund.pdf', height = 10, width = 15, units = 'in', useDingbats = FALSE)


