########################################
## SH3-dependendy tests on PCA data 
## Figure S4
## author : Pascale Lemieux           ##
## Date : 2023-01-12                  ##
########################################

library(reshape2)
library(tidyverse)
library(magrittr)
library(ggpubr)

firstup <- function(x) {
  substring(x, 2) <- tolower(substring(x, 2))
  x
}

# Importation of the data set for MYO3&MYO5 
MTX2_data <- readRDS('~/ancSH3_paper/SupplementaryMaterial/Data/complete_paralog_data/MTX2data_complete_paralog.rds')

# Use the sh3_sequence nomenclature chosen for the paper
sh3.sequence <-
  factor(
    unique(MTX2_data$strain),
    levels = c(
      "Myo3_Myo5_Duplication_1.MYO3p" ,
      "Myo3_Myo5_Duplication_1.MYO5p" ,
      "Myo3_Myo5_Duplication_3.MYO3p" ,
      "Myo3_Myo5_Duplication_3.MYO5p" ,
      "Myo3_Myo5_post-WGD_1.MYO3p"    ,
      "Myo3_Myo5_post-WGD_1.MYO5p"    ,  
      "Myo3_Myo5_Saccharomycetaceae_1.MYO3p",
      "Myo3_Myo5_Saccharomycetaceae_1.MYO5p",
      "Myo3_Myo5_Saccharomycetales_1.MYO3p",
      "Myo3_Myo5_Saccharomycetales_1.MYO5p",
      "Myo3.MYO3p",
      "Myo3.MYO5p" ,
      "Myo5.MYO3p"          ,
      "Myo5.MYO5p"          ,
      "MYO3p_MYO3SH3.MYO3p" ,
      "MYO3p_MYO5SH3.MYO3p" ,
      "MYO5p_MYO3SH3.MYO5p" ,
      "MYO5p_MYO5SH3.MYO5p" ,
      "MYO3_WT.MYO3p",
      "MYO5_WT.MYO5p"       ,
      "MYO3s.MYO3p"         ,
      "MYO5s.MYO5p"
    ), 
    labels = c('AncA', 'AncA','AncA_3', 'AncA_3','AncB', 'AncB','AncC', 'AncC','AncD', 'AncD',
               'optMyo3', 'swap-optMyo3','swap-optMyo5', 'optMyo5','controlMyo3',
               'swapMyo5', 'controlMyo5','swapMyo3', 'extantMyo3', 'extantMyo5', 'SH3-depleted', 'SH3-depleted')
  )


MTX2_data$sh3_sequence <- 
  factor(
    MTX2_data$strain,
    levels = c(
      "Myo3_Myo5_Duplication_1.MYO3p" ,
      "Myo3_Myo5_Duplication_1.MYO5p" ,
      "Myo3_Myo5_Duplication_3.MYO3p" ,
      "Myo3_Myo5_Duplication_3.MYO5p" ,
      "Myo3_Myo5_post-WGD_1.MYO3p"    ,
      "Myo3_Myo5_post-WGD_1.MYO5p"    ,  
      "Myo3_Myo5_Saccharomycetaceae_1.MYO3p",
      "Myo3_Myo5_Saccharomycetaceae_1.MYO5p",
      "Myo3_Myo5_Saccharomycetales_1.MYO3p",
      "Myo3_Myo5_Saccharomycetales_1.MYO5p",
      "Myo3.MYO3p",
      "Myo3.MYO5p" ,
      "Myo5.MYO3p"          ,
      "Myo5.MYO5p"          ,
      "MYO3p_MYO3SH3.MYO3p" ,
      "MYO3p_MYO5SH3.MYO3p" ,
      "MYO5p_MYO3SH3.MYO5p" ,
      "MYO5p_MYO5SH3.MYO5p" ,
      "MYO3_WT.MYO3p",
      "MYO5_WT.MYO5p"       ,
      "MYO3s.MYO3p"         ,
      "MYO5s.MYO5p"
    ), 
    labels = c('AncA', 'AncA','AncA_3', 'AncA_3','AncB', 'AncB','AncB=C', 'AncC','AncD', 'AncD',
               'optMyo3', 'swap-optMyo3','swap-optMyo5', 'optMyo5','controlMyo3',
               'swapMyo5', 'controlMyo5','swapMyo3', 'extantMyo3', 'extantMyo5', 'SH3-depleted', 'SH3-depleted')
  )

MTX2_data$strain <- str_c(MTX2_data$sh3_sequence, MTX2_data$Bait.Standard_name, sep = '.')

MTX2_data$Bait.Standard_name <- 
  as.factor(firstup(gsub(x = MTX2_data$Bait.Standard_name, pattern = 'p', replacement = '', fixed = T)))

prey <- 
  read.csv('~/ancSH3_paper/SupplementaryMaterial/Data/Prey.Standard_names.csv',
           header = T)[, -1]
xnam <- readRDS('~/ancSH3_paper/Reviews/SupplementaryMaterial/Data/ynam.RDS')

MTX2_data <- 
  merge(MTX2_data, 
        prey, 
        by.x = 'orf', 
        by.y = 'Prey.Systematic_name')

sub_data <- subset(MTX2_data,
                   subset = sh3_sequence %in% c('extantMyo3', 'extantMyo5','SH3-depleted'))

sh3_dep <- 
ggplot(sub_data, aes(x = Prey.Standard_name, y = norm, color = sh3_sequence))+
  facet_grid(rows = vars(Bait.Standard_name))+
  geom_boxplot()+
  scale_color_manual(values = c('#3366CC', 'orangered', 'black'), labels = c('Myo3', 'Myo5', 'SH3-depleted'))+
  stat_compare_means(method = 'wilcox.test', hide.ns = T, label = 'p.signif', 
                     method.args = list(p.adjust.methods = 'BH'), label.y = c(0.97))+
  xlab('Preys')+
  ylab('PPI score')+
  theme_bw()+
  geom_vline(xintercept = 23.5, linetype = 2)+
  geom_vline(xintercept = 14.5, linetype = 2)+
  scale_x_discrete(limits = xnam)+
  scale_y_continuous(limits = c(0, 1))+
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust = 1,size =12),
        panel.grid = element_blank(), legend.position = 'bottom',
        axis.text.y = element_text(size =12),
        legend.text = element_text(size = 14),
        axis.title = element_text(size =14),
        strip.text = element_text(size = 16))+
  guides(color = guide_legend(title = '', nrow = 1,
                              override.aes = aes(label = "")))
#ggsave('~/ancSH3_paper/Reviews/SH3-dep.svg', height = 8)

sub_data <- subset(MTX2_data,
                   subset = sh3_sequence %in% c('extantMyo3', 'extantMyo5'))

para_dep <- 
  ggplot(sub_data, aes(x = Prey.Standard_name, y = norm, color = sh3_sequence))+
  geom_boxplot()+
  scale_color_manual(values = c('#3366CC', 'orangered'), labels = c('Myo3', 'Myo5'))+
  stat_compare_means(hide.ns = T, label = 'p.signif', method.args = list(p.adjust.methods = 'BH'), label.y = c(0.97))+
  xlab('Preys')+
  ylab('PPI score')+
  theme_bw()+
  geom_vline(xintercept = 23.5, linetype = 2)+
  geom_vline(xintercept = 14.5, linetype = 2)+
 scale_x_discrete(limits = xnam)+
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust = 1,size =12),
        panel.grid = element_blank(), legend.position = 'bottom', 
        axis.text.y = element_text(size =12),
        legend.text = element_text(size = 14), 
        axis.title = element_text(size =14),
        strip.text = element_text(size = 16))+
  guides(color = guide_legend(title = '', nrow = 1, 
                              override.aes = aes(label = "")))
#ggsave('~/ancSH3_paper/Reviews/para_dif.png', height = 6)
signif_para <- c('Abp1', 'Arc19', 'Bzz1', 'Cdc10', 'Hek2', 'Las17', 'Lsb3', 'Mid2', 
                 'Rvs167', 'Sfk1', 'Sla1', 'Ssb1', 'Ssb2', 'Syp1', 'Vrp1', 'Ygl242c', 'Yor1', 'Ysc84', 'Cmd1', 'Pkh2')


sub_data <- subset(MTX2_data,
                   subset = (sh3_sequence %in% c('swap-optMyo3', 'optMyo5', 'swap-optMyo5', 'optMyo3')))

sub_data$lab <- 
factor(sub_data$sh3_sequence, 
       levels = c('swap-optMyo3', 'optMyo5', 'swap-optMyo5', 'optMyo3'), 
       labels = c('optMyo3', 'optMyo5', 'optMyo5', 'optMyo3'))

swapM5 <- 
  ggplot(sub_data, aes(x = Prey.Standard_name, y = norm, color = lab))+
    facet_grid(rows = vars(Bait.Standard_name))+
  geom_boxplot()+
  scale_color_manual(values = c('#3366CC','orangered'),
                     labels = c('optMyo3', 'optMyo5'))+
  stat_compare_means(method = 'wilcox.test', label.y = 0.97,
    hide.ns = T, label = 'p.signif', method.args = list(p.adjust.methods = 'BH'))+
  #facet_grid(cols = vars(Bait.Standard_name))+
  geom_vline(xintercept = 23.5, linetype = 2)+
  geom_vline(xintercept = 14.5, linetype = 2)+
  scale_x_discrete(limits = xnam)+
  xlab('Preys')+
  ylab('PPI score')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust = 1,size =12),
        panel.grid = element_blank(), legend.position = 'bottom', 
        axis.text.y = element_text(size =12),
        legend.text = element_text(size = 14), 
        axis.title = element_text(size =14),
        strip.text = element_text(size = 16))+
  guides(color = guide_legend(title = '', nrow = 1, 
                              override.aes = aes(label = "")))

sub_data <- subset(MTX2_data,
                   subset = (sh3_sequence %in% c('AncA', 'optMyo3', 'optMyo5')))

swapAM3 <- 
  ggplot(sub_data, aes(x = Prey.Standard_name, y = norm, color = sh3_sequence))+
  geom_boxplot()+
  scale_color_manual(values = c('#9bcd9bff', '#3366CC', 'orangered'),
                     labels = c('AncA', 'optMyo3', 'optMyo5'))+
  stat_compare_means(method = 'wilcox.test', label.y = 0.97,
                     hide.ns = T, label = 'p.signif', method.args = list(p.adjust.methods = 'BH'))+
  facet_grid(rows = vars(Bait.Standard_name))+
  geom_vline(xintercept = 23.5, linetype = 2)+
  geom_vline(xintercept = 14.5, linetype = 2)+
  scale_x_discrete(limits = xnam)+
  xlab('Preys')+
  ylab('PPI score')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust = 1,size =12),
        panel.grid = element_blank(), legend.position = 'bottom', 
        axis.text.y = element_text(size =12),
        legend.text = element_text(size = 14), 
        axis.title = element_text(size =14),
        strip.text = element_text(size = 16))+
  guides(color = guide_legend(title = '', nrow = 1, 
                              override.aes = aes(label = "")))



library(cowplot)
plot_grid(para_dep+theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),axis.text.x = element_blank(), plot.margin = unit(c(0.2, 0, 0, 0), "cm")),
          sh3_dep+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title.x = element_blank()),
          #swapM3+theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")),
          swapM5+theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),axis.text.x = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")),
          swapAM3+theme(legend.position = 'bottom', plot.margin = unit(c(0, 0, 0, 0.2), "cm")),
          #swapAM5+theme(legend.position = 'bottom', plot.margin = unit(c(0, 0, 0, 0), "cm")), 
          nrow = 4, labels = 'AUTO', align = 'v', axis = 'rl', rel_heights = c(1,1.6,1.6,1.7),
          label_size = 16, label_fontface = 'plain') #rel_heights = c(1,0.8,0.8), greedy = F)

ggsave('~/ancSH3_paper/Reviews/SupplementaryMaterial/FigureSX_RD.svg', width = 12, height = 20)


signif_para <- c('Abp1', 'Arc19', 'Bzz1', 'Cdc10', 'Hek2', 'Las17', 'Lsb3', 'Mid2', 
                 'Rvs167', 'Sfk1', 'Sla1', 'Ssb1', 'Ssb2', 'Syp1', 'Vrp1', 'Ygl242c', 'Yor1', 'Ysc84', 'Cmd1', 'Pkh2')

saveRDS(signif_para, '~/ancSH3_paper/Reviews/SupplementaryMaterial/Data/signif_para.rds')

# Function does a krusal-wallis test, if significant it computes  
# both hypothesis (greater and less) pairwise wilcoxon test with 
# Benjamini & Hochberg correction. The more significant hypothesis is 
# considered and its p-value is conserved. Return a data frame with
# the adjusted p-value for each comparison.
wilcoxorf <-
  function(sh3.sequence, subdata) {
    f <- vector(mode = 'list', length = length(unique(subdata$orf)))
    names(f) <- unique(subdata$orf)
    
    for (k in unique(subdata$orf)) {
      subsubdata <- subset(subdata,
                           subset = orf == k)
      
      ktest <- kruskal.test(norm ~ strain, data = subsubdata)
      
      if (ktest$p.value < 0.05) {
        l <-
          pairwise.wilcox.test(
            subsubdata$norm,
            subsubdata$strain,
            p.adjust.method = "BH",
            alternative = 'less'
          ) 
        
        g <-
          pairwise.wilcox.test(
            subsubdata$norm,
            subsubdata$strain,
            p.adjust.method = "BH",
            alternative = 'greater'
          )
        
      l <- as.data.frame(l$p.value)
      l <- cbind('strain' = rownames(l), l)
      
      l <- reshape(
        l,
        direction = 'long',
        varying = list(colnames(l[,-1])),
        v.names = 'p.value.adj',
        times = colnames(l[,-1]),
        timevar = 'strain2'
      )
      l <- na.omit(l)
      
      l$orf <- rep(paste(k), length.out = nrow(l))
      l$comp <- 'less'
      l$id <- NULL
      
      
      g <- as.data.frame(g$p.value)
      g <- cbind('strain' = rownames(g), g)
      
      g <- reshape(
        g,
        direction = 'long',
        varying = list(colnames(g[,-1])),
        v.names = 'p.value.adj',
        times = colnames(g[,-1]),
        timevar = 'strain2'
      )
      g <- na.omit(g)
      
      g$orf <- rep(paste(k), length.out = nrow(g))
      g$comp <- 'greater'
      g$id <- NULL
      
      x <- merge(l,
                 g,
                 by = c('strain', 'strain2', 'orf'))
      
      x$p.value.adjusted <- ''
      x$comp <- ''
      
      for (i in 1:nrow(x)) {
        a.x <- abs(x[i, 'p.value.adj.x'])
        a.y <- abs(x[i, 'p.value.adj.y'])
        
        if (a.x > a.y) {
          x[i, c('p.value.adjusted', 'comp')] <-
            x[i, c('p.value.adj.y', 'comp.y')]
          
        } else if (a.x < a.y) {
          x[i, c('p.value.adjusted', 'comp')] <-
            x[i, c('p.value.adj.x', 'comp.x')]
        }
      }
      
      
      f[[k]] <- x[,-c(4:7)]
      
    }
  }  
    
    x <- unlist(lapply(f, FUN = is.null))
    f <- f[!x]
    
    f <- do.call(rbind.data.frame, f)
    
    colnames(f) <-
      c('strain1',
        'strain2',
        'DHFR3',
        'p.value.adjusted',
        'str1.vs.str2')
    rownames(f) <- NULL
    
    
    liste <- strsplit(f$strain1, split = '.', fixed = T)
    
    x <-
      data.frame(matrix(unlist(liste), nrow = length(liste), byrow = TRUE))
    colnames(x) <- c('sh3.sequence1', 'para1')
    
    liste <- strsplit(f$strain2, split = '.', fixed = T)
    y <-
      data.frame(matrix(unlist(liste), nrow = length(liste), byrow = TRUE))
    colnames(y) <- c('sh3.sequence2', 'para2')
    
    f <- cbind(f, x, y)
    
    
    return(f)
  }

padj.MYO3.MYO5 <- wilcoxorf(sh3.sequence, MTX2_data)

# Data transformation
padj.MYO3.MYO5$p.value.adjusted <- 
  as.numeric(padj.MYO3.MYO5$p.value.adjusted)
padj.MYO3.MYO5 %<>%
  mutate('log10.pvalue.adjusted' = log10(p.value.adjusted))
padj.MYO3.MYO5$`log10.pvalue.m.adjusted` <- 0


# Transform each log10(p-value adjusted) according to its tag (less or greater)
# Values tagged with 'less' are kept negative and values tagged with 'greater' are
# transformed into positive values 
for (i in 1:nrow(padj.MYO3.MYO5)) {
  
  if(padj.MYO3.MYO5[i, 'str1.vs.str2'] == 'greater'){
    padj.MYO3.MYO5[i, 'log10.pvalue.m.adjusted'] <- padj.MYO3.MYO5[i, 'log10.pvalue.adjusted']
  } else if(padj.MYO3.MYO5[i, 'str1.vs.str2'] == 'less'){
    x <- padj.MYO3.MYO5[i, 'log10.pvalue.adjusted']
    padj.MYO3.MYO5[i, 'log10.pvalue.m.adjusted'] <- -(as.numeric(x))
  }
  
}


# Keep only comparison within a paralog, comparing its variants
psig.MYO3 <- subset(padj.MYO3.MYO5, 
                    subset = ( para1 == 'MYO3p' & para2 == 'MYO3p'))

psig.MYO5 <- subset(padj.MYO3.MYO5, 
                    subset = ( para1 == 'MYO5p' & para2 == 'MYO5p'))


# Keep only comparisons in which one of the variants compared 
# is the codon optimized extant SH3 (optSH3)
psig.MYO3 <- subset(psig.MYO3, 
                    subset = (sh3.sequence1 == 'optMyo3' | sh3.sequence2 == 'optMyo3'))

psig.MYO5 <- subset(psig.MYO5, 
                    subset = (sh3.sequence1 == 'optMyo5' | sh3.sequence2 == 'optMyo5'))


# Reorganizing the data set for the sh3.sequence1 = the reference SH3 variant
# This function assigns to new columns the reorganized sh3 variants and the 
# log10.pvalue,m.adjusted correct sign for the new comparison order
org <- 
  function(data, ref){
    data$sh3.1 <- ''
    data$sh3.2 <- ''
    for (i in 1:nrow(data)) {
      
      if(data[i, 'sh3.sequence1'] != ref ){
        data[i, 'sh3.1'] <- data[i, 'sh3.sequence2']
        data[i, 'sh3.2'] <- data[i, 'sh3.sequence1']
        data[i, 'log10.pvalue.m.adjusted'] <- -(data[i, 'log10.pvalue.m.adjusted'])
        
      } else if (data[i, 'sh3.sequence1'] == ref){
        data[i, 'sh3.1'] <- data[i, 'sh3.sequence1']
        data[i, 'sh3.2'] <- data[i, 'sh3.sequence2']
      }
    }
    return(data)
  }

# 
psig.MYO3 <- 
  org(psig.MYO3, 'optMyo3')

psig.MYO5 <- 
  org(psig.MYO5, 'optMyo5')

psig.MYO3.MYO5 <- 
  rbind(psig.MYO3, psig.MYO5)

# Select significant tests
subpadj <- subset(psig.MYO3.MYO5, 
                  subset = p.value.adjusted < 0.05)


# Determination of SH3 dependent interactions (gain and loss)
# Creation of a list containing 4 categories of interaction partners according
# to the p-value obtained with the pairwise wilcoxon test between SH3 WT, SH3 optimized and
# SH3 stuffed strains

submyo3 <- 
   subset(padj.MYO3.MYO5,
                  subset = (para1 == 'MYO3p' & para2 == 'MYO3p' &
                              sh3.sequence1 %in% c( 'optMyo3', 'SH3-depleted') &
                              sh3.sequence2 %in% c( 'optMyo3', 'SH3-depleted')))


submyo5 <- subset(padj.MYO3.MYO5,
                  subset = (para1 == 'MYO5p' & para2 == 'MYO5p' &
                              sh3.sequence1 %in% c( 'optMyo5', 'SH3-depleted') &
                                 sh3.sequence2 %in% c( 'optMyo5', 'SH3-depleted')))


# SH3-dependency attribution

submyo3[submyo3$p.value.adjusted < 0.05 & submyo3$str1.vs.str2 == 'less', "SH3.dep" ] <- T
submyo5[submyo5$p.value.adjusted < 0.05 & submyo5$str1.vs.str2 == 'less', "SH3.dep" ] <- T

dependency <-
  rbind(submyo3, submyo5)
dependency$DHFR12.strain <-
  matrix(unlist(strsplit(dependency$strain1, '.', fixed = T)), ncol = 2, 
         byrow = T, nrow = nrow(dependency))[, 2]


# Merge SH3-dependency with data

MTX2_data <- 
merge(MTX2_data, 
      dependency[, c(3,12,13)], 
      by.x = c('orf', 'Bait.Standard_name'), 
      by.y = c('DHFR3', 'DHFR12.strain'), 
      all = T)

# Merge opt SH3 difference

subpadj$optSH3_dif <- T

MTX2_data <- 
merge(MTX2_data, 
      subpadj[, c(3,9,13,14)], 
      by.x = c('orf','Bait.Standard_name', 'sh3_sequence'),
      by.y = c('DHFR3', 'para2','sh3.2'), 
      all.x = T)


# replace all NAs (genereated because no differences are observed
# between baits-preys pair, kruskall-wallis) by FALSE (no SH3 dependency)
library(gtools)
MTX2_data$SH3.dep <- 
  na.replace(MTX2_data$SH3.dep, replace = FALSE)

MTX2_data$optSH3_dif <- 
  na.replace(MTX2_data$optSH3_dif, replace = FALSE)

# Export data (used for SupplementaryData.R)
saveRDS(file = '~/ancSH3_paper/SupplementaryMaterial/Data/complete_paralog_data/MTXdata_sh3dependency.rds', MTX2_data)
