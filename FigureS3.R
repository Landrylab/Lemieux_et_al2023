######################################
## Supp Figure 2                    ##
## author : Pascale Lemieux         ##
## Date : 2023-01-12                ##
######################################
library(ggplot2)
library(tidyverse)
library(stringr)
library(ggpubr)
library(plyr)
library(ggcorrplot)
library(ggrepel)
library(magrittr)
library(viridis)
library(cowplot)

firstup <- function(x) {
  substring(x, 2) <- tolower(substring(x, 2))
  x
}
# colony area distribution, from PCA_analysis_complete_paralog.R
p1 <- readRDS(file = '~/ancSH3_paper/Reviews/FigurePanels/SuppFig2A.rds')

# look for expression biais in raw data, from PCA_analysis_complete_paralog.R
p2 <- readRDS(file = '~/ancSH3_paper/Reviews/FigurePanels/SuppFig2b.rds')+
        theme(legend.position = 'none')

# technical replicates comparison, from PCA_analysis_complete_paralog.R
p3 <- readRDS(file = '~/ancSH3_paper/Reviews/FigurePanels/SuppFig2C.rds')+
      xlab('PPI score replicate 1')+
      ylab('PPI score replicate 2')

# Cytometry expression data

SuppFigG <- readRDS(file = '~/ancSH3_paper/Reviews/FigurePanels/SuppFig2G.rds')


# Read Table S4 with liquid validation data
gc_out <- 
  read_xlsx('~/ancSH3_paper/Reviews/SupplementaryMaterial/TableS4.xlsx')

# compute median PPI score for liquid PCA

gc_out %>%
group_by(Prey.Systematic_name, sh3_sequence, Bait.Standard_name) %>%
  dplyr::mutate(med.PPI_score_gc = median(PPI_score_gc, na.rm = T)) -> gc_out

# Read TableS1 : complete paralog PCA result to compare 
# with liquid PCA result
solid_MTX <- 
  read_xlsx('~/ancSH3_paper/Reviews/SupplementaryMaterial/TableS1.xlsx')

solid_MTX$sh3_sequence <- 
  as.character(solid_MTX$sh3_sequence)

solid_MTX[(grepl(solid_MTX$sh3_sequence, pattern = 'swap-optMyo3')), "sh3_sequence"] <- 'optMyo3'
solid_MTX[(grepl(solid_MTX$sh3_sequence, pattern = 'swap-optMyo5')), "sh3_sequence"] <- 'optMyo5'

comp_MTX <- 
  merge(gc_out, 
        unique(solid_MTX[, c(1:4,16)]), 
        by = c('Prey.Systematic_name', 'sh3_sequence', 'Bait.Standard_name'),
        all= F, 
        sort = F)
# Correlation solid vs liquid PCA results
# Remove non inoculated well data
comp_MTX <- comp_MTX[!(comp_MTX$med.PPI_score>0.75 & comp_MTX$med.PPI_score_gc < 10),]

# keep only interesting columns
comp_MTX <- 
  unique(comp_MTX[, c(1:3, 16,18)])

comp_MTX <- na.omit(comp_MTX)

# Labelled each Bait strain as an experiment condition
# swap.optSH3

comp_MTX$exp <- ''

#comp_MTX[(grepl(comp_MTX$sh3_sequence, pattern = 'swap-opt')), "exp"] <- 'swap.optSH3'

comp_MTX[(comp_MTX$sh3_sequence =='optMyo3' & comp_MTX$Bait.Standard_name == 'Myo5'), "exp"] <- 'swap.optSH3'
comp_MTX[(comp_MTX$sh3_sequence =='optMyo5' & comp_MTX$Bait.Standard_name == 'Myo3'), "exp"] <- 'swap.optSH3'

# AncD
comp_MTX[(comp_MTX$sh3_sequence =='AncD'), "exp"] <- 'AncD'

#opt.SH3
comp_MTX[(comp_MTX$sh3_sequence =='optMyo5' & comp_MTX$Bait.Standard_name == 'Myo5'), "exp"] <- 'optSH3'
comp_MTX[(comp_MTX$sh3_sequence =='optMyo3' & comp_MTX$Bait.Standard_name == 'Myo3'), "exp"] <- 'optSH3'

comp_MTX[(comp_MTX$sh3_sequence =='SH3-depleted'), "exp"] <- 'SH3-depleted'


# SuppFig2G
p7 <- 
ggplot(unique(comp_MTX))+
  facet_grid(cols = vars(Bait.Standard_name))+
  geom_point(aes(med.PPI_score_gc, med.PPI_score, color= exp), size = 2.5, alpha = 0.8)+
  stat_cor(aes(med.PPI_score_gc, med.PPI_score, color = exp), method = 'spearman', size = 4.5, cor.coef.name = c('r'), 
           label.y = c(0.05, 0.15, 0.25, 0.35), label.x = c(26))+
  stat_cor(aes(med.PPI_score_gc, med.PPI_score), method = 'spearman', size = 4.5, cor.coef.name = c('r'))+
  scale_color_manual(values = c('#551A8B', 'darkcyan','grey30', '#A20056'))+
  theme_bw()+
  ylim(0,1)+
  xlim(0,40)+
  ylab('med. solid PPI score')+
  xlab('med. liquid PPI score')+
  theme(legend.position = 'bottom', legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        #legend.title = element_text(size = 16),
        axis.title = element_text(size =14), 
        axis.text = element_text(size =12),
        strip.text.x = element_text(size = 16), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(
         color = guide_legend(nrow = 1, title = "",
                              override.aes = aes(label = "")))


# FigSupp2H : codon optimization effect

#p6_l <- 
#  readRDS('~/ancSH3_paper/SupplementaryMaterial/FigurePanels/SuppFig2H_l.rds')
#p6_r <- 
#  readRDS('~/ancSH3_paper/SupplementaryMaterial/FigurePanels/SuppFig2H_r.rds')
  
# SuppFig2D : Result comparison with BioGrid reference dataset

library(ggupset)
# Read file with the comparison information
#exp_collapsed <- 
#  readRDS('~/ancSH3_paper/SupplementaryMaterial/Data/BioGridcomp.rds')

#exp_collapsed$DHFR12.strain <- 
#  firstup(exp_collapsed$DHFR12.strain)
# Read the file with the preys used as abundance controls
abondance_c <-
  c(
    "YBL032W",
    "YCR002C",
    "YCR088W"   ,
    "YER165W"   ,
    "YER177W"   ,
    "YGL234W" ,
    "YGR162W",
    "YGR192C",
    "YGR240C",
    "YHL034C" ,
    "YHR064C ",
    "YLL026W" ,
    "YNL209W",
    "YPL240C"
  )
  
# SuppFig2D & E :Test for expression biais after normalization & SH3-depletion impact on PPI scores
##here

subdata <- 
unique(
  subset(solid_MTX, 
       subset = sh3_sequence %in% c('SH3-depleted', 'extantMyo3', 'extantMyo5'), 
       select = c('Prey.Systematic_name', 'Bait.Standard_name', 'sh3_sequence', 'med.PPI_score', 'SH3_dep')))

subdata$sh3_sequence <- 
factor(subdata$sh3_sequence, 
       levels = c('extantMyo3', 'extantMyo5', 'SH3-depleted'), 
       labels = c('extantSH3', 'extantSH3','SH3-depleted'))


subdata$abond_c <- 
subdata$Prey.Systematic_name %in% abondance_c


# Test to validate there is no expression biais betweem the paralog after normalization
pval <- 
  compare_means(med.PPI_score ~ Bait.Standard_name, subdata, method = 'wilcox.test', p.adjust.method = 'BH')

# SuppFig2D
p3.1 <- 
ggplot(subdata[subdata$abond_c & subdata$sh3_sequence == 'extantSH3', ])+
  geom_boxplot(aes(x = as.factor(Bait.Standard_name), 
                  y =med.PPI_score), width = 0.3, alpha = 0.3, color = 'grey30', fill = 'grey65', outlier.color = 'transparent')+
  geom_jitter(aes(x = as.factor(Bait.Standard_name), 
                 y =med.PPI_score), 
              width= 0.1)+
  stat_compare_means(aes(x = as.factor(Bait.Standard_name), 
                         y =med.PPI_score) ,size = 4.5, label.x = 1.2, label.y = 0.98)+
  #facet_grid(cols = vars(recovered), scales = 'free', drop = T)+
  theme_bw() +
  #xlim(0,1)+  
  ylim(0,1)+
  ylab('med. PPI score')+
  xlab('WT bait')+
   scale_color_manual(values = c('grey30'), 
                     labels = c(''))+
  theme( legend.position = 'none',
         legend.text = element_text(size = 14),
         axis.title = element_text(size = 14), 
         axis.text = element_text(size =12), 
         legend.title = element_text(size =16),
         strip.text.x = element_text(size = 16), 
         panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(title = "",
                              override.aes = aes(label = ""), nrow = 1), 
         fill =guide_legend(title = '', nrow = 1))  


## Supplementary Figure 2 assembly

top <- 
plot_grid(p1, p2, 
          #p3+theme(axis.title.y = element_text()),
          p3.1+theme(legend.position = 'none')+guides(color = guide_legend(nrow =2, title = ''), 
                            fill = 'none'), 
          ncol =3, align = 'h', axis = 'b', 
          labels = c('A', 'B', 'C'), label_fontface = 'plain', 
          rel_widths = c(1.2,1,1))

middle <- 
  plot_grid(p3,
            SuppFigG+theme(axis.title.x = element_text(vjust = -2)), 
            align = 'h', axis = 'tb',
            ncol =2, rel_widths = c(0.5,0.9),
            labels = c('D', 'E'), label_fontface = 'plain')

 bottom <- 
   plot_grid(  p7,
             nrow = 1,ncol = 1, 
             labels = c('F'), label_fontface = 'plain')


plot_grid(top, 
          middle,
          bottom, 
          nrow = 3, labels = '', rel_heights = c(1,1.1,1), 
          align = 'v', axis = 'rl')
library(svglite)
ggsave('~/ancSH3_paper/Reviews/SupplementaryMaterial/FigureS2.svg', 
       height = 13, width = 12)

