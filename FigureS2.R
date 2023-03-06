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
# colony area distribution
p1 <- readRDS(file = '~/ancSH3_paper/SupplementaryMaterial/FigurePanels/SuppFig2A.rds')

# look for expression biais in raw data
p2 <- readRDS(file = '~/ancSH3_paper/SupplementaryMaterial/FigurePanels/SuppFig2b.rds')+
        theme(legend.position = 'none')

# technical replicates comparison
p3 <- readRDS(file = '~/ancSH3_paper/SupplementaryMaterial/FigurePanels/SuppFig2C.rds')+
      xlab('PPI score replicate 1')+
      ylab('PPI score replicate 2')

# Cytometry expression data

SuppFigG <- readRDS(file = '~/ancSH3_paper/SupplementaryMaterial/FigurePanels/SuppFig2G.rds')


# Read Table S4 with liquid validation data
gc_out <- 
  read.csv('~/ancSH3_paper/SupplementaryMaterial/TableS4.csv')[, -1]

# compute median PPI score for liquid PCA

gc_out %<>%
  dplyr::group_by(Prey.Systematic_name, sh3_sequence, Bait.Standard_name)%>%
  dplyr::mutate(n(), med.PPI_score_gc = median(PPI_score_gc))


# Read TableS1 : complete paralog PCA result to compare 
# with liquid PCA result
solid_MTX <- 
  read.csv(file = '~/ancSH3_paper/SupplementaryMaterial/TableS1.csv')[, -1]

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
  unique(comp_MTX[, c(1:3, 15:17)])

comp_MTX <- na.omit(comp_MTX)

# Labelled each Bait strain as an experiment condition
# swap.optSH3

comp_MTX$exp <- ''

#comp_MTX[(grepl(comp_MTX$sh3_sequence, pattern = 'swap-opt')), "exp"] <- 'swap.optSH3'

comp_MTX[(comp_MTX$sh3_sequence =='optMyo3' & comp_MTX$Bait.Standard_name == 'Myo5'), "exp"] <- 'swap.optSH3'
comp_MTX[(comp_MTX$sh3_sequence =='optMyo5' & comp_MTX$Bait.Standard_name == 'Myo3'), "exp"] <- 'swap.optSH3'

# Ancc
comp_MTX[(comp_MTX$sh3_sequence =='AncC'), "exp"] <- 'AncC'

#opt.SH3
comp_MTX[(comp_MTX$sh3_sequence =='optMyo5' & comp_MTX$Bait.Standard_name == 'Myo5'), "exp"] <- 'optSH3'
comp_MTX[(comp_MTX$sh3_sequence =='optMyo3' & comp_MTX$Bait.Standard_name == 'Myo3'), "exp"] <- 'optSH3'

comp_MTX[(comp_MTX$sh3_sequence =='SH3-depleted'), "exp"] <- 'SH3-depleted'


# SuppFig2G
p7 <- 
ggplot(unique(comp_MTX))+
  geom_point(aes(med.PPI_score_gc, med.PPI_score, color= exp, shape= Bait.Standard_name), size = 2.5, alpha = 0.8)+
  stat_cor(aes(med.PPI_score_gc, med.PPI_score, color = exp), method = 'spearman', size = 4.5, cor.coef.name = c('r'), 
           label.y = c(0.05, 0.15, 0.25, 0.35), label.x = c(28))+
  stat_cor(aes(med.PPI_score_gc, med.PPI_score), method = 'spearman', size = 4.5, cor.coef.name = c('r'))+
  scale_color_manual(values = c('#551A8B', '#3366CC','grey20', 'orangered'))+
  theme_bw()+
  ylim(0,1)+
  ylab('med. PPI score')+
  xlab('med. liquid PPI score')+
  theme(legend.position = 'bottom', legend.title = element_blank(), 
        legend.text = element_text(size = 14),
        #legend.title = element_text(size = 16),
        axis.title = element_text(size =14), 
        axis.text = element_text(size =12),
        strip.text.x = element_text(size = 16), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(shape = guide_legend(nrow = 2), 
         color = guide_legend(nrow = 2, title = "",
                              override.aes = aes(label = "")))


# FigSupp2H : codon optimization effect

#p6_l <- 
#  readRDS('~/ancSH3_paper/SupplementaryMaterial/FigurePanels/SuppFig2H_l.rds')
#p6_r <- 
#  readRDS('~/ancSH3_paper/SupplementaryMaterial/FigurePanels/SuppFig2H_r.rds')
  
# SuppFig2D : Result comparison with BioGrid reference dataset

library(ggupset)
# Read file with the comparison information
exp_collapsed <- 
  readRDS('~/ancSH3_paper/SupplementaryMaterial/Data/BioGridcomp.rds')

exp_collapsed$DHFR12.strain <- 
  firstup(exp_collapsed$DHFR12.strain)
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

# Identify the comparison qith the abundance control information
exp_collapsed[exp_collapsed$orf %in% abondance_c, "recovered"] <-'Abundance control' 
exp_collapsed$DHFR12.strain <- 
  firstup(exp_collapsed$DHFR12.strain)

# Summary table with detected, undetected PPI vs the reference database
table(exp_collapsed$recovered)

# Remove NAs
exp_collapsed <- 
  exp_collapsed[!is.na(exp_collapsed$orf), ]

# proportion of detected PPI compared to the reference database
120/(89+120)

# Verify number of new interaction excluding abundance control
exp_collapsed[exp_collapsed$method_list == 'Not reported in BioGrid', "method_list"] <- 'Not reported in BioGRID'

# Count PPI previously detected by two methods that were detected in the PCA experiment 
length(exp_collapsed[lapply(exp_collapsed$method_list, 'length') >=2, "recovered"] == 'Detected')
#78/95

# Known interactors included in the array
length(unique(exp_collapsed[!exp_collapsed$method_list == 'Not reported in BioGRID', 'orf']))
# 93

p4 <- 
  ggplot(exp_collapsed, aes(method_list, fill = recovered))+
  facet_grid(row = vars(DHFR12.strain))+
  geom_bar()+
  theme_bw()+
  scale_x_upset(n_intersections = 20)+
  ylab('Reported in BioGRID')+
  xlab(NULL)+
  scale_fill_manual(values = c('grey20','darkcyan', 'grey65'))+
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(), legend.position = 'bottom', 
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size =12), 
        strip.text = element_text(size =16), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(fill = guide_legend(nrow = 2))+
  theme_combmatrix(combmatrix.panel.line.size = 1, 
                  combmatrix.panel.point.size = 2, 
                  combmatrix.label.text = element_text(size = 12), 
                  combmatrix.panel.point.color.fill = 'grey20', 
                  combmatrix.panel.line.color = 'grey20')


  
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

subdata <- 
merge(subdata, 
      exp_collapsed[, c(1:3)], 
      by.x = c('Prey.Systematic_name', 'Bait.Standard_name'), 
      by.y = c('orf', 'DHFR12.strain'))

subdata$recovered <- factor(subdata$recovered, 
                            levels = c('Abundance control', 'Detected'),
                            labels = c('Abundance control', 'Specific PPIs'))


# Test to validate there is no expression biais betweem the paralog after normalization
pval <- 
  compare_means(med.PPI_score ~ Bait.Standard_name, subdata, method = 'wilcox.test', p.adjust.method = 'BH')

# SuppFig2D
p3.1 <- 
ggplot(subdata[subdata$recovered == 'Abundance control' & subdata$sh3_sequence == 'extantSH3', ])+
  geom_boxplot(aes(x = as.factor(Bait.Standard_name), 
                  y =med.PPI_score), width = 0.3, alpha = 0.3, color = 'grey30', fill = 'grey65', outlier.color = 'transparent')+
  geom_jitter(aes(x = as.factor(Bait.Standard_name), 
                 y =med.PPI_score, color = recovered), 
              width= 0.1)+
  stat_compare_means(aes(x = as.factor(Bait.Standard_name), 
                         y =med.PPI_score) ,size = 4.5, label.x = 1.2, label.y = 0.98)+
  facet_grid(cols = vars(recovered), scales = 'free', drop = T)+
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


#p3.2 <- 
  # ggplot(subdata[subdata$recovered == 'Specific PPIs', ])+
  # geom_boxplot(aes(x = as.factor(sh3_sequence), 
  #                  y =PPI_score), width = 0.3, alpha = 0.3, color = 'grey30', fill = 'grey65', outlier.color = 'transparent')+
  # geom_jitter(aes(x = as.factor(sh3_sequence), 
  #                 y =PPI_score, color = SH3_dep), 
  #             width= 0.1)+
  # stat_compare_means(aes(x = as.factor(sh3_sequence), 
  #                        y =PPI_score) ,size = 4.5, label.x = 1, label.y = 0.98)+
  # facet_grid(cols = vars(recovered), scales = 'free', drop = T)+
  # theme_bw() +
  # #xlim(0,1)+  
  # ylim(0,1)+
  # ylab('PPI score')+
  # xlab('paralog variant')+
  # scale_color_manual(values = c('grey65', 'darkcyan'), 
  #                    labels = c('SH3-independent', 'SH3-dependent'))+
  # scale_fill_manual(values = c('grey65', 'darkcyan'), 
  #                   labels = c('SH3-independent', 'SH3-dependent'))+
  # theme( legend.position = 'bottom',
  #        legend.text = element_text(size = 14),
  #        axis.title = element_text(size = 14), 
  #        axis.text = element_text(size =12), 
  #        legend.title = element_text(size =16),
  #        strip.text.x = element_text(size = 16), 
  #        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  # guides(color = guide_legend(title = "",
  #                             override.aes = aes(label = ""), nrow = 2), 
  #        fill =guide_legend(title = '', nrow = 1))  

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
            p4+theme(axis.title.y = element_text(hjust = 2.2)), 
            #SuppFigG+theme(legend.text = element_text(size = 14), 
            #               axis.title.y = element_text(hjust = 4)),
            ncol =2, rel_widths = c(0.5,0.9),
            labels = c('D', 'E'), label_fontface = 'plain')

bottom <- 
  plot_grid(SuppFigG+theme(legend.text = element_text(size = 14), 
                          axis.title.y = element_text(hjust = 4)),
            p7, rel_widths = c(1, 1),
            nrow = 1, ncol = 2,  margin = 't', align = 'h', axis = 'b', 
            labels = c('F','G'), label_fontface = 'plain')


plot_grid(top, 
          middle,
          bottom, 
          nrow = 3, labels = '', rel_heights = c(1,1.1,1.2), 
          align = 'v', axis = 'rl')

ggsave('~/ancSH3_paper/SupplementaryMaterial/FigureS2.png', 
       height = 13, width = 14)

