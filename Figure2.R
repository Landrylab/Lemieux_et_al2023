######################################
## Figure 2                         ##
## author : Pascale Lemieux         ##
## Date : 2023-06-02                ##
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
library(magick)
library(cowplot)

# Import data from TableS1 
data <-
  readRDS('~/ancSH3_paper/Reviews/SupplementaryMaterial/SupplementaryData1/SupplementaryData1_complete_paralog.rds')


# verification of the optimized sequence effect on PPI scores
wt <- subset(data[, -c(5:8, 9:15)],
             subset = sh3_sequence %in% c('extantMyo3', 'extantMyo5'))
opt <- subset(data[, -c(5:8, 9:15)],
              subset = sh3_sequence %in% c('optMyo3', 'optMyo5'))

datat <- rbind(wt, opt)

seq_opt <- 
unique(merge(wt, 
      opt, 
      by = c('Bait.Standard_name', 
             'Prey.Standard_name', 'SH3_dep')))


seq_opt$sh3_sequence.y <-   
  matrix(unlist(strsplit(as.character(seq_opt$sh3_sequence.y), split = 't', fixed = T)), nrow = nrow(seq_opt), ncol = 2, byrow = T)[, 2]

seq_opt <- 
  seq_opt[seq_opt$Bait.Standard_name == seq_opt$sh3_sequence.y, ]

m3 <- 
  seq_opt[seq_opt$Bait.Standard_name == 'Myo3',]
m5 <- 
  seq_opt[seq_opt$Bait.Standard_name == 'Myo5',]

compare_means(med.PPI_score ~ sh3_sequence, data = unique(datat[datat$Bait.Standard_name =='Myo3', ]), paired = T, 
              method = 'wilcox.test')
# .y.       group1  group2          p p.adj p.format p.signif method  
# <chr>     <chr>   <chr>       <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 PPI_score optMyo3 extantMyo3 0.0431 0.043 0.043    *        Wilcoxon
compare_means(med.PPI_score ~ sh3_sequence, data = unique(datat[datat$Bait.Standard_name =='Myo5', ]), paired = T, 
              method = 'wilcox.test')
#.y.       group1  group2              p     p.adj p.format p.signif method  
#<chr>     <chr>   <chr>           <dbl>     <dbl> <chr>    <chr>    <chr>   
#  1 PPI_score optMyo5 extantMyo5 0.00000353 0.0000035 3.5e-06  ****     Wilcoxon

# optimize sequence induces a change in PPI score, we have to use optSH3 as control for the ancestral state

# Figure 2B : Comparison of PPI score between paralogs


data_para <- subset(data,
                    subset = sh3_sequence %in% c('extantMyo3', 'extantMyo5'))

data_para$Low_rep <- data_para$n.bio_rep <=2

data_para <- 
  pivot_wider(unique(data_para[,c(1,2,4, 16:19)]), 
              values_from = c(med.PPI_score, SH3_dep, optSH3_dif, Low_rep),
              names_from = sh3_sequence)

signif_para <- 
  readRDS('~/ancSH3_paper/Reviews/SupplementaryMaterial/Data/signif_para.rds')


Fig2B <- 
  ggplot(data_para, 
         aes(x = med.PPI_score_extantMyo5, y = med.PPI_score_extantMyo3, 
             color = as.factor(SH3_dep_extantMyo3 | SH3_dep_extantMyo5)))+
  geom_point(size = 4, aes(shape = Prey.Standard_name %in% signif_para))+
  #geom_point(data = data_para[data_para$Low_Rep_extantMyo3 | data_para$Low_Rep_extantMyo5, ], 
  #           mapping = aes(med.PPI_score_extantMyo5, med.PPI_score_extantMyo3), 
  #           color = 'black', size = 4, shape = 18)+
  theme_bw()+
  xlim(0,1)+
  ylim(0,1)+
  stat_cor( aes(color = as.factor(SH3_dep_extantMyo3 | SH3_dep_extantMyo5)), 
    size = 4.5, method = 'spearman', label.y = c( 0.95, 0.89),
    cor.coef.name = c('r'))+
  ylab('med. Myo3 PPI score')+
  xlab('med. Myo5 PPI score')+
  geom_text_repel(
    data = data_para, 
    aes(label = Prey.Standard_name),
    color = 'black',
    size = 4,
    box.padding = unit(0.4, "lines"),
    point.padding = unit(0.3, "lines"),
    min.segment.length = 0,
    segment.curvature = -0.1,
    segment.ncp = 3,
    segment.angle = 20
  )+
  scale_color_manual(values = c('grey65','darkcyan'), 
                     labels = c('SH3-independent','SH3-dependent'))+
  scale_shape_discrete(labels = c('equivalent', 'divergent'), solid = T)+
  theme(legend.position = 'bottom',  
        strip.text.y = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        axis.title = element_text(size =16), 
        axis.text = element_text(size = 14), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(title = "",
                              override.aes = aes(label = ""), nrow = 2), 
         shape = guide_legend(title = "",
                              override.aes = aes(label = ''), nrow = 2))

#ggsave('~/ancSH3_paper/Reviews/FigurePanels/Fig2B.png')

#Figure 2D & E:  show differences between DupSH3, optSH3s and swapSH3s
swapm3 <- subset(data,
               subset = sh3_sequence %in% c('optMyo3', 'swapMyo5'))
swapm5 <- subset(data,
                 subset = sh3_sequence %in% c('optMyo5', 'swapMyo3'))
dupm3 <- subset(data,
                 subset = sh3_sequence %in% c('optMyo3', 'AncA'))
dupm5 <- subset(data,
                subset = sh3_sequence %in% c('optMyo5', 'AncA'))

swapm3 <- 
  pivot_wider(unique(swapm3[swapm3$Bait.Standard_name == 'Myo3' , c(1,2,3,4, 16:17)]), 
              names_from = sh3_sequence, 
              values_from = c(med.PPI_score))

swapm5 <- 
  pivot_wider(unique(swapm5[swapm5$Bait.Standard_name == 'Myo5' ,  c(1,2,3,4, 16:17)]), 
              names_from =  sh3_sequence, 
              values_from = c(med.PPI_score))


dupm3 <- 
  pivot_wider(unique(dupm3[dupm3$Bait.Standard_name == 'Myo3' , c(1,2,3,4, 16:17)]), 
              names_from = sh3_sequence, 
              values_from = c(med.PPI_score))

dupm5 <- 
  pivot_wider(unique(dupm5[dupm5$Bait.Standard_name == 'Myo5' , c(1,2,3,4, 16:17)]), 
              names_from = sh3_sequence, 
              values_from =  c(med.PPI_score))

sm3_dep <- c('Bzz1', 'Sla1')
sm5_dep <- c('Bzz1', 'Ysc84')
ancAm3_dep <- c('Bzz1', 'Lsb3', 'Ysc84', 'Sla1', 'Syp1', 'Pkh2', 'Ste20', 'Yor1')
ancAm5_dep <- c('Mid2','Myo5', 'Srv2','Ysc84', 'Syp1', 'Rvs167', 'Osh2')

swapm3$sh3_spe <- swapm3$Prey.Standard_name %in% sm3_dep
swapm5$sh3_spe <- swapm5$Prey.Standard_name %in% sm5_dep
dupm3$sh3_spe <- dupm3$Prey.Standard_name %in% ancAm3_dep
dupm5$sh3_spe <- dupm5$Prey.Standard_name %in% ancAm5_dep


ld <- 
  ggplot(swapm3[swapm3$Prey.Standard_name %in% signif_para, ], aes(optMyo3, swapMyo5)) +
  geom_point(data = data_para[data_para$Prey.Standard_name %in% signif_para, ],
             aes(x = med.PPI_score_extantMyo3, y = med.PPI_score_extantMyo5), 
             color = 'black', size = 3 , shape = 17, alpha = 0.5)+
  geom_point(size = 3, color = '#A20056', aes(shape = sh3_spe))+
  stat_cor(data = swapm3[swapm3$Prey.Standard_name %in% signif_para, ], 
             size = 4.5, method = 'spearman', cor.coef.name = c('r'), 
           color = '#A20056', label.y = 0.97)+
  stat_cor(data = data_para[data_para$Prey.Standard_name %in% signif_para, ], 
           aes(x = med.PPI_score_extantMyo3, y = med.PPI_score_extantMyo5),
           size = 4.5, method = 'spearman', cor.coef.name = c('r'), 
           color = 'grey20', label.y = 0.87)+
  facet_grid(cols = vars(Bait.Standard_name))+
  theme_bw() +
  scale_x_continuous(labels = c(0, 0.25, 0.5, 0.75, 1))+
  scale_y_continuous(labels = c(0, 0.25, 0.5, 0.75, 1))+
  xlab('med. optMyo3 PPI score') +
  ylab('med. swapMyo5 PPI score') +
  lims(x = c(0,1), y = c(0,1))+  
    theme(plot.title = element_text(hjust = 0.5, size =16),
          legend.position = 'none', 
          legend.text = element_text(size = 14), 
          axis.title = element_text(size =14), 
          strip.text.x = element_text(size = 16), 
          axis.text = element_text(size =12),
          strip.text.y = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    guides(color = guide_legend(title = "",
                                override.aes = aes(label = "")))


rd <- 
  ggplot(swapm5[swapm5$Prey.Standard_name %in% signif_para, ]) +
  geom_point(data = data_para[data_para$Prey.Standard_name %in% signif_para, ],
             aes(x = med.PPI_score_extantMyo5, y = med.PPI_score_extantMyo3), 
             color = 'black', size = 3 , shape = 17, alpha = 0.5)+
  geom_point(size = 3, color = '#A20056', aes(optMyo5, swapMyo3, shape = sh3_spe))+
  stat_cor(data = data_para[data_para$Prey.Standard_name %in% signif_para, ], 
           aes(x = med.PPI_score_extantMyo5, y = med.PPI_score_extantMyo3),
           size = 4.5, method = 'spearman', cor.coef.name = c('r'), 
           color = 'grey20', label.y = 0.87)+
  stat_cor(data = swapm5[swapm5$Prey.Standard_name %in% signif_para, ],
           aes(x= optMyo5, y = swapMyo3),
           size = 4.5, method = 'spearman', cor.coef.name = c('r'), 
           color = '#A20056', label.y = 0.97)+
  facet_grid(cols = vars(Bait.Standard_name))+
  theme_bw() +
  scale_x_continuous(labels = c(0, 0.25, 0.5, 0.75, 1))+
  scale_y_continuous(labels = c(0, 0.25, 0.5, 0.75, 1))+
  xlab('med. optMyo5 PPI score') +
  ylab('med. swapMyo3 PPI score') +
  lims(x = c(0,1), y = c(0,1))+  
  theme(plot.title = element_text(hjust = 0.5, size =16),
        legend.position = 'none', 
        legend.text = element_text(size = 14), 
        axis.title = element_text(size =14), 
        strip.text.x = element_text(size = 16), 
        axis.text = element_text(size =12),
        strip.text.y = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(title = "",
                              override.aes = aes(label = "")))


le <- 
  ggplot(dupm3[dupm3$Prey.Standard_name %in% signif_para, ])+
  geom_point(size = 3, color = '#9bcd9bff', aes(optMyo3, AncA, shape = sh3_spe))+
  stat_cor(data = dupm3[dupm3$Prey.Standard_name %in% signif_para, ],
           aes(x= optMyo3, y = AncA),
           size = 4.5, method = 'spearman', cor.coef.name = c('r'), 
           color = '#9bcd9bff', label.y = 0.97)+
  facet_grid(cols = vars(Bait.Standard_name))+
  theme_bw() +
  scale_x_continuous(labels = c(0, 0.25, 0.5, 0.75, 1))+
  scale_y_continuous(labels = c(0, 0.25, 0.5, 0.75, 1))+
  xlab('med. optMyo3 PPI score')+
  ylab('med. AncA PPI score') +
  lims(x = c(0,1), y = c(0,1))+  
  theme(plot.title = element_text(hjust = 0.5, size =16),
        legend.position = 'none', 
        legend.text = element_text(size = 14), 
        axis.title = element_text(size =14), 
        strip.text.x = element_text(size = 16), 
        axis.text = element_text(size =12),
        strip.text.y = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(title = "",
                              override.aes = aes(label = "")))

re <- 
ggplot(dupm5[dupm5$Prey.Standard_name %in% signif_para, ]) +
  geom_point(size = 3, color = '#9bcd9bff', aes(optMyo5, AncA, shape = sh3_spe))+
  stat_cor(data = dupm5[dupm5$Prey.Standard_name %in% signif_para, ],
           aes(x= optMyo5, y = AncA),
           size = 4.5, method = 'spearman', cor.coef.name = c('r'), 
           color = '#9bcd9bff', label.y = 0.97)+
  facet_grid(cols = vars(Bait.Standard_name))+
  theme_bw() +
  scale_x_continuous(labels = c(0, 0.25, 0.5, 0.75, 1))+
  scale_y_continuous(labels = c(0, 0.25, 0.5, 0.75, 1))+
  xlab('med. optMyo5 PPI score') +
  ylab('med. AncA PPI score') +
  lims(x = c(0,1), y = c(0,1))+  
  theme(plot.title = element_text(hjust = 0.5, size =16),
        legend.position = 'none', 
        legend.text = element_text(size = 14), 
        axis.title = element_text(size =14), 
        strip.text.x = element_text(size = 16), 
        axis.text = element_text(size =12),
        strip.text.y = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text.y = element_blank(), 
        panel.grid.minor = element_blank())+
  guides(color = guide_legend(title = "",
                              override.aes = aes(label = "")))

# assembly of Figure 2 E & D

d <- 
plot_grid(ld+theme(legend.position = 'none'), rd, labels = 'D', label_fontface = 'plain', label_size = 16)

e <- 
plot_grid(le, re, rel_widths = c(1, 0.9),labels = 'E', label_fontface = 'plain', label_size = 16)

Figure2ED <- 
  plot_grid(d,e, rel_widths = c(1, 0.9))



# Figure 2A : protein abundance measured by cytometry, Table S6
GFP_data <- 
  read.csv('~/ancSH3_paper/Reviews/SupplementaryMaterial/TableS6.csv', header = T)[, -1]

GFP_data$strain <- str_c(GFP_data$sh3_sequence, GFP_data$GFP_tagged_protein, sep ='.')

# Keep a subset of the data
test <- 
  subset(GFP_data, 
         subset = strain %in% c('extantSH3.Myo3', 'SH3-depleted.Myo3', 'extantSH3.Myo5', 'SH3-depleted.Myo5', 'neg.control.neg.control'))

test$sh3 <- factor(test$strain, 
                   levels = c('neg.control.neg.control', 'extantSH3.Myo3', 'SH3-depleted.Myo3', 'extantSH3.Myo5', 'SH3-depleted.Myo5'), 
                   labels = c('neg.control', 'extantMyo3', 'SH3-depleted Myo3', 'extantMyo5', 'SH3-depleted Myo5'))
test$stuf <- 
  grepl(test$sh3, pattern = 'depleted') 

test$label <- 
  factor(test$sh3, 
         levels = c('neg.control', 'extantMyo3', 'SH3-depleted Myo3', 'extantMyo5', 'SH3-depleted Myo5'), 
         labels = c('neg.control', 'extant', 'SH3-depleted', 'extant', 'SH3-depleted'))


my_comparison = list(c('extantMyo3', 'SH3-depleted Myo3'), c('extantMyo5', 'SH3-depleted Myo5'), c('extantMyo3', 'extantMyo5'))

# Student't t-tests
pval <-
 compare_means(med_rep ~ sh3, test, method = 't.test', p.adjust.method = 'BH')[c(5,6,10), ]

pval$my.pval <- paste0('p = ', pval$p.adj)

Fig2A <- 
  ggplot(test, aes(y = log2(med_rep), x = sh3 , color = label))+
  geom_jitter( alpha=1, width = 0.2, size = 3)+
  ylab(expression(paste('med. fluorescent intensity(', log[2], ' scale)')))+
  ylim(0, 2.3)+
  scale_color_manual(values = c('grey30', 'grey65', '#57A337'), labels = c('control', 'extant', 'SH3-depleted'))+
  stat_pvalue_manual(pval, label = "my.pval", y.position= c(2, 2.2, 2), size = 4.5)+
  theme_bw()+
  scale_x_discrete(limits =c('neg.control', 'extantMyo3', 'SH3-depleted Myo3', 'extantMyo5', 'SH3-depleted Myo5'))+
  theme(legend.position = 'bottom', 
        legend.title = element_blank(), 
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 14), 
        axis.text = element_text(size =12),
        axis.text.x = element_text(hjust = 1, angle = 30),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = 'none')



# Figure Supplementary 2G : All cytometry data
median(GFP_data[ GFP_data$GFP_tagged_protein =='Myo5', "med_rep"])
#3.46

median(GFP_data[ GFP_data$GFP_tagged_protein =='Myo3', "med_rep"])
#2.5

x <-
  compare_means(med_rep~sh3_sequence, group.by = 'GFP_tagged_protein', data = GFP_data, method = 't.test', ref.group = 'extantSH3', 
                   p.adjust.method = 'BH')

y <-
compare_means(med_rep~GFP_tagged_protein, data = GFP_data[GFP_data$sh3_sequence %in% c('extantSH3', 'control'), ], method = 't.test',  
              p.adjust.method = 'BH')[-1, ]




SuppFig2GFP <- 
ggplot(GFP_data, aes(x = sh3_sequence, y = med_rep, color = GFP_tagged_protein))+
  geom_jitter(alpha=0.8, width = 0.2, size = 3)+
  scale_color_manual(values = c( '#3366CC', 'orangered','grey30'))+
  geom_hline(yintercept = c((2.5), (3.46)), linetype = 'dashed')+
   xlab('SH3 domain')+
   ylab(expression(paste('med. fluorescent intensity(', log[2], ' scale)')))+
   scale_y_continuous(trans = 'log2', breaks = c(1.5, 2,2.5,3, 3.5, 4))+
   scale_x_discrete(limits = c('neg.control', 'extantSH3','optSH3', 'AncA', 'AncA_3', 'AncB', 'AncC', 'AncD',  'swap.optSH3' , 'swap.control', 'swapSH3', 'SH3-depleted'), 
                    labels = c('neg.control', 'extantSH3','optSH3', 'AncA', 'AncA_3', 'AncB', 'AncC', 'AncD',  'swap.optSH3' , 'swap.control', 'swapSH3', 'SH3-depleted'))+
   theme_bw()+
   geom_text(data = x, aes(group2, y = rep(c(2.8, 4), each = 11), label = p.adj), 
             show.legend = F, size = 3)+
  #geom_text(data = y, aes(x = group2, y = rep(c(2.8, 4), each = 2), label = p.adj),
   #         show.legend = F, color = 'grey30', size = 3)+
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(hjust = 1, angle = 30),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(nrow = 1))

saveRDS(SuppFig2GFP, file = '~/ancSH3_paper/Reviews/FigurePanels/SuppFig2G.rds')

# Assembly of Figure 2

Fig2C <- c('~/ancSH3_paper/Reviews/FigurePanels/Fig2C.png')
Fig2C <-ggdraw() + #create blank canvas
  draw_image(Fig2C,scale=1)

top <- 
plot_grid(Fig2A+theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size =12)), 
          Fig2B+theme(legend.text = element_text(size = 14),
                      axis.title = element_text(size = 14),
                      axis.text = element_text(size =12)),
           rel_widths = c(0.75, 1.2),
          align = 'h', axis = 'tb',
          ncol = 2, labels = c('A', 'B'), label_fontface = 'plain', label_size = 16)

top2<-
  plot_grid(top, Fig2C, rel_widths = c(2, 1.2),
            ncol = 2, labels = c('', 'C'), label_fontface = 'plain', label_size = 16, align = 'h', 
            axis = 't')


plot_grid(top2, Figure2ED,
          ncol = 1, nrow = 2, 
          label_fontface = 'plain', label_size = 16, 
          rel_heights = c(1, 0.6))
# Save Figure 2
ggsave('~/ancSH3_paper/Reviews/Figure2.svg', 
       height = 8, width = 12)



