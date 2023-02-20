######################################
## Figure 1 & Supp Figure 3         ##
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
library(gtools)
library(xlsx)

firstup <- function(x) {
  substring(x, 2) <- tolower(substring(x, 2))
  x
}

# Import sheet1 (full-lenght PCA) SupplementaryData1
data <-
  read.csv('~/ancSH3_paper/SupplementaryMaterial/TableS1.csv')[,-1]

# Import abundance controls information
abondance_c <-
  c(
    "Hek2",
    "Cdc10",
    "Abp1"   ,
    "Pab1"   ,
    "Bmh1"   ,
    "Ade57" ,
    "Tif4631",
    "Tdh3",
    "Pfk1",
    "Sbp1" ,
    "Ssz1",
    "Hsp104" ,
    "Ssb2",
    "Hsp82"
  )


# Modification of data by columns to fit figures annotation
data$Bait.Standard_name <- as.factor(data$Bait.Standard_name)

# confirmed proline motif on prey by PCA
motif_conf <-
  c('Pkh2',
    'Sla1',
    'Myo5',
    'Osh2',
    'Ste20',
    'Lsb3',
    'Bzz1',
    'Syp1',
    'Myo3')

# Figure 1B: SH3-depletion effect on PPI score
stuf <- subset(data,
               subset = sh3_sequence %in% c('extantMyo3', 'extantMyo5', 'SH3-depleted'))

stuf$sh3_sequence <- 
  factor(stuf$sh3_sequence, 
         levels = c('extantMyo3', 'extantMyo5', 'SH3-depleted'), 
         labels = c('WT', 'WT', 'SH3-depleted'))


wilcox.test(PPI_score ~ sh3_sequence, stuf)

data_stuf <-
  pivot_wider(unique(stuf[, c(1:4, 16,17,19)]),
              names_from = sh3_sequence,
              values_from = c(med.PPI_score, SH3_dep, Low_Rep))

# Remove data containing NAs values for SH3-depleted PPI score
data_stuf <- 
  data_stuf[!is.na(data_stuf$`med.PPI_score_SH3-depleted`), ]

# column with label used in plots
data_stuf$label <- ''
 
data_stuf[data_stuf$SH3_dep_WT | data_stuf$`Low_Rep_SH3-depleted`, "label"] <-
 data_stuf[data_stuf$SH3_dep_WT | data_stuf$`Low_Rep_SH3-depleted`, "Prey.Standard_name"]


pvalm3 <- 
  wilcox.test(PPI_score~sh3_sequence,
              data = stuf[stuf$Bait.Standard_name == 'Myo3', ])

pvalm5 <- 
  wilcox.test(PPI_score~sh3_sequence,
              data = stuf[stuf$Bait.Standard_name == 'Myo5', ])



# Figure 1C
p1 <-
  ggplot(data_stuf[!data_stuf$`Low_Rep_SH3-depleted`, ], aes(med.PPI_score_WT, `med.PPI_score_SH3-depleted`)) +
  facet_grid(cols = vars(Bait.Standard_name), scales = 'free') +
  geom_point(aes(color = SH3_dep_WT, shape = `Low_Rep_SH3-depleted`),size = 3) +
  geom_point(data = data_stuf[data_stuf$`Low_Rep_SH3-depleted`, ],
             aes(med.PPI_score_WT, `med.PPI_score_SH3-depleted`, shape = `Low_Rep_SH3-depleted`), size = 3, color = 'grey20')+   
  theme_bw() +
  stat_cor(aes(med.PPI_score_WT, `med.PPI_score_SH3-depleted`, color = SH3_dep_WT), size = 5, method = 'spearman', 
           cor.coef.name = c('r')) +
  ylab('med. SH3-depleted PPI score') +
  xlab('med. WT PPI score') +
  xlim(0,1)+
  ylim(0,1)+
  geom_text_repel(
    data = data_stuf,
    aes(label = label),
    color = 'black',
    size = 4,
    box.padding = unit(0.4, "lines"),
    point.padding = unit(0.4, "lines"),
    segment.curvature = 0.1,
    segment.ncp = 4,
    segment.angle = 20,
    min.segment.length = 0, 
    nudge_x = 0.07 
  ) +
  scale_color_manual(values = c('#CCCCCC', 'darkcyan'), labels = c('SH3-independent', 'SH3-dependent')) +
  theme(legend.position = 'bottom',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=16), 
        axis.text = element_text(size =14),
        legend.text = element_text(size = 16), 
        axis.title = element_text(size =16)) +
  guides(color = guide_legend(title = "",
                            override.aes = aes(label = "")), 
       shape = 'none')

#Figure 1D : proline motif prediction

pred_motif <- 
  read.csv('~/ancSH3_paper/SupplementaryMaterial/TableS5.csv')[, -1]

data_stuf %<>% 
  mutate(deltaPPI = med.PPI_score_WT - `med.PPI_score_SH3-depleted`)

pred_motif$Prey.Standard_name <- 
  firstup(pred_motif$Prey.Standard_name)

pred_motif <- 
merge(data_stuf, 
      pred_motif, 
      by.x = c('Bait.Standard_name','Prey.Standard_name'), 
      by.y = c('Paralog', 'Prey.Standard_name'))

pred_motif <- 
  as_tibble(pred_motif)

colnames(pred_motif)[16] <- 'p-value'

p2 <- 
ggplot(pred_motif,
       aes(Max_MSS, deltaPPI, color = `p-value`))+
  facet_grid(cols = vars(Bait.Standard_name))+
  theme_bw()+
  geom_point(size = 3)+
  #scale_fill_gradient()
  #scale_color_viridis(option = 'mako', discrete = F)+
  stat_cor(size = 5, method = 'spearman', cor.coef.name = c('r'))+
  scale_color_gradientn(colours = c('#0B0405FF', '#2E1E3CFF', '#413D7BFF', '#37659EFF', '#348FA7FF', '#40B7ADFF', '#8AD9B1FF', '#DEF5E5FF'),
                        trans = 'log', breaks=c(0,0.001, 0.01,0.1, 1))+
  ylab(bquote(atop(''*Delta~'('~PPI[WT]~ - ~PPI[depleted]~')', 'med. score')))+
  xlab('Max MSS score')+
  ylim(0,1)+
  xlim(0,1)+
  #scale_x_continuous(trans = 'log10', breaks = c(0.3,0.5,0.6, 0.7, 0.8, 0.9, 1))+
  geom_smooth(method = 'glm', color = 'orangered', fill = 'grey65')+
  geom_text_repel(data=pred_motif[pred_motif$Prey.Standard_name %in% motif_conf, ],
    aes(label = Prey.Standard_name),
    color = 'black',
    size = 4,
    box.padding = unit(0.4, "lines"),
    point.padding = unit(0.4, "lines"),
    segment.curvature = 0.1,
    segment.ncp = 4,
    segment.angle = 20,
    min.segment.length = 0, 
    nudge_x = 0.07
  )+
  theme(legend.position = 'bottom',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
    strip.text.x = element_blank(), 
    axis.text = element_text(size =12),
    legend.text = element_text(size = 12), 
    axis.title = element_text(size =14), 
    legend.title = element_text(size =16), 
    legend.key.height  = unit(0.4, 'cm'),
     legend.key.width = unit(1.2, 'cm'))#+
  #guides(color = guide_legend(title.hjust = -0.3))

viridis(n=n, option = LETTERS[7])


# Figure 1E: proline motif confirmation PCA

data_motif <- 
  read.csv('~/ancSH3_paper/SupplementaryMaterial/TableS2.csv')[, -1]

data_motif <- 
subset(data_motif, 
       subset = sh3_sequence != 'UNI1')

data_motif_w <- 
pivot_wider(unique(data_motif[, c(1,2,3,4,5,20)]), 
                   values_from = med.PPI_score, 
            names_from = motif_deletion) 

colnames(data_motif_w)[c(5,6)] <- c('extantPrey', 'motifdepletedPrey')

data_motif_w %<>%
  mutate(dPPI.score = extantPrey - motifdepletedPrey)

## Save file used for Figure 1B 
# write.csv2(data_motif_w, file = '~/ancSH3_paper/SupplementaryMaterial/Data/Fig1A_Cytoscape_input.csv')
  
data_fig <- 
  subset(data_motif, sh3_sequence %in% c('extantMyo3', 'extantMyo5'))

p3 <- 
ggplot(data_fig, aes(x=Prey.Standard_name, y = PPI_score, 
                     fill = motif_deletion))+
  facet_grid(col = vars(Bait.Standard_name))+
  geom_boxplot(outlier.colour = 'grey20', outlier.size = 1, position=position_dodge(width=1))+
  stat_compare_means(aes(group = motif_deletion), 
                     label = "p.signif", 
                     bracket.size = 'none', 
                     hide.ns = T,
                     method.args = list(alternative = "less"), 
                     size = 4, label.y = 0.97)+
  theme_bw()+
  xlab('Preys')+
  ylab('PPI score')+
  scale_fill_manual(values = c('grey65', '#E5CA28'), labels = c('WT prey', bquote('motif'*Delta~'prey')))+
  theme(legend.position = 'bottom', 
        strip.text.x = element_blank(), 
        axis.text.x = element_text(size =14, hjust = 1, angle=30),
        axis.text.y = element_text(size =14),
        legend.text = element_text(size = 16), 
        axis.title = element_text(size =16), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
    guides(fill = guide_legend(title = "",
                              override.aes = aes(label = "")), 
         color = 'none')


## Panel assembly for Figure 1
library(cowplot)

# Network created with cytoscape and modified by inkscape
network <- c('~/ancSH3_paper/SupplementaryMaterial/FigurePanels/Fig1A.png')
a <- 
  ggdraw()+
  draw_image(network)


left <- plot_grid(p1+theme(legend.text = element_text(size = 14),
                           legend.title = element_text(size = 16),
                           axis.title = element_text(size = 14),
                           axis.text = element_text(size =12)), 
                  p2+theme(legend.text = element_text(size = 12),
                           legend.title = element_text(size = 16),
                           axis.title = element_text(size = 14),
                           axis.text = element_text(size =12),
                           axis.title.y = element_text(hjust = 0)), 
                  p3+theme(legend.text = element_text(size = 14),
                           legend.title = element_text(size = 16),
                           axis.title = element_text(size = 14),
                           axis.text = element_text(size =12)), 
                  nrow = 3, labels = c('C', 'D', 'E'), rel_heights = c(-1,-1,-1),
          align = 'v', axis = 'lr', label_size = 16, label_fontface = 'plain')


b <- ggdraw()+
  draw_image('~/ancSH3_paper/SupplementaryMaterial/FigurePanels/Fig1B.png')

right <- 
  plot_grid(a, b, 
                   labels = c('A', 'B'), 
                   rel_heights = c(4, 1), 
                   label_size = 16, label_fontface = 'plain', 
            ncol =1)

plot_grid(right, left, 
          ncol = 2, labels = c('', ''), 
          label_size = 16, label_fontface = 'plain', 
          align = 'h',axis = 't', rel_widths = c(1, 1.1))

ggsave('~/ancSH3_paper/Figure1.png', 
       height = 13, width = 15)
# Panel labels (A, B) are moved manually on Inkscape 


# Supplementary Figure 3C : proline motif impact on the 
# interaction with other SH3 sequence

data_motif$sh3_sequence <- 
  factor(data_motif$sh3_sequence, 
         levels = c('extantMyo3', 'extantMyo5', 'optMyo3', 'optMyo5', 'AncC', 'SH3-depleted'), 
         labels = c('extantSH3', 'extantSH3', 'optSH3', 'optSH3', 'AncC', 'SH3-depleted') )

sf3c<-
  ggplot(data_motif[data_motif$sh3_sequence %in% c('AncC', 'SH3-depleted'),],
         aes(x=Prey.Standard_name, y = PPI_score, fill = motif_deletion))+
  facet_grid(rows= vars(Bait.Standard_name), 
             cols = vars(sh3_sequence))+
  geom_boxplot(outlier.colour = 'grey20', outlier.size = 1, position=position_dodge(width=1))+
    stat_compare_means(aes(group = motif_deletion), 
                     label = "p.signif", 
                     bracket.size = 'none', 
                     hide.ns = T,
                     size = 4, label.y = 0.96)+
  theme_bw()+
  xlab('Preys')+
  ylab('PPI score')+
  scale_fill_manual(values = c('grey65', '#E5CA28'), labels = c('WT prey', bquote('motif'*Delta~'prey')))+
  #scale_color_manual(values = c('grey65', '#9BCD9B'))+
  theme(legend.position = 'bottom', 
        #strip.text.x = element_blank(), 
        axis.text.x = element_text(size =12, hjust = 1, angle=30),
        axis.text.y = element_text(size =12),
        legend.text = element_text(size = 14), 
        axis.title = element_text(size =14),
        strip.text = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  guides(fill = guide_legend(title = "",
                             override.aes = aes(label = "")), 
         color = 'none')

# Supplementary Figure 3B : proline motif position on prey
library('drawProteins')
library(ggplot2)

# get the protein sequences
#bzz1, lsb3, Myo5, Osh2, pkh2, ste20, sla1

prot <- 
  get_features('P38822 P43603 Q04439 Q12451 Q12236 Q03497 P32790')

# create the sequence annotation for pxxp motifs
# predicted max motif
prey <- 
  toupper(c('bzz1', 'lsb3', 'Myo5', 'Osh2', 'pkh2', 'ste20', 'sla1'))

#566-575, 196-205, 1072-1081, 771-780, 633-642, 474-483, 620-629

begin <- c(566, 196, 1072, 771, 633, 474, 620)
end <- c(575, 205, 1081, 780, 642, 483, 629)

max_motif <- 
  data.frame(prey, begin, end)

#loop to create the annotations in protein object
for (i in 1:length(prot)) {
  beg <- gregexpr(pattern = '.P..P.', text = prot[[i]]$sequence, fixed = F)
  end <- beg[[1]]+4
  
  for (j in 1:length(beg[[1]])) {
    
    mot <-   
      list(list('type' = 'motif', 
                'category' = 'xx', 
                'begin'= '', 
                'end' = '', 
                'molecule' = '', 
                'evidence' = list()))
    
    mot[[1]]$begin <- beg[[1]][j]
    mot[[1]]$end <- end[j]
    
    prot[[i]]$features <-
      append(prot[[i]]$features, mot)
    
  } 
  
  
  prey_n <- unlist(strsplit(prot[[i]]$entryName, split = '_', fixed = T))[1]
  
  mot_max <-   
    list(list('type' = 'motif_max', 
              'category' = 'xx', 
              'begin'= '', 
              'end' = '', 
              'molecule' = 'max', 
              'evidence' = list()))
  
  r <- grep(prey_n, max_motif$prey)
  
  mot_max[[1]]$begin <-  max_motif[r, "begin"]
  mot_max[[1]]$end <- max_motif[r, "end"]
  
  prot[[i]]$features <-
    append(prot[[i]]$features, mot_max)
  
}

data_prot <- 
  feature_to_dataframe(prot)

data_prot[data_prot$type == 'motif', "description" ] <- 'proline motif'
data_prot[data_prot$type == 'motif_max', "description" ] <- 'max motif prediction'

data_prot[data_prot$type %in% c('motif', 'motif_max'), "type" ] <- 'MOTIF'

data_prot <- 
  data_prot[!data_prot$description == 'FFAT',  ] 


p <- 
  draw_canvas(data_prot)+
  theme_bw()

p <- 
  draw_chains(p, data_prot, fill = 'grey20', 
              labels = c('Sla1', 'Bzz1', 'Lsb3', 'Ste20', 'Myo5', 'Pkh2', 'Osh2'))

data_prot[data_prot$type == 'MOTIF', "end"] <- 
data_prot[data_prot$type == 'MOTIF', "end"]+3

data_prot[data_prot$type == 'MOTIF', "begin"] <- 
  data_prot[data_prot$type == 'MOTIF', "begin"]-3

sf3b <- 
  draw_motif(p, data_prot)+
  scale_fill_manual(values = c('#e5ca28ff', 'grey65'), labels = c('Max MSS motif', 'proline motif (PXXP)'))+
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank()) +
  theme(axis.ticks = element_blank(), 
        axis.text.y = element_blank()) +
  theme(panel.border = element_blank())+
  theme(legend.position = 'bottom', 
        strip.text.x = element_blank(), 
        legend.title = element_blank(), 
        axis.text.x = element_text(size =14),
        #axis.text.y = element_text(size =14),
        legend.text = element_text(size = 16), 
        axis.title = element_text(size =16))

# Supplementary Figure 3A : correlation PPI score and prediction score
sf3a <- 
ggplot(pred_motif,
       aes(Max_MSS, med.PPI_score_WT, color = `p-value`))+
  facet_grid(cols = vars(Bait.Standard_name))+
  theme_bw()+
  geom_point(size = 3)+
  stat_cor(size = 5, method = 'spearman', cor.coef.name = c('r'))+
  scale_color_gradientn(colours = c('#0B0405FF', '#2E1E3CFF', '#413D7BFF', '#37659EFF', '#348FA7FF', '#40B7ADFF', '#8AD9B1FF', '#DEF5E5FF'),
                        trans = 'log', breaks=c(0,0.001, 0.01,0.1, 1))+
  ylab('med. PPI score')+
  xlab('Max MSS score')+
  ylim(0,1)+
  xlim(0,1)+
  geom_smooth(method = 'glm', color = 'orangered', fill = 'grey65')+
  geom_text_repel(data=pred_motif[pred_motif$Prey.Standard_name %in% motif_conf, ],
                  aes(label = Prey.Standard_name),
                  color = 'black',
                  size = 4,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"),
                  segment.curvature = 0.1,
                  segment.ncp = 4,
                  segment.angle = 20,
                  min.segment.length = 0, 
                  nudge_x = 0.07
  )+
  
  theme(legend.position = 'right',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x = element_blank(), 
        axis.text = element_text(size =12),
        legend.text = element_text(size = 14), 
        axis.title = element_text(size =14), 
        legend.title = element_text(size =16), 
        legend.key.width = unit(1, 'cm'))

# Assembly of Figure S3


side <-
  plot_grid(sf3a+theme(legend.key.height  = unit(0.4, 'cm'),
                       legend.key.width = unit(1.2, 'cm'),
                       legend.position = 'bottom'),
            sf3b+theme(legend.position = 'bottom', 
                       legend.text = element_text(size = 14)), 
                  ncol = 1, nrow = 2, labels = c('A', 'B'), 
                  label_size = 16, label_fontface = 'plain', 
                  #align = 'v', axis = 'l',
                  rel_heights = c(1, 1))

plot_grid(side, 
          sf3c+theme(legend.position = 'bottom'),
          ncol = 2, label_size = 16, label_fontface = 'plain', 
          labels = c('', 'C'), rel_widths = c(1, 1))


ggsave('~/ancSH3_paper/SupplementaryMaterial/FigureS3.png', 
       width = 13, height = 7.5)

