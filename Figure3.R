######################################
## Figure 3                         ##
## author : Pascale Lemieux         ##
## Date : 2023-01-12                ##
######################################

library(ggplot2)
library(dplyr)
library(gtools)
library(viridis)
library(cowplot)

firstup <- function(x) {
  substring(x, 2) <- tolower(substring(x, 2))
  x
}

# Read abundance control information
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

# Read Table S1
PCA_complete <- 
  read.csv('~/ancSH3_paper/SupplementaryMaterial/TableS1.csv')[, -1]

# Verify that all preys have standard names

PCA_complete[PCA_complete$Prey.Standard_name == '', 'Prey.Standard_name'] <- 
  firstup(PCA_complete[PCA_complete$Prey.Standard_name == '', 'Prey.Systematic_name'])


# Figure 3C : all SH3s PPI score heatmap

# arrange data for heatmap
# keep significant comparison vs optSH3s 
signif <- 
  unique(na.omit(PCA_complete[PCA_complete$optSH3_dif, c(2,3,4)]))

signif$signif <- T

# reorder the sh3_sequence factor and keep relevant sh3_sequence data
# optSH3s + ancestral SH3s
PCA_complete$sh3_sequence <- 
  as.character(PCA_complete$sh3_sequence)
data <- 
  subset(PCA_complete, 
         subset = sh3_sequence %in% c(
           'SH3-depleted',
           'AncC',
           'AncB',
           'AncA',
           'DupSH3',
           'optMyo3',
           'optMyo5',
           'extantMyo3', 
           'extantMyo5'

         ))

data$sh3_sequence <- 
  ordered(as.factor(data$sh3_sequence), levels=c('SH3-depleted',
                                                 'AncC',
                                                 'AncB',
                                                 'AncA',
                                                 'DupSH3','extantMyo3', 
                                                 'extantMyo5',
                                                 'optMyo3', 
                                                 'optMyo5'
                                                 
                                                 
  ))


# Remove the opt-swapSH3 f
data <- 
  data[!c(data$sh3_sequence == 'optMyo5' & data$Bait.Standard_name =='Myo3'), ]
data <- 
  data[!c(data$sh3_sequence == 'optMyo3' & data$Bait.Standard_name =='Myo5'), ]

# Set y axis label order

name <- 
unique(PCA_complete[order(PCA_complete[ , "med.PPI_score"], decreasing=T),
                    c(1:4, 16:19)])
name <- 
  name[name$sh3_sequence =='SH3-depleted',]

sh3dep <- 
  unique(name[name$SH3_dep, 2])
sh3indep <- 
  unique(name[!name$SH3_dep, 2])
sh3indep <- 
  sh3indep[!(sh3indep %in% sh3dep)]
sh3indep <- 
  sh3indep[!(sh3indep %in% abondance_c)]

abond<- 
  unique(name[name$Prey.Standard_name %in% abondance_c, 2])


ynam <- unique(c(abond, sh3indep, sh3dep))


hm <- 
ggplot(data) +
  facet_grid(
    cols = vars(Bait.Standard_name),
    scales = 'free',
    drop = T,
    as.table = F,
    space = c('free_x')
  ) +
  geom_tile(
    data = data,
    aes(
    x = sh3_sequence,
    y = Prey.Standard_name ,
    fill = (med.PPI_score)
  ))+
  geom_point(data = data[data$optSH3_dif, ], 
    aes(x = sh3_sequence, y = Prey.Standard_name), shape = 4, color = 'orangered')+
  scale_fill_viridis(option = 'mako', na.value = 'grey')+
  scale_y_discrete(limits = na.omit(ynam)) +
  #scale_x_discrete(limits=c('SH3-depleted', 'AncC', 'AncB', 'AncA', 'DupSH3', 'extantMyo3', 'extantMyo5', 'optMyo3', 'optMyo5'))+
  annotate(
    "segment",
    x = 0.5,
    xend = 7.5,
    y = 23.5,
    yend = 23.5,
    colour = "#343434",
    size = 0.8,
    linetype = 2
  ) +
  annotate("segment",
           x = 0.5,
           xend = 7.5,
           y = 14.5,
           yend = 14.5,
           colour = "#343434",
           size = 0.8,
           linetype = 2)+
  annotate(
    "rect",
    xmin = 6.5, 
    xmax = 7.5, 
    ymin = 0.5, 
    ymax = 48.5,
    colour = "orangered",
    size = 0.65,
    linetype = 2, 
    alpha = 0
  )+
  theme_bw() +
  theme(
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size =14),
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(size =16),
    legend.key.width = unit(1, "cm")
  ) +
  labs(fill = 'med. PPI score', x = 'SH3 domains', y = 'Preys')


# Assemble Figure 3

a <- 
ggdraw()+
  draw_image('~/ancSH3_paper/SupplementaryMaterial/FigurePanels/Fig3A.png', scale = 1)

b <- ggdraw()+
  draw_image('~/ancSH3_paper/SupplementaryMaterial/FigurePanels/Fig3B.png', scale = 1)

left <- 
plot_grid(a, hm, labels = c('A', 'C'), label_fontface = 'plain', label_size = 16,
          rel_heights = c(1,3.5), ncol = 1)
p3 <- 
plot_grid(left, 
          b, ncol = 2, labels = c('', 'B'), rel_widths = c(3, 1.2),
          label_fontface = 'plain', label_size = 16)


ggsave(p3, file = '~/ancSH3_paper/Figure3.svg', 
      height = 13, width = 13)

# Manual annotation are added for the missing data (grey tiles) and to highlight the
# type of prey one the right y axis
