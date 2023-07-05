######################################
## Figure 3                         ##
## author : Pascale Lemieux         ##
## Date : 2023-06-02                ##
######################################

library(ggplot2)
library(dplyr)
library(gtools)
library(viridis)
library(cowplot)
library(readxl)

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
  read_xlsx("~/ancSH3_paper/Reviews/SupplementaryMaterial/TableS1.xlsx")

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
           'AncD',
           'AncC',
           'AncB',
           'AncA',
           
           'optMyo3',
           'optMyo5',
           'extantMyo3', 
           'extantMyo5'

         ))

data$sh3_sequence <- 
  ordered(as.factor(data$sh3_sequence), levels=c('SH3-depleted',
                                                 'AncD',
                                                 'AncC',
                                                 'AncB',
                                                 'AncA',
                                                 'extantMyo3', 
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
                    c(2,4, 16:18)])
name <- 
  unique(name[name$sh3_sequence =='SH3-depleted', c(1,4)])

sh3dep <- 
  unique(name[name$SH3_dep, 1])
sh3indep <- 
  unique(name[!name$SH3_dep, 1])
sh3indep <- 
  sh3indep[!(sh3indep %in% sh3dep)]
sh3indep <- 
  unlist(sh3indep)[!(unlist(sh3indep) %in% abondance_c)]

abond<- 
  unique(unlist(name[,1])[name$Prey.Standard_name %in% abondance_c])


ynam <- unique(unlist(c(abond, sh3indep, sh3dep)))
ynam <- c(ynam[1:22], 'Cmd1', ynam[27:46], 'Ste20', ynam[47:48], 'Rvs167', 'Osh2')

saveRDS(ynam, file = '~/ancSH3_paper/Reviews/SupplementaryMaterial/Data/ynam.RDS')

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
  draw_image('~/ancSH3_paper/Reviews/Fig3A.png', scale = 1)

b <- ggdraw()+
  draw_image('~/ancSH3_paper/Reviews/Fig3B.png', scale = 1)

left <- 
plot_grid(a, hm, labels = c('A', 'C'), label_fontface = 'plain', label_size = 16,
          rel_heights = c(1,3.5), ncol = 1)
p3 <- 
plot_grid(left, 
          b, ncol = 2, labels = c('', 'B'), rel_widths = c(3, 1.2),
          label_fontface = 'plain', label_size = 16)


ggsave(p3, file = '~/ancSH3_paper/Reviews/Figure3.svg', 
      height = 13, width = 13)

# Manual annotation are added for the missing data (grey tiles) and to highlight the
# type of prey one the right y axis
