######################################
## Figure 5 & SuppFigure5           ##
## author : Pascale Lemieux         ##
## Date : 2023-01-16                ##
######################################

library(gtools)
library(tidyverse)
library(viridis)
library(ggpubr)
library(magrittr)
library(scales)
library(ggrepel)
library(cowplot)

## Analysis of free SH3 PCA :

# find best concentration for each prey and keep the PPI scores
# at those concentrations
sdata <- 
  readRDS('~/ancSH3_paper/SupplementaryMaterial/Data/SupplementaryData1_freeSH3.rds')

colnames(sdata)[6] <- 'estradiol'

# Keep data of wt preys, discard delta motifs preys
udata <- 
  unique(sdata[!sdata$motif_deletion, c(1:6, 20)])

# Create the dataframe with estradiol concentration as columns

difdata <- 
  pivot_wider(udata, names_from = `estradiol`, values_from = PPI_score, names_sort = T)

# Keep data of delta motif preys
ddata <- 
  unique(sdata[sdata$motif_deletion, c(1:6, 20)])

ddata$DHFR12_tag <- factor(ddata$DHFR12_tag , 
                             levels = c('C', 'N'), 
                             labels = c('C-tag', 'N-tag'))

# Find the estradiol concentration that correspond to the median
# max slope of PPI score ~ [estradiol]
j = 'difdata'
z <- get(j)
y <- matrix()

# Compute the slope between [estradiol] and [estradiol + 2levels]
for (i in 6:13) {
  x <- z[, (i + 2)] - z[, i]
  y <- cbind(y, x)
  
}

# Set the colnames corresponding to [estradiol + 1level]
colnames(y)[-1] <- c(10, 20, 30, 40, 60, 80)

# Prepare the dataframe to compute the max median [estradiol] 
# for which the slope is max across the preys
z <-
  cbind(z[, 1:5], y[,-1])
z$max_dif <-
  colnames(z[,-c(1:5)])[apply(z[-c(1:5)], 1, which.max)]
z$max_dif <-
  as.numeric(z$max_dif)
library(dplyr)

z %>%
  dplyr::group_by(Prey.Standard_name) %>%
  mutate(est_prey = median(max_dif)) -> z


# Make sure that the median [estradiol] corresponds to
# the existing conditions tested
z[!(z$est_prey %in% c(0, 10, 20, 30 , 40 , 60, 80, 100)), "est_prey"] <-
  z[!(z$est_prey %in% c(0, 10, 20, 30 , 40 , 60, 80, 100)), "est_prey"] +5

colnames(z)[13] <- 'estradiol'

# Merge with the data to keep only the optimal 
# [estradiol] PPI score for each prey
i <-
  merge(
    udata,
    z[,-c(6:12)],
    by = c(
      'Prey.Systematic_name',
      'Prey.Standard_name',
      'sh3_sequence',
      'DHFR12_tag',
      'estradiol', 
      'motif_deletion'
    ),
    all.y = T
  )

assign(x = j, value = i)

difdata$DHFR12_tag <- factor(difdata$DHFR12_tag, 
                           levels = c('C', 'N'), 
                           labels = c('C-tag', 'N-tag'))

# write the csv file for this dataset
write.csv2(difdata, file = '~/ancSH3_paper/SupplementaryMaterial/Data/free_sh3_data/ppi_score_est_adj.csv', 
           fileEncoding = 'UTF-8')


# Create the wide form of the data frame to contain wt.prey 
# and dmotif.prey scores in different columns
dif_pxxp <- 
merge(difdata, 
      ddata, 
      by = c('Prey.Systematic_name', 
             'Prey.Standard_name',
             'sh3_sequence', 
             'estradiol',
             'DHFR12_tag'), 
      suffixes = c('.wt', '.delta'), 
      all = F)

# Compute the difference between wt.prey and dmotif.prey scores
dif_pxxp %<>%
  mutate(dif_motif = PPI_score.wt -PPI_score.delta)

write.csv2(dif_pxxp, file = '~/Maitrise/PCA_DHFR_array/PCA_SH3_libre/ppi_score_est_adj.csv', 
           row.names = F, fileEncoding = 'UTF-8')

# Merge all data points with the optimal [estradiol]

sdata$DHFR12_tag <- factor(sdata$DHFR12_tag, 
                             levels = c('C', 'N'), 
                             labels = c('C-tag', 'N-tag'))

cdata <- 
  merge(sdata, 
        difdata, 
        by = c('Prey.Systematic_name', 
               'Prey.Standard_name', 
               'sh3_sequence', 
               'DHFR12_tag'), 
        suffixes = c('.all', '.cor'))


# Verify the noise in the experiment obtain with the 
# emptyLP tagged in C and in N
difdata[difdata$sh3_sequence == 'emptyLP', ] %>%
  group_by(DHFR12_tag) %>%
  mutate(noise.tag = median(PPI_score, na.rm = T))
#noise C tag = 0.0388 N tag 0.0687

# Prepare cdata to format the figure

cdata$sh3_sequence <- factor(cdata$sh3_sequence, 
                    levels = c('emptyLP', 'AncC', 'AncB', 
                               'AncA', 'Dup3', 'Dup2', 'DupSH3', 'optMyo3', 'optMyo5', 
                               'extantMyo3', 'extantMyo5'))
cdata$estradiol.all <- 
  as.numeric(paste(cdata$estradiol.all))

# Compare the noise in the experiment with the 
# emptyLP tagged in C and in N
compare_means(PCA_score ~ DHFR12_tag, cdata[cdata$sh3_sequence == 'emptyLP', c(3,4,19)])
#.y.       group1 group2        p   p.adj p.format p.signif method  
#<chr>     <chr>  <chr>     <dbl>   <dbl> <chr>    <chr>    <chr>   
#  1 PCA_score C-tag  N-tag  1.06e-11 1.1e-11 1.1e-11  ****     Wilcoxon

cdata[cdata$estradiol.all == 0, "estradiol.all"] <- 0.99
#na.replace(cdata$estradiol.all, replace = 0.99)


Fig5SuppB <- 
ggplot(cdata[!cdata$motif_deletion.all, ])+
  facet_grid(vars(Prey.Standard_name), vars(DHFR12_tag))+
  geom_vline(aes(xintercept = as.numeric(paste(estradiol.cor))), alpha = 0.4, size = 2)+
  geom_line(aes(x = as.numeric(estradiol.all),
                y = PPI_score.all, 
                color = sh3_sequence,
                group = sh3_sequence), size =0.8)+
  scale_color_viridis(option = 'inferno', discrete = T)+
  theme_bw()+
  scale_x_continuous(trans = 'log2', breaks = c(0,10,20,30,40,60,80,100))+
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x))
                #labels = trans_format("log10", math_format(10^.x)))+
  scale_y_continuous(limits = c(0, 1), breaks = c(0.5, 1))+
  ylab('med. PPI score')+
  expand_limits(x = 0)+
  xlab(expression(paste('[Estradiol] ', mu, 'M', sep = '')))+
  theme(
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size =14),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    strip.text = element_text(size = 16), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )+
  guides(color = guide_legend(title = 'free SH3'), 
         size  = guide_legend(legend.position = ''))


# Supplementary Figure 5B : Comparison of PCA score from complete paralog vs free sh3 

# Import and adapt PCA complete paralog
dcomp_para <- readRDS('~/ancSH3_paper/SupplementaryMaterial/Data/SupplementaryData1_motif_confirmation_PCA.rds')
dcomp_para <- 
  pivot_wider(unique(dcomp_para[dcomp_para$sh3_sequence != 'UNI1' , c(2,3,4,5, 20)]), 
                          names_from = motif_deletion, 
                          values_from = PPI_score, names_prefix = c('extantPrey', 'difmotif'))
pcomp <- 
  pivot_wider(dcomp_para[dcomp_para$sh3_sequence %in% c('optMyo3', 'optMyo5') &
                           dcomp_para$Prey.Standard_name %in% difdata$Prey.Standard_name, c(1,3,4)], 
              values_from = extantPreyFALSE, 
              names_from = sh3_sequence)

pcomp$exp <- 'paralogue'

# adapt data from free sh3 PCA
comp <- 
  pivot_wider(difdata, 
              values_from = PPI_score, 
              names_from = sh3_sequence)

comp <- comp[comp$DHFR12_tag == 'C-tag', c("Prey.Standard_name", "optMyo3", "optMyo5")]
comp$exp <- 'free'

# combine both data frame
comp <- 
  rbind(comp, pcomp)

# Set factor levels labeling the PCA experiment 
comp$exp <- 
  factor(comp$exp, 
         levels = c('free', 'paralogue'), 
         labels = c(1,2))


# Figure 5A
Fig5B <- 
ggplot(comp, aes(optMyo5, optMyo3))+
  geom_point(aes(color = exp), size = 2.5)+
  lims(x = c(0,1), y = c(0,1))+
  stat_cor(data = comp[comp$exp == 1, ], label.y = 0.95, color = '#A20056', size = 4.5, cor.coef.name = c('r'), method = 'spearman')+
  stat_cor(data = comp[comp$exp == 2, ], label.y = 0.87, color = 'grey20', size = 4.5, cor.coef.name = c('r'), method = 'spearman')+
  scale_color_manual(values = c('#A20056', 'grey20'), 
                     labels = c('free SH3s', 'complete paralog'))+
  theme_bw()+
  labs(x = 'med. optMyo5 PPI score', y = 'med. optMyo3 PPI score')+
  geom_text_repel(data = comp,
    aes(label = Prey.Standard_name),
    seed= 40,
    color = 'black',
    size = 4,
    box.padding = unit(0.4, "lines"),
    point.padding = unit(0.4, "lines"),
    segment.curvature = -0.1,
    segment.ncp = 2,
    segment.angle = 10,
    min.segment.length = 0,
    nudge_x = 0.07
  )+ 
theme(
legend.position = 'bottom',
legend.text = element_text(size = 14),
legend.title = element_text(size = 16),
axis.title = element_text(size = 14),
axis.text = element_text(size =12),
axis.text.x = element_text(hjust = 0.7),
strip.text = element_text(size = 16), 
panel.grid.major = element_blank(), panel.grid.minor = element_blank()
)+
  guides(color = guide_legend(title = 'PCA experiment', 
                            title.position = 'left', title.hjust = 0.5, nrow = 2))
  


# Demonstrating that the proline motif deletion affects in the same
# way the PPI between the free SH3s and the preys vs SH3s-paralogues 
dcomp_para %>%
mutate(dPPI_score = extantPreyFALSE - difmotifTRUE) ->dcomp_para

# Merge the free SH3s data and the SH3-paralogue data
comp <- 
  merge(dif_pxxp, 
        dcomp_para, 
        by = c( 'Prey.Standard_name','sh3_sequence'))

Fig5SuppC <- 
ggplot(comp[comp$sh3_sequence %in% c('optMyo3', 'optMyo5'), ],
       aes(y = dif_motif, x = dPPI_score))+
  facet_grid(cols = vars(sh3_sequence), rows = vars(DHFR12_tag))+
  stat_cor(method = 'spearman', cor.coef.name = c('r'), size = 4)+
  ylab(expression(atop(paste(Delta,'(WT prey - motif',Delta, ' prey)'), paste('free SH3 med. PPI score'))))+
  xlab(expression(paste(Delta,'(WT prey - motif',Delta, ' prey) paralogue med. PPI score')))+
  geom_point(aes(color = Prey.Standard_name), size = 2.5)+
  xlim(-0.1, 0.6)+
  ylim(-0.1, 0.6)+
  theme_bw()+
  scale_color_viridis(option = 'mako', discrete = T)+
  geom_text_repel(
    aes(label = Prey.Standard_name),
    seed= 32,
    color = 'black',
    size = 4,
    box.padding = unit(0.4, "lines"),
    point.padding = unit(0.4, "lines"),
    segment.curvature = -0.1,
    segment.ncp = 4,
    segment.angle = 20,
    min.segment.length = 0,
    nudge_x = 0.07
    #nudge_y = 0.07
  )+
  theme(
    legend.position = 'none',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size =14),
    axis.text.x = element_text(hjust = 0.7),
    strip.text = element_text(size = 16), 
     #strip.text.y = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )+
  guides(color =  guide_legend(title = 'Prey', ncol = 1))

# Significant correlation observed 
# r = 0.89 p=0.012 optMyo3 free vs optMyo3-Myo3p
# r = 0.64 p=0.14 optMyo5 free vs optMyo5-Myo5p
#ggsave('~/Maitrise/PCA_DHFR_array/PCA_SH3_libre/Figure_ancSH3_paper/effect_dmotifvspara.png', 
 #      height = 4.5, width = 6)



#heatmap ppi score for mat supp
compare_means(PPI_score ~ sh3_sequence, data = difdata[difdata$sh3_sequence != 'emptyLP' &
                                                      difdata$DHFR12_tag == 'C-tag', ], 
              p.adjust.method = 'BH', method = 'kruskal.test')
# difference non-significant between SH3s

FigSuppD <- 
ggplot(difdata) +
  facet_grid(
    cols = vars(DHFR12_tag),
    scales = 'free',
    drop = T,
    as.table = F,
    space = 'free_x'
  ) +
  geom_tile(aes(
    x = sh3_sequence,
    y = Prey.Standard_name ,
    fill = (PPI_score)
  )) +
  scale_x_discrete(limits = c('emptyLP', 'AncC', 'AncB', 'AncA', 'DupSH3', 'optMyo3', 'optMyo5'))+
  scale_fill_viridis(option = 'mako')+
  theme_bw() +
  theme(
    legend.position = 'right',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size =14),
    axis.text.x = element_text(angle = 40, hjust = 1),
    strip.text = element_text(size = 16),
    #legend.key.width = unit(1, "cm"), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  ) +
  labs(fill = 'med.\nPPI score', x = 'free SH3', y = 'Prey')



# Figure 5C :  PCA score comparison between optSH3s and AncC


figdata <- 
cdata[cdata$sh3_sequence %in% c('AncC', 'optMyo3', 'optMyo5') &
                cdata$DHFR12_tag== 'C-tag', c(1:6,19, 21)]

figdata[!figdata$motif_deletion.all, ] %>%
dplyr::filter(estradiol.all == estradiol.cor)-> figdata


figdata$sh3_sequence<- 
  factor(figdata$sh3_sequence, 
         levels = c( 'AncC', 'optMyo3', 'optMyo5'))

Fig5C <- 
ggplot(figdata, aes(Prey.Standard_name, PCA_score))+
  geom_boxplot(aes(fill = sh3_sequence), width = 0.5, alpha = 0.8)+
  scale_fill_manual(values= c('#551A8B', '#3366CC', 'orangered'))+
  stat_compare_means(aes(group = sh3_sequence), label = 'p.signif', label.y = 0.92)+
  ylim(0,1)+
  ylab('PPI score')+
  xlab('Prey')+
  theme_bw()+
  theme(
    legend.position = 'bottom',
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size =12),
    #axis.text.x = element_blank(),
    strip.text = element_text(size = 16), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )+
  guides(fill = guide_legend(title = 'free SH3', nrow = 2))


# Figure 5 assembly
Fig5A <- 
  ggdraw()+
  draw_image('~/ancSH3_paper/SupplementaryMaterial/FigurePanels/Fig5A.png')


fig5 <- 
plot_grid(Fig5A, Fig5B, Fig5C,labels = 'AUTO', label_fontface = 'plain', 
          label_size = 16, ncol = 3, rel_widths = c(0.8,1,1))

ggsave(fig5, file='~/ancSH3_paper/Figure5.png', width = 12.5, height = 4.5)

# Supplementary Figure 5 assembly

FigSuppA <- 
  readRDS(file = '~/ancSH3_paper/SupplementaryMaterial/FigurePanels/SuppFig5A.rds')


left <- 
plot_grid(Fig5SuppC+theme(legend.text = element_text(size = 14),
                           legend.title = element_text(size = 16),
                           axis.title = element_text(size = 14),
                           axis.text = element_text(size =12), 
                          panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
          FigSuppD+theme(legend.text = element_text(size = 14),
                          legend.title = element_text(size = 16),
                          axis.title = element_text(size = 14),
                          axis.text = element_text(size =12), 
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
          ncol =1, nrow = 2, labels = c( 'C', 'D'), 
          label_fontface = 'plain', label_size = 16,
          align = 'hv', axis = 'lr', greedy = T)

tleft <- 
plot_grid(FigSuppA+ylab('PPI score replicate 2')+
            xlab('PPI score replicate 1')+
            theme(legend.text = element_text(size = 14),
                           legend.title = element_text(size = 16),
                           axis.title = element_text(size = 14),
                           axis.text = element_text(size =12),
                          panel.grid.major = element_blank(), panel.grid.minor = element_blank()), 
          left, labels = c('A', ''), label_fontface = 'plain',
          nrow = 2, rel_heights = c(1,2), 
          align = 'h', axis = 'l')



plot_grid(tleft, Fig5SuppB+theme(legend.text = element_text(size = 14),
                                  legend.title = element_text(size = 16),
                                  axis.title = element_text(size = 14),
                                  axis.text = element_text(size =12)),
          ncol = 2, nrow = 1, labels = c('', 'B'), 
          label_fontface = 'plain', label_size = 16)


ggsave('~/ancSH3_paper/SupplementaryMaterial/SupplementaryFigure5.png', width = 18, height = 14)



