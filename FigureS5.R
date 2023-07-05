# Prey orthologs analysis
# Figure S5
# author : Pascale Lemieux
# 20-06-2023



library(seqinr)
library(tidyverse)
library(httr)
library(jsonlite)
library(xml2)
library(ggplot2)


# SH3 dependent partners
yeast_id <- readRDS('~/ancSH3_paper/Reviews/SupplementaryMaterial/Data/sh3_dep_prey.rds')

#read ref species

setwd('~/ancSH3_paper/SupplementaryMaterial/SupplementaryData2/')

myo_ortho <- read.fasta('ortholog_sequences.fa')
species_ortho <- names(myo_ortho)

species_ortho_myo <- 
  gsub(names(myo_ortho), pattern = '/.{0,10}$', replacement = '')
  
d_species <- duplicated(species_ortho_myo)
species_ortho_myo <- gsub(species_ortho_myo, pattern = '.*_YEAST', replacement = 'YEAST')


# compara dataset june 2021
setwd('~/ancSH3_paper/Reviews/SupplementaryMaterial/SupplementaryData6/')

compara <- 
  read.table('Compara.103.protein_default.homologies.tsv', header = T)
stand_name <- read.csv('~/ancSH3_paper/SupplementaryMaterial/Data/Prey.Standard_names.csv')[, -1]


# read tree myo3 myo5
library(ggplot2)
library(ggtree)
library(treeio)
myotree <- 
  read.newick('~/ancSH3_paper/SupplementaryMaterial/SupplementaryData2/ortholog_phylogeny.newick')

# modify tree label to fit compara dataset
myotree$tip.label[grep('YEAST', myotree$tip.label)] <- 'saccharomyces_cerevisiae'
label <- vector(mode = 'character', length = length(myotree$tip.label))
specie <- unique(compara$homology_species)

for(i in specie){

  x <- grepl(i, myotree$tip.label)
  label[x] <- i

}

label[label == ''] <- gsub(myotree$tip.label[label == ''], pattern = '_[A-Za-z]+_[0-9]+_[0-9]+', replacement = '')


compara <- 
  compara[compara$gene_stable_id %in% yeast_id, ]

sub_compara <- 
 compara[compara$homology_species %in% label, ]

sub_compara$homology_type <- gsub(sub_compara$homology_type, pattern = 'ortholog_', replacement = '')

sub_compara <- 
  merge(sub_compara, 
        stand_name, 
        by.x = 'gene_stable_id', 
        by.y = 'Prey.Systematic_name')



library(tidytree)
library(TDbook)
library(tidyverse)
library(gtools)

x <- sub_compara[, c(1,5,8)]
x$seen <- 1
x_compara <- pivot_wider(x, names_from = homology_type, values_from = seen, values_fn = unique)

x_compara <- 
  na.replace(x_compara, replace = 0)

label[duplicated(label)] <- paste0(label[duplicated(label)], '_2')
myotree$tip.label <- label


# identify ancestral SH3 node
Node <- c('N156', 'N149', 'N140', 'N105')
title <- factor(c('AncA', 'AncB', 'AncC', 'AncD'), levels = c('AncA', 'AncB', 'AncC', 'AncD'))

myotree$node.label %in% Node

tree <- as.treedata(myotree)

x_compara %>%
group_by(homology_species)%>% 
  mutate(sum(one2one, one2many, within_species_paralog, other_paralog))%>% select(c(2,8)) %>% unique() -> plot_compara


colnames(plot_compara)<- c('label', 'trait')

y <- grepl(pattern = '_2$', myotree$tip.label)
todup <- 
  gsub(myotree$tip.label[y], pattern = '_2$', replacement = '')

to_dup <- 
cbind(paste0(unlist(plot_compara[plot_compara$label %in% todup, 'label']), '_2'), plot_compara[plot_compara$label %in% todup, 'trait'])

colnames(to_dup)[1] <- 'label'

plot_compara <- 
  rbind(plot_compara, to_dup)


myotree$tip.label[!(myotree$tip.label %in% plot_compara$label)]


y <- 
  full_join(myotree, plot_compara, by = 'label')

# Orthologs information with the tree
ptree <- 
ggtree(y, branch.length = 'none')+
  geom_tippoint(aes(color = trait))+
  scale_colour_viridis_c(option = 'H')+
  theme(legend.position = 'bottom', 
        legend.title = element_blank(), 
        plot.margin = margin(t = 0, r = 0, b = 3, l = 0, unit = "cm"))


ptree <- ggtree(y, branch.length = 'none') +
  geom_tippoint(aes(color = trait)) +
  scale_colour_viridis_c(option = 'H') +
  theme_tree2()+
  theme(legend.position = 'left',
        legend.title = element_blank())
ptree


y@phylo[["tip.label"]] %in% sub_compara$homology_species 

colnames(sub_compara)[8] <- 'label'
p_compara <- 
  sub_compara[, c(8,5,16)]


p_compara[-c(which(duplicated(p_compara))), ]

x <- 
  pivot_wider(p_compara[-c(which(duplicated(p_compara[,c(1,3)]))), ], names_from = Prey.Standard_name, values_from = homology_type, values_fill = NA, id_expand = T)

x <- as.data.frame(x)
rownames(x) <- x$label

x <- x %>% replace(.=="NULL", NA)
miss <- y@phylo[["tip.label"]][!(y@phylo[["tip.label"]] %in% x$label)]

for(i in miss){
  u <- gsub(pattern = '_2', '', i)
 z <- x[x$label == u, ]
 if(nrow(z) == 0){
   next
 }
 z$label <- i
 rownames(z) <- i
 x <- rbind(x, z)
  
}


gheatmap(ptree, x[, -1], offset=0.5, width=0.6, font.size=18, 
         colnames_angle=90, hjust=1)+
  scale_fill_manual(values = c('#422CB2' , '#2C85B2', '#85B22C','#B22C2C', '#B26F2C'), 
                    na.value = 'black')+
  scale_x_ggtree()+
  theme(legend.position = 'bottom', 
        #legend.margin=margin(-10, 0, 0, 0),
        legend.text = element_text(size = 10), 
         axis.title = element_text(size =16), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
        #plot.margin = margin(t = 0, r = 0, b = 3, l = 0, unit = "cm"))

library(svglite)
#ggsave('FigureS5A.svg', width = 10, height = 14)


# Motif conservation result
setwd('~/ancSH3_paper/Reviews/SupplementaryMaterial/SupplementaryData6/motif_conservation_results/')

# Vizualisation of motif prediction

#get file list
motif_file <- list.files()
df_all <- data.frame()

for (i in motif_file) {
  
 df <- read.csv(i)
 df$Prey.Standard_name <- df[1,1]

df_all <- rbind(df_all, df)   
}

yeast_gene <- c('Bzz1', 'Lsb3', 'Myo5', 'Osh2', 'Pkh2', 'Sla1', 'Ste20')

df_all[df_all$Species %in% yeast_gene, 'Species'] <- 'saccharomyces_cerevisiae'
df_all[df_all$Taxonomy %in% yeast_gene, 'Taxonomy'] <- 'Saccharomyces'



test <- 
  pivot_wider(df_all[, -c(2:4)], values_from = Max_MSS, names_from = Prey.Standard_name)

test <- 
  as.data.frame(test)

lab <- 
  myotree$tip.label[grepl(myotree$tip.label, pattern = '_2$')]


#here duplicate species name to fit with tip label

to_d <- unlist(strsplit(lab, split = '_2'))

to_d <- test[test$Species %in% to_d, ]

to_d$Species <- str_c(to_d$Species, '_2')

test <- rbind(test, to_d)


rownames(test) <- test$Species


ptree <- ggtree(y, branch.length = 'none') +
  geom_tippoint(aes(color = trait)) +
  #scale_colour_viridis_c(option = 'H') +
  theme_tree2()+
  geom_nodelab()+
  #geom_highlight(node = 105)+
  #xlim(-.1, 10)+
  theme(legend.position = 'left',
        legend.title = element_blank())
#plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "cm"))
ptree

for(i in 4:ncol(test)){
  
  test[, i] <- as.numeric(unlist(test[, i]))
  
} 

x <- gheatmap(ptree, test[,-1], offset=0.5, width=0.6, font.size=18, 
         colnames_angle=90, hjust=1)+
  scale_fill_viridis_c(option="G", name="continuous\nvalue", limits = c(0,1), breaks=c(0,0.5,1))+
  scale_x_ggtree()+
  theme(legend.position = 'bottom', 
        #legend.margin=margin(-10, 0, 0, 0),
        legend.text = element_text(size = 10), 
        axis.title = element_text(size =16), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#ggsave(plot = x, '~/ancSH3_paper/Reviews/SupplementaryMaterial/FigureS5B.svg', width = 8, height = 14)
#manualy assemble orthology for the preys and the motifs with Inkscape