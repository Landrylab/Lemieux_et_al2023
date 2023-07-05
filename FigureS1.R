######################################
## Supp Figure 1                    ##
## author : Pascale Lemieux         ##
## Date : 2023-01-19                ##
######################################
library(ggplot2)
library(tidyverse)
library(ggmsa)

# Read posterior probability for reconstructed amino acid at each node
anc_p <- 
read.table('~/ancSH3_paper/SupplementaryMaterial/SupplementaryData2/Ancestral_MaxMarginalProb_Char_Indel.txt', header = T)

# Keep nodes of interest (Dup, AncA, AncB, AncC)
Node <- c('N156', 'N149', 'N140', 'N105')
title <- factor(c('AncA', 'AncB', 'AncC', 'AncD'), levels = c('AncA', 'AncB', 'AncC', 'AncD'))

Anc <- 
  anc_p[anc_p$Node %in% c('N156', 'N149', 'N140', 'N105'),  ]
# remove empty positions in the reconstructed sequence
Anc <- 
  Anc[Anc$Char != '-', ]

ref <- 
data.frame(Node, title)

Anc <- 
merge(Anc, 
      ref,
      by = 'Node')

prob <- 
ggplot(Anc)+
  geom_histogram(aes(x = CharProb), binwidth = 0.05, color = 'black', fill = 'darkcyan')+
  facet_wrap(~title, nrow = 2)+
  theme_bw()+
  ylab('Frequency')+
  xlab('Probability')+
  theme( strip.text.x = element_text(size=16), 
         axis.text = element_text(size =12),
         legend.text = element_text(size = 12), 
         axis.title = element_text(size =14),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(prob, filename = '~/ancSH3_paper/SupplementaryMaterial/FigurePanels/SuppFig1C.png', height = 5, width = 5)

# Create MSA visualization for ancestral sequence
library(ggmsa)
library(seqinr)

sh3_sequences <-paste0("~/ancSH3_paper/SupplementaryMaterial/SupplementaryData2/myo3.myo5.sh3sequence.fasta")
protein_sequences <- system.file("extdata", "sample.fasta", package = "ggmsa")

# MSA ancestral sequence
p3 <- 
ggmsa(sh3_sequences, color = "Chemistry_AA", font = "DroidSansMono", char_width = 0.5, seq_name = TRUE )+
  geom_msaBar()

ggsave(p3, filename = '~/ancSH3_paper/SupplementaryMaterial/FigurePanels/SuppFig1E.png', width = 10, height = 2, device = 'png')


# MSA ortholog from the DupSH3 node
#sh3_sequences <-read.fasta("~/ancSH3_paper/SupplementaryMaterial/SupplementaryFiles1/SupplementaryFiles1_ortholog_sequences.fa")
#k <- read.table("~/ancSH3_paper/SupplementaryMaterial/SupplementaryFiles1/ortholog_dup_node.txt")

#write.fasta(sh3_sequences[names(sh3_sequences)%in% k$V1], names = names(sh3_sequences[names(sh3_sequences)%in% k$V1]),
#            file.out = '~/ancSH3_paper/SupplementaryMaterial/SupplementaryFiles1/ortholog_dup_MSA.fa')

sh3_sequences <- '~/ancSH3_paper/SupplementaryMaterial/SupplementaryData2/ortholog_dup_MSA.afa'
protein_sequences <- system.file("extdata", "sample.fasta", package = "ggmsa")


p1 <- 
ggmsa(sh3_sequences, font = "DroidSansMono", color = "Chemistry_AA",char_width = 0.5, seq_name = TRUE )+
  geom_msaBar()
ggsave(p1, filename = '~/ancSH3_paper/SupplementaryMaterial/FigurePanels/SuppFig1D.png', width = 12, height = 2, device = 'png')


# Assemble Supplementary Figure 1 manually
# Add secondary structure annotation on the MSA
# Add node label on the tree

