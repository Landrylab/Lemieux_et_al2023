# PPI Databases complementary analysis
# Figure S2 
# author : Pascale Lemieux
# 5-07-2023
library(ggplot2)
library(ggupset)
library(gtools)
library(tidyverse)

firstup <- function(x) {
  substring(x, 2) <- tolower(substring(x, 2))
  x
}


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

# Identify the comparison with the abundance control information
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
  scale_fill_manual(values = c('grey20','darkcyan', 'grey65'), labels = c('Abundance control', 'Detected', 'Undetected in this study'))+
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


# comparison Biogrid Myo3 vs Myo5
compM35 <- table(exp_collapsed[, 1:2])
compM35 <- as.matrix(compM35)

spe_M35 <-  compM35[(compM35[, 1] == 0 | compM35[, 2] == 0), ]

spe_M35 <- as.data.frame(spe_M35)    
table(spe_M35[, 2:3])

# 6 specific to Myo3, and 22 specific to Myo5
detect_spe <- 
merge(spe_M35, 
      exp_collapsed, 
      by = c('orf', 'DHFR12.strain'))

# all paralog specific interaction undetected in the PCA experiment


## validation with string database
## only score above 700 have high confidence
string <- read.delim('~/ancSH3_paper/Reviews/SupplementaryMaterial/Data/4932.protein.links.v11.5_string.txt', sep = ' ')

string <- lapply(string, gsub, pattern = '4932.', replacement = '', fixed = T)

string <- 
  as.data.frame(string)
# select only Myo3 and Myo5 PPI
sub_string <- 
unique(subset(string, 
       subset = (protein1 %in% c('YKL129C', 'YMR109W') | protein2 %in% c('YKL129C', 'YMR109W'))))

sub_string$combined_score <- 
  as.numeric(sub_string$combined_score)

ggplot(sub_string)+
  geom_density(aes(x = combined_score))

# reorganisation of the bait vs prey
sub_string$DHFR12.strain <- NA
sub_string[sub_string$protein1 %in% c('YKL129C', 'YMR109W'), 4] <- sub_string[sub_string$protein1 %in% c('YKL129C', 'YMR109W'), 1]
sub_string[sub_string$protein2 %in% c('YKL129C', 'YMR109W'), 4] <- sub_string[sub_string$protein2 %in% c('YKL129C', 'YMR109W'), 2]

sub_string$DHFR12.strain <- 
  gsub(sub_string$DHFR12.strain, pattern = 'YMR109W', replacement = 'Myo5')
sub_string$DHFR12.strain <- 
  gsub(sub_string$DHFR12.strain, pattern = 'YKL129C', replacement = 'Myo3')

sub_string$orf <- NA
sub_string[!(sub_string$protein1 %in% c('YKL129C', 'YMR109W')), 5] <- sub_string[!(sub_string$protein1 %in% c('YKL129C', 'YMR109W')), 1]
sub_string[!(sub_string$protein2 %in% c('YKL129C', 'YMR109W')), 5] <- sub_string[!(sub_string$protein2 %in% c('YKL129C', 'YMR109W')), 2]


sub_string$present_string <- T

sub_string <-  sub_string[, c(3:6)] %>%
                dplyr::group_by(DHFR12.strain, orf)%>% unique() 

# specific PPI for paralogs
string_spe <- as.matrix(table(sub_string[, 2:3]))
strin_spe<-  string_spe[, (string_spe[1,] == 0 | string_spe[2,] == 0)]

strin_spe <- 
  as.data.frame(strin_spe)

table(strin_spe[strin_spe$Freq!=0, 1])
#  specific to myo3 :83, specific to myo5:121

specific_s <- sub_string

# filter with tested PPI only
prey_tested <- unique(read.csv('~/ancSH3_paper/SupplementaryMaterial/Data/plate1536_complete_paralog.csv')[, 4])[-1]
sub_string <- sub_string[sub_string$orf %in% prey_tested, ]


exp_collapsed <- 
merge(sub_string, 
      exp_collapsed, 
      by = c('orf', 'DHFR12.strain'),
      all = T)


length(unique(exp_collapsed$orf))
sum(prey_tested %in% exp_collapsed[exp_collapsed$present_string, 'orf'])
# 146/167 prey tested present in string


#MINT database
library(stringr)
mint <- read.delim('~/ancSH3_paper/Reviews/SupplementaryMaterial/Data/species_yeastMINT.txt', header = F)

#keep only Myo3 and Myo5 PPIs
para <- c('P36006', 'Q04439')
para <- str_c('uniprotkb:', para)

para_mint <- mint[mint$V1 %in% para | mint$V2 %in% para, ]

# retrieve systematic orf name 
orf.6 <- 
  str_extract_all(pattern = 'uniprotkb:[Y][A-Z]+.+(locus name)', para_mint[, 6])

orf.6 <- lapply(orf.6, strsplit, split = ':')
orf.6 <- lapply(orf.6, unlist)
orf.6 <- lapply(orf.6, '[[', 2)
orf.6 <- lapply(orf.6, strsplit, split = '(', fixed = T)
orf.6 <- lapply(orf.6, unlist)
orf.6 <- lapply(orf.6, '[[', 1)

orf.5 <- 
  str_extract_all(pattern = 'uniprotkb:[Y][A-Z]+.+(locus name)', para_mint[, 5])

orf.5 <- lapply(orf.5, strsplit, split = ':')
orf.5 <- lapply(orf.5, unlist)
orf.5 <- lapply(orf.5, '[[', 2)

orf.5[unlist(lapply(orf.5, is.null))] <- 'none('

orf.5 <- lapply(orf.5, strsplit, split = '(', fixed = T)
orf.5 <- lapply(orf.5, unlist)
orf.5 <- lapply(orf.5, '[[', 1)

# verfication of orf systematic name
orf.6[nchar(orf.6) != 7]
orf.5[nchar(orf.5) != 7]

orf.5 <- lapply(orf.5, gsub, pattern = 'YSC84', replacement = 'YHR016C')
orf.6<- lapply(orf.6, gsub, pattern = 'YSC84', replacement = 'YHR016C')
orf.6 <- lapply(orf.6, gsub, pattern = 'YPK1', replacement = 'YKL126W')

# create PPIs dataframe
para_mint <- 
data.frame('p1' = unlist(orf.5),
           'p2' = unlist(orf.6),
           'mi_score' = as.numeric(gsub(para_mint$V15, pattern = 'intact-miscore:', replace = '')))

para_mint$DHFR12.strain <- NA
para_mint[para_mint$p1 %in% c('YKL129C', 'YMR109W'), 4] <- para_mint[para_mint$p1 %in% c('YKL129C', 'YMR109W'), 1]
para_mint[para_mint$p2 %in% c('YKL129C', 'YMR109W'), 4] <- para_mint[para_mint$p2 %in% c('YKL129C', 'YMR109W'), 2]

para_mint$DHFR12.strain <- 
  gsub(para_mint$DHFR12.strain, pattern = 'YMR109W', replacement = 'Myo5')
para_mint$DHFR12.strain <- 
  gsub(para_mint$DHFR12.strain, pattern = 'YKL129C', replacement = 'Myo3')


para_mint$orf <- NA
para_mint[!(para_mint$p1 %in% c('YKL129C', 'YMR109W')), 5] <- para_mint[!(para_mint$p1 %in% c('YKL129C', 'YMR109W')), 1]
para_mint[!(para_mint$p2 %in% c('YKL129C', 'YMR109W')), 5] <- para_mint[!(para_mint$p2 %in% c('YKL129C', 'YMR109W')), 2]

for(i in 1:nrow(para_mint)) {
  if (is.na(para_mint[i, 'orf'])) {
    if (para_mint[i, 'p1'] == para_mint[i, 'p2']) {
      para_mint[i, 'orf'] <- para_mint[i, 'p1']
    }
    else if (para_mint[i, 'DHFR12.strain'] == 'Myo3') {
      if (para_mint[i, 'p1'] == 'YKL129C') {
        para_mint[i, 'orf'] <- para_mint[i, 'p2']
      } else if (para_mint[i, 'p2'] == 'YKL129C') {
        para_mint[i, 'orf'] <- para_mint[i, 'p1']
      }  
    } else if (para_mint[i, 'DHFR12.strain'] == 'Myo3') {
      if (para_mint[i, 'p1'] == 'YMR109W') {
        para_mint[i, 'orf'] <- para_mint[i, 'p2']
      } else if(para_mint[i, 'p2'] == 'YMR109W') {
        para_mint[i, 'orf'] <- para_mint[i, 'p1']
    }
  }
}  
}

# keep unique interaction
para_mint <-  
  para_mint[, c(3:5)] %>%
  dplyr::group_by(DHFR12.strain, orf)%>% unique() 

# establish specific ppi
mint_spe <- as.matrix(table(para_mint[, 2:3]))
mint_spe<-  mint_spe[, (mint_spe[1,] == 0 | mint_spe[2,] == 0)]

mint_spe <- 
  as.data.frame(mint_spe)

table(mint_spe[mint_spe$Freq!=0, 1])
#  specific to myo3 :12, specific to myo5:43
para_mint$present_mint<- T


#add to summary data frame
exp_collapsed <- 
merge(exp_collapsed, 
      para_mint, 
      by = c('DHFR12.strain', 'orf'), 
      all = T)

exp_collapsed <- exp_collapsed[exp_collapsed$orf %in% prey_tested, ]

exp_collapsed$present_biogrid <-!is.na(exp_collapsed$recovered)
exp_collapsed[, c(4,8,9)] <- 
  na.replace(exp_collapsed[, c(4,8,9)], replace = F)
exp_collapsed$recovered <- 
  na.replace(exp_collapsed$recovered, replace = 'Undetected')
exp_collapsed <- 
  exp_collapsed[, c(1,2,5,6,9,4,3,8,7)]

colnames(exp_collapsed)[4] <- 'biogrid_method'
colnames(exp_collapsed)[c(7,9)] <- c('string_comb.score', 'mint_score')

#create summary column with DB present information (similar as biogrid method)

db <- exp_collapsed[, c(5,6,8)]

db$present_biogrid <- gsub(db$present_biogrid, pattern = TRUE, replacement ='BioGRID')
db$present_string <- gsub(db$present_string, pattern = TRUE, replacement ='STRING')
db$present_mint <- gsub(db$present_mint, pattern = TRUE, replacement ='MINT')

x <- apply(as.matrix(db), 1, str_c, simplify = F)

x <- lapply(x, str_remove_all, pattern = 'FALSE')

x <- 
  lapply(x, str_subset, pattern = '.+')

exp_collapsed$database <- x

# Summary detected in string
dose.labs <- c('STRING')
names(dose.labs) <- c(T)

string_p <- 
ggplot(exp_collapsed[exp_collapsed$present_string & exp_collapsed$recovered != 'Abundance control',])+
  facet_grid(rows = vars(recovered), cols = vars(present_string), 
             labeller = labeller(present_string = dose.labs))+
  geom_histogram(aes(x = string_comb.score, fill = recovered))+
  scale_fill_manual(values = c('darkcyan', 'grey65'))+
  theme_bw()+
  xlab('STRING score')+
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(), legend.position = 'none', 
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size =12), 
        strip.text = element_text(size =16),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides()

# MINT score distribution
dose.labs <- c('MINT')
names(dose.labs) <- c(T)
mint_p <- 
ggplot(exp_collapsed[exp_collapsed$present_mint & exp_collapsed$recovered != 'Abundance control',])+
  facet_grid(rows = vars(recovered), cols = vars(present_mint), 
            labeller = labeller(present_mint = dose.labs))+
  geom_histogram(aes(x = mint_score, fill =recovered), na.rm = T)+
  scale_fill_manual(values = c('darkcyan', 'grey65'))+
  theme_bw()+
  xlab('MINT score')+
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(), legend.position = 'none', 
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size =12), 
        strip.text = element_text(size =16), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides()

# Database comparison
db_comp <- 
  ggplot(exp_collapsed, aes(database, fill = recovered))+
  facet_grid(row = vars(DHFR12.strain))+
  geom_bar()+
  theme_bw()+
  scale_x_upset(n_intersections = 20)+
  ylab('Reported in PPI databases')+
  xlab(NULL)+
  scale_fill_manual(values = c('grey20','darkcyan', 'grey65'), 
                    labels = c('Abundance control', 'Detected', 'Undetected in this study'))+
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

x <- 
subset(exp_collapsed, 
       subset = present_biogrid == T & present_mint == T & present_string ==T)

table(x$recovered)

#83/135 detected in all 3 DB
#Assembly of Figure S2
library(cowplot)

b <- plot_grid(string_p, mint_p, ncol = 2, rel_widths = c(1, 1.05), labels = c('C', 'D'), label_fontface = 'plain')

t <-
  plot_grid( db_comp+guides(fill = guide_legend(nrow = 3)), 
              p4+theme( 
                     legend.position = 'none'),
                ncol = 2, rel_widths = c(1,1.2), labels = 'AUTO', label_fontface = 'plain' , axis = 'tb', align = 'h')


plot_grid(t, b, nrow = 2, labels = '', label_fontface = 'plain', rel_heights = c(1, 0.7))
ggsave('~/ancSH3_paper/Reviews/SupplementaryMaterial/FigureS2.png', height = 8.5, width = 10)

# verify specific PPI in each database

strin_spe$string <- T
mint_spe$mint <- T
detect_spe$bio_grid <- T

strin_spe <- strin_spe[strin_spe$Freq == 1, ]
mint_spe <- mint_spe[mint_spe$Freq == 1, ]

strin_spe$PPI <- str_c(strin_spe$DHFR12.strain, strin_spe$orf)
mint_spe$PPI <- str_c(mint_spe$DHFR12.strain, mint_spe$orf)
detect_spe$PPI <- str_c(detect_spe$DHFR12.strain, detect_spe$orf)

unique(c(strin_spe$PPI, mint_spe$PPI, detect_spe$PPI))

spe <- 
merge(detect_spe[, -c(3, 7)], 
      mint_spe[, -c(3,5)], 
      by = c('DHFR12.strain', 'orf'),
      all = T)

spe <- 
merge(spe, 
      strin_spe[, -c(3,5)], 
      by = c('DHFR12.strain', 'orf'),
      all = T)


table(spe$DHFR12.strain)
#Myo3 = 92, #Myo5 = 143

spe <- 
merge(spe, 
      specific_s[, c(1:3)], 
      by = c('orf', 'DHFR12.strain'), 
      all.x = T)

spe$tested <- (spe$orf %in% prey_tested)

ggplot(spe)+
  geom_histogram(aes(x = combined_score, color = tested), na.rm = T)

x <- spe[spe$combined_score > 700 & !is.na(spe$combined_score), ]
table(x$DHFR12.strain)

spe[spe$bio_grid & !is.na(spe$bio_grid) & spe$combined_score > 700, ]


x

