#############################################
## proline motif confirmation PCA analysis ##
## author : Pascale Lemieux                ##
## Date : 2023-01-10                       ##
#############################################

library(tidyverse)
library(magrittr)
library(gtools)
library(mixtools)

##############################################
## Import raw data generated by pyphe-quant ##
##############################################

#  For prey array plates, final diploid selection plates and final MTX plates #
# Read array plan in 1536 format
plate1536 <-
  read.csv2('~/ancSH3_paper/SupplementaryMaterial/Data/plate1536_motif_confirmation.csv')

# Set working directory with raw data from pyphe
setwd("~/ancSH3_paper/SupplementaryMaterial/Data/motif_confirmation_RD/MTX2/")
file = dir(pattern = ".csv")

# Assign names to each data frame object
j = 1
for (i in file) {
  assign(x = paste0('plate', j), read.csv(i))
  j = j + 1
  
}

# merge each data frame with array plan
for (i in 1:8) {
  
  p <- get(paste0('plate', i))
  assign('p',   
         merge(p, 
               plate1536, 
               by.x = c('row', 'column'), 
               by.y = c('row', 'col'), 
               all = T)
  )
  
  p$plate <- i
  assign(paste0('plate', i), p)
}

# assemble all data frame
# Due to incubation problems, the MTX2 plates 6 and 8 had to be removed
all_plate <-
  rbind(plate1, plate2, plate3, plate4, plate5, plate7)
# write file with all data
filename = '~/ancSH3_paper/SupplementaryMaterial/Data/motif_confirmation_data/MTX2_data.csv'# MTX or diploid selection
  write.csv(all_plate,
             file = filename,
             row.names = T,
             quote = T)

######################################
## Data verification  and filtering ##
######################################

#  Removal of non-growing position in the prey array and
#  diploid selection for the MTX data  #
#  Because the full-lenght paralog PCA experiment showed growth for 
#  all the prey tested here the DHFR3 array was not analysed

# Set working directory where all the output files will be generated
setwd('~/ancSH3_paper/SupplementaryMaterial/Data/motif_confirmation_data/')

# Read assembled data 
MTX2_data <- 
  read.csv(file = './MTX2_data.csv')[, -1]

diploid2_data <- 
  read.csv(file = './diploid2_data.csv')[, -1]


# Verification of no growth positions in diploid2 selection plates
# Histogram of colony area distribution of diploid selection
ggplot(data = diploid2_data)+
  geom_histogram(mapping = aes(x = area),
                 binwidth = 20)+
  xlim(0, 8000)+
  geom_vline(xintercept = 500)+
  theme_minimal()
# no low growth positions to remove


# Replace NAs in the MTX colony size by 0
MTX2_data$area <- 
  na.replace(MTX2_data$area, 0)

# Histrogram of the distribution of MTX colony area
ggplot(data = MTX2_data)+
  geom_histogram(mapping = aes(x = area),
                 binwidth = 20)+
  theme_minimal()

### Export clean MTX2 data
write.csv(x = MTX2_data, 
          file = '~/ancSH3_paper/SupplementaryMaterial/Data/motif_confirmation_data/MTX2_data_clean.csv', 
          quote = TRUE, 
          fileEncoding = 'UTF-8')


#############################
## MTX data transformation ##
#############################
library(viridis)

MTX2.data <- 
  read.csv(file = '~/ancSH3_paper/SupplementaryMaterial/Data/motif_confirmation_data/MTX2_data_clean.csv', 
           fileEncoding = 'UTF-8', )

# Columns removal and transformation
MTX2.data$X <- NULL
MTX2.data$X.1 <- NULL
MTX2.data$plate <- as.factor(MTX2.data$plate)

# Adjust area = 0 to area = 1 for logarithm transformation
MTX2.data %<>% mutate(adjust = area + 1)

# Transform colony adjusted area in log2 scale
MTX2.data$log2.area <- log2(MTX2.data$adjust)

#Distribution of log2(colony area) per plate
#Log2.distribution <- 
ggplot(data = MTX2.data) +
  geom_density(mapping = aes(x = log2.area, color = plate))+ 
  theme_minimal()+
  theme(legend.position = 'none')+
  theme_minimal()+
  geom_vline(xintercept = 3, color = 'orangered')+
  theme(legend.position = 'none')+
  xlim(0, 14)+
  theme_bw()+
  ylab('Density')+
  scale_color_viridis(option = 'G', discrete = T)+
  theme(legend.position = 'none')

# Remove the background
MTX2.data <-  subset(MTX2.data, 
                     MTX2.data$log2.area >= 3)

# Save the data frame object
saveRDS(MTX2.data, file = '~/ancSH3_paper/SupplementaryMaterial/Data/motif_confirmation_data/MTX2.data.rds')

###############################
## Normalization : PCA score ##
###############################
library(ggpubr)

# Important to compare the plates with each other
MTX2.list <- readRDS('~/ancSH3_paper/SupplementaryMaterial/Data/motif_confirmation_data/MTX2.data.rds')

#Withdrawal of the border

  MTX2.list %<>% filter(preys != 'border')

# Calculation of min/max for each plate
  MTX2.list %<>% 
    group_by(plate) %>% 
    mutate(max = max(log2.area, na.rm = TRUE))%>%
    ungroup()
  
  MTX2.list %<>% 
    group_by(plate) %>% 
    mutate(min = min(log2.area, na.rm = TRUE))%>%
    ungroup() 

# Calculation of standardized data using the min/max calculated 
 MTX2.list %<>% 
  mutate(norm = ((log2.area - min)/(max-min)))


# Verification of the correlation between technical replicates
# technical replicate plate combination : 1&5, 2&6, 3&7, 4&8

 dfrep1 = MTX2.list %>% filter(plate == 3)
 dfrep2 = MTX2.list %>% filter(plate == 7)
 
 dfrep1 %<>% mutate(baitXprey_1 = interaction(baits, preys)) %>% select(baitXprey_1, norm, log2.area)
 
 dfrep2 %<>% mutate(baitXprey_2 = interaction(baits, preys)) %>% select(baitXprey_2, norm, log2.area)
 
 crep <-  merge(dfrep1,
                dfrep2,
                by.x="baitXprey_1",
                by.y="baitXprey_2")
   
   ggplot(data = crep, aes(x = norm.x, y = norm.y)) +
     geom_point(aes(x = norm.x, y = norm.y), alpha = 0.2)+
     xlab('technical replicate 1')+
     ylab('technical replicate 2')+
     stat_cor()+
     theme_bw()+
     theme(legend.position = 'none')

# Export normalized data
saveRDS(MTX2.list, file = '~/ancSH3_paper/SupplementaryMaterial/Data/motif_confirmation_data/MTX2_norm.rds')   

###########################
## PPI score computation ##
###########################

MTX2.list <- readRDS(file = '~/ancSH3_paper/SupplementaryMaterial/Data/motif_confirmation_data/MTX2_norm.rds')
MTX2.list.rep <- vector(mode = 'list', length = 1)

# Compute the median of biological replicates (unique prey-bait combination)  
MTX2.list.rep <-  
  MTX2.list%>%
  group_by(baits, preys) %>%
  dplyr::summarise(n(), medsize = median(adjust ,na.rm = TRUE ), 
            medrep = median(norm, na.rm = TRUE ))

# Transform columns
MTX2.list.rep%<>% mutate(strain= str_c(baits, preys, sep = '.'))
MTX2.list.rep %<>% mutate(log2medsize = log2(medsize))

# Create object with PCA and PPI scores
MTX2.list.rep <-   
  merge(MTX2.list,
        MTX2.list.rep[, c('baits', 'preys', 'strain', 'medsize','log2medsize','medrep','n()')], 
        by = c('baits', 'preys'))

# Establishing the threshold for median colony area size 
# Bellow the threshold is noise 

x = 250 #colony size threshold

ggplot(data = MTX2.list.rep ) +
  geom_density(aes(x = medsize, color = plate)) +
  theme_minimal()+
  geom_vline(xintercept = x)+
  theme(legend.position = 'none')


# we try to identify PPI loss very low data points are conserved, 
# because we already know that PPI is detected for the WT preys

# Export the data (used for Supplementary Material, Figures)
saveRDS(MTX2.list.rep , file = '~/ancSH3_paper/SupplementaryMaterial/Data/motif_confirmation_data/MTX2data_rep.rds')



