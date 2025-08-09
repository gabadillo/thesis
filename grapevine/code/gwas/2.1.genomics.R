# 23-nov 2022 - after second meeting with Paola
# In this script:
# - load raw vcf (801 x 49210)
# - compute G matrix

library(tidyverse)
#numbers

grape_012 <- read_table("../../data/raw/grape_final.012", col_names = FALSE)
grape_012$X1 = NULL
grape_012 = as.matrix(grape_012)

#marker mapping and id
grape_markers = read_table("../../data/raw/grape_final.012.pos", col_names = FALSE)
grape_markers$X0 = paste("M", grape_markers$X1, grape_markers$X2, sep="_")
colnames(grape_markers) = c("CHROM","POS","ID")
grape_markers = grape_markers[,c(3,1,2)]

#geno names
grape_GIDs <- read_table("../../data/raw/grape_final.012.indv", col_names = FALSE)
grape_GIDs = paste0('G_',grape_GIDs$X1)

# prepare M
M_grape = grape_012
colnames(M_grape) = grape_markers$ID
rownames(M_grape) = grape_GIDs

Mmat = cbind(data.frame(GID = rownames(M_grape), M_grape))
rownames(Mmat) = NULL
write_csv(Mmat, file = '../../data/clean/Mmatrix.csv')

mapping = grape_markers

rm(grape_012, grape_markers, grape_GIDs)

# filter M -> done by humberto
# table(M_grape)
#        0        1        2
# 19773172 11178728  8465310

#fix names
mapping$CHROM[is.na(mapping$CHROM)] = "Un"
for (i in 1:nrow(mapping)){
  tmp = as.character(mapping[i,]$CHROM)
  if (nchar(tmp)<2){
    tmp = paste("0",tmp,sep="")
  }
  mapping[i,]$CHROM = tmp
}
mapping$ID = paste("M",mapping$CHROM,mapping$POS,sep="_")
mapping = as.data.frame(mapping) #no tibble pls

write_csv(mapping, file = '../../data/clean/mapping.csv')

#G matrix
G_grape = AGHmatrix::Gmatrix(M_grape, maf = 0.05)
Gmat = cbind(data.frame(GID = rownames(G_grape)), G_grape)
rownames(Gmat) = NULL
write_csv(Gmat, file = '../../data/clean/Gmatrix.csv')

