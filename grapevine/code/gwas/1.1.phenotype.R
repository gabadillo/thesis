# 24-nov 2022 after second meeting with Paola
# In this script:
# split traits in four groups:
# - SEED: seed_number, seed_fresh_weight and seed_dry_weight
# - BERRY: weight, height, witdh, shape (ratio)
# - HARVEST: cluster, rachis
# - POST: cluster, cluster_loss, rachis, rachis_loss

library(tidyverse)
library(readxl)

#load raw data
# segregantes
seg.pre <- read_excel("../../data/raw/phenotyping_plb_quantitative_raw.xlsx",
                      sheet = "espaldera_5")
seg.berry <- read_excel("../../data/raw/phenotyping_plb_quantitative_raw.xlsx",
                        sheet = "espaldera_5_berry_")
seg.post <- read_excel("../../data/raw/phenotyping_plb_quantitative_raw.xlsx",
                       sheet = "espaldera_5_post")

# jardin
jar.pre <- read_excel("../../data/raw/phenotyping_plb_quantitative_raw.xlsx",
                      sheet = "jvarieties")
jar.berry <- read_excel("../../data/raw/phenotyping_plb_quantitative_raw.xlsx",
                        sheet = "jvarieties_berry_size")
jar.post <- read_excel("../../data/raw/phenotyping_plb_quantitative_raw.xlsx",
                       sheet = "jvarieties_post")

##### harmonize segregants and jardin
# pre
seg.pre$plant = 1
jar.pre$family = "jardin"

seg.pre = seg.pre[,c(1:3,13,5,6:12)]
colnames(seg.pre) = c("season","family","gid","plant",
                      "replicated","H_cluster","B_weight",
                      "soluble_solids","S_number","S_fresh","S_dry","H_rachis")
jar.pre = jar.pre[,c(1,13,2,4,5,6:12)]
colnames(jar.pre) = c("season","family","gid","plant",
                      "replicated","H_cluster","B_weight",
                      "soluble_solids","S_number","S_fresh","S_dry","H_rachis")

full.pre = rbind(jar.pre, seg.pre)
full.pre$gid = paste("G",full.pre$gid,sep="_")
full.pre$season = as.numeric(full.pre$season)
rm(jar.pre, seg.pre)

# berry
jar.berry$family = "jardin"
seg.berry$plant = 1

jar.berry = jar.berry[,c(1,9,2,3,4,5:8)]
seg.berry = seg.berry[,c(1,2,3,9,4,5:8)]

colnames(jar.berry) = colnames(seg.berry) = c(
  "season","family","gid","plant",
  "replicated","B_height","B_width","S_dry","B_shape")
full.berry = rbind(jar.berry, seg.berry)
full.berry$gid = paste("G",full.berry$gid,sep="_")
full.berry$season = as.numeric(full.berry$season)
rm(jar.berry, seg.berry)

# post
seg.post$plant = 1
jar.post$family = "jardin"

seg.post = seg.post[,c(1,2,3,10,5,6:9)]
jar.post = jar.post[,c(1,10,2,4,5,6:9)]

colnames(jar.post) = colnames(seg.post) = c(
  "season","family","gid","plant","replicated",
  "H_cluster","P_cluster_weight","P_cluster_loss","P_rachis_weight"
)

full.post = rbind(jar.post, seg.post)
full.post$gid = paste("G",full.post$gid, sep = "_")
full.post$season = as.numeric(full.post$season)
rm(seg.post, jar.post)
#done

##### generate dataframes for each trait with appropiate covariables
#Rasgos de semilla: Sin covariable
#Rasgos de baya:  seed dry weight + SS si sirve (para ratio)
#Rasgos de racimo cosecha: SS + seed dry weight
#Rasgos de racimo a postcosecha: SS

phenoT = list()

### SEED - NO COVARIATES
# seed number
tmp = as.data.frame(full.pre[,c(1:5,9)])
phenoT[["S_number"]] = tmp

# seed fresh
tmp = as.data.frame(full.pre[,c(1:5,10)])
phenoT[["S_fresh"]] = tmp

# seed dry
tmp = as.data.frame(full.pre[,c(1:5,11)])
phenoT[["S_dry"]] = tmp

### BERRY - BOTH COVARIATES
# date is also added!
# berry weight #since they come from pre, covariates are already here
tmp = as.data.frame(full.pre[,c(1:5,7,8,11)])
phenoT[["B_weight"]] = tmp

# berry height
joinset = full.pre[,c(1:5,8,11)]
# get trait
tmp = as.data.frame(full.berry[,c(1:5,6,8)])
# add soluble solids and dry (dry were already joined in the original excel)
# i will use as a check
tmp2 = left_join(tmp,joinset)
tmp2 = unique(tmp2)
phenoT[["B_height"]] = tmp2

## berry width
tmp = as.data.frame(full.berry[,c(1:5,7:8)])
tmp2 = left_join(tmp,joinset)
tmp2 = tmp2[,c(1:6,8,7)]
tmp2 = unique(tmp2)
phenoT[["B_width"]] = tmp2

## berry shape
tmp = as.data.frame(full.berry[,c(1:5,8:9)])
tmp2 = left_join(tmp,joinset)
tmp2 = tmp2[,c(1:5,7,8,6)]
tmp2 = unique(tmp2)
phenoT[["B_shape"]] = tmp2

### HARVEST - BOTH COVARIATES
#since they come from pre, covariates are already here
## cluster
tmp = as.data.frame(full.pre[,c(1:5,6,8,11)])
phenoT[["H_cluster"]] = tmp

## rachis
tmp = as.data.frame(full.pre[,c(1:5,12,8,11)])
phenoT[["H_rachis"]] = tmp

### POST - only soluble solids
joinset = full.pre[,c(1:5,8)]
# here we have joined already H_cluster to compute H_cluster_loss
# i need to do the same for rachis
# loss formula = initial - final / initial
# post-cluster-weight
tmp = as.data.frame(full.post[,c(1:5,7)])
tmp2 = left_join(tmp,joinset)
tmp2 = unique(tmp2)
phenoT[["P_cluster_weight"]] = tmp2

# post-cluster-loss
tmp = as.data.frame(full.post[,c(1:5,8)])
tmp2 = left_join(tmp,joinset)
tmp2 = unique(tmp2)
phenoT[["P_cluster_loss"]] = tmp2

# post-rachis-weight
tmp = as.data.frame(full.post[,c(1:5,9)])
tmp2 = left_join(tmp,joinset)
tmp2 = unique(tmp2)
phenoT[["P_rachis_weight"]] = tmp2

# post-rachis-loss
# reuse tmp2
tmp3 = left_join(tmp2,phenoT[["H_rachis"]])
tmp3 = unique(tmp3)
tmp3$P_rachis_loss = 100*(tmp3$H_rachis-tmp3$P_rachis_weight)/tmp3$H_rachis

tmp3$P_rachis_weight = tmp3$H_rachis = NULL
tmp3 = tmp3[,c(1:5,8,6,7)]
# set negative to zero
tmp3$P_rachis_loss[tmp3$P_rachis_loss<0] = 0
phenoT[["P_rachis_loss"]] = tmp3

## DONE
for (phenotype in names(phenoT)){
  write_csv(phenoT[[phenotype]], file = sprintf('../../data/clean/phenotype/%s.csv', phenotype))
}

