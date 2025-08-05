rm(list = ls())

#FIRST STEP
#gxe BLUEs - different equation for each environment, 690 BLUEs
#g BLUEs - using trial as factor, same equation for all, 230 BLUEs
#23-03-22

library(tidyverse)
library(sommer)
library(lme4)

#Phenotype data
FHB.GDD.Traits = read_csv('../data/clean/FHB.GDD.Traits.csv')

aux.df = FHB.GDD.Traits[grep("TS",FHB.GDD.Traits$GID),]
GIDs = factor(sort(unique(aux.df$GID)))

# genotype-in-environment BLUEs - needs to be done for every trait

BOKU.gxes = data.frame(Env = "BOKU", GID = GIDs)
NMBU.gxes = data.frame(Env = "NMBU", GID = GIDs)
Secobra.gxes = data.frame(Env = "Secobra", GID = GIDs)


traits = colnames(FHB.GDD.Traits)[9:12]
for (trait in traits){
  #BOKU 
  aux2.df = aux.df %>% filter(Env == "BOKU")
  BOKU.formula = as.formula(sprintf('%s ~ -1 + GID + Year + (1|Year:Rep:Row) + (1|Year:Rep:Col)', trait))
  BOKU.BLUEs = lmer(BOKU.formula, data = aux2.df)
  BOKU.gxes[,trait] = c(scale(BOKU.BLUEs@beta[1:230],scale=F))
  
  #NMBU
  aux2.df = aux.df %>% filter(Env == "NMBU")
  NMBU.formula = as.formula(sprintf('%s ~ -1 + GID + Year + (1|Year:Rep:Row) + (1|Year:Rep:Col)', trait))
  NMBU.BLUEs = lmer(NMBU.formula, data = aux2.df)
  NMBU.gxes[,trait] = c(scale(NMBU.BLUEs@beta[1:230],scale=F))
  
  #Secobra
  aux2.df = aux.df %>% filter(Env == "Secobra")
  Secobra.formula = as.formula(sprintf('%s ~ -1 + GID + Year', trait))
  Secobra.BLUEs = lm(Secobra.formula, data = aux2.df)
  Secobra.gxes[,trait] = c(scale(Secobra.BLUEs$coefficients[1:230],scale=F))
}

#############################
gxe.BLUEs = rbind(BOKU.gxes,
                  NMBU.gxes,
                  Secobra.gxes)

############################################
boxplot(gxe.BLUEs$trait.AUDPC ~ gxe.BLUEs$Env)
write_csv(gxe.BLUEs, file = '../data/clean/gxe.BLUES.csv')

#Single model for all environments (g.BLUEs)
g.BLUEs = data.frame(GID = GIDs)
for (trait in traits){
  NoEnv.BLUEs = lmer(as.formula(sprintf('%s ~ -1 + GID + Trial + (1|GID:Trial)', trait)) ,data = aux.df)
  g.BLUEs[,trait]= c(scale(NoEnv.BLUEs@beta[1:230],scale=F))
}

write_csv(g.BLUEs, file = '../data/clean/g.BLUES.csv')
