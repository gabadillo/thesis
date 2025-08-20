library(tidyverse)
library(BGLR)

pheno = as.data.frame(read_csv('../data/clean/pheno.csv'))
load('../data/clean/ETA.RData')
files = dir('../data/folds', full.names = TRUE)

traits = colnames(pheno)[-c(1:2)]

# BGLR params
nIter = 12000
burnIn = 2000
thin = 5

output = list()
dir.create('logs')

for (file in files){
  rep = str_match(file, '_(.+).csv')[,2]
  folds = as.data.frame(read_csv(file))
  for (i in 3:ncol(folds)){
    sc = colnames(folds)[i]
    print(sc)
    for (trait in traits){
      # set BGLR model
      cat(paste0(trait,':'))
      y = pheno[,trait]
      yNA = y
      yNA[folds[,i] != 2] = NA
      fm1 = BGLR(y = yNA, ETA = ETA[1:2], 
                 nIter = nIter, burnIn = burnIn, thin = thin, 
                 verbose = FALSE, saveAt = 'logs/fm1_')
      cat(' fm1')
      fm2 = BGLR(y = yNA, ETA = ETA[1:3], 
                 nIter = nIter, burnIn = burnIn, thin = thin, 
                 verbose = FALSE, saveAt = 'logs/fm2_')
      cat(' fm2')
      fm3 = BGLR(y = yNA, ETA = ETA[1:4], 
                 nIter = nIter, burnIn = burnIn, thin = thin, 
                 verbose = FALSE, saveAt = 'logs/fm3_')
      cat(' fm3\n')
      testSet = which(folds[,sc] == 1)
      output[[sprintf('%s_%s_%s', trait, sc, rep)]] = data.frame(Trait = trait, Rep = rep, Scheme = sc, 
                       EID = pheno$EID[testSet], GID = pheno$GID[testSet], yObs = y[testSet], 
                       yFm1 = fm1$yHat[testSet], yFm2 = fm2$yHat[testSet], yFm3 = fm3$yHat[testSet])
    }
  }
}

unlink('logs', recursive = TRUE)
output = bind_rows(output)
write_csv(output, '../output/predictions.csv')
