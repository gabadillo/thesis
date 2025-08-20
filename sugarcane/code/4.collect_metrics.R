library(tidyverse)
source('Tiezzi_functions.R')

predictions = read_csv('../output/predictions.csv')

# correlation
envCorr = predictions %>% 
  group_by(Trait, Rep, Scheme, EID) %>% 
  summarize(n = n(),
    rFm1 = cor(yObs, yFm1),
    rFm2 = cor(yObs, yFm2),
    rFm3 = cor(yObs, yFm3)) %>% 
  group_by(Trait, Scheme, EID) %>% 
  summarize(n = mean(n),
            rFm1 = mean(rFm1),
            rFm2 = mean(rFm2),
            rFm3 = mean(rFm3))

allCorr = envCorr %>% 
  group_by(Trait, Scheme) %>% 
  summarize(
            rFm1 = Tiezzi(rFm1, n),
            rFm2 = Tiezzi(rFm2, n),
            rFm3 = Tiezzi(rFm3, n))

# RMSE
envRMSE = predictions %>% 
  group_by(Trait, Rep, Scheme, EID) %>% 
  summarize(n = n(),
            rmseFm1 = RMSE(yObs, yFm1),
            rmseFm2 = RMSE(yObs, yFm2),
            rmseFm3 = RMSE(yObs, yFm3)) %>% 
  group_by(Trait, Scheme, EID) %>% 
  summarize(n = mean(n),
            rmseFm1 = mean(rmseFm1),
            rmseFm2 = mean(rmseFm2),
            rmseFm3 = mean(rmseFm3))

allRMSE = envRMSE %>% 
  group_by(Trait, Scheme) %>% 
  summarize(
    rmseFm1 = mean(rmseFm1),
    rmseFm2 = mean(rmseFm2),
    rmseFm3 = mean(rmseFm3))

#write_csv(envCorr, '../output/envCorr.csv')
#write_csv(allCorr, '../output/allCorr.csv')
#write_csv(envRMSE, '../output/envRMSE.csv')
#write_csv(allRMSE, '../output/allRMSE.csv')
