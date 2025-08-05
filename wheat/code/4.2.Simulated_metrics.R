rm(list = ls())

library(tidyverse)

# compute metrics (AUDPC, ANGLE, GDD50, MAXVAR) from masked simulated values.

#############################################
#
obs.df = read_csv('../data/clean/simulated_observations.csv')
simulate.df = read_csv('../data/clean/simulated_parameters.csv')
# metrics functions
source('1.0.Phenotype_functions.R')

h2s = c(1,0.8, 0.5, 0.2)
Year = c(2020,2021)
Rep = c(1,2)
GIDs = unique(obs.df$GID)

Scenarios = list( # which observations are used to compute the metrics?
  all = 1:6,
  firsts = 1:3,
  mids = 3:5,
  evens = c(2,4,6),
  odds = c(1,3,5),
  limits = c(1,2,6),
  pair1 = c(2,5),
  pair2 = c(3,6),
  first = 1,
  fourth = 4,
  last = 6
)

########
phenos = c()
for (h in h2s){
  for (ye in Year){
    for (re in Rep){
      print(paste(h,ye,re,sep="."))
      for (GI in GIDs){
        full_obs = obs.df %>%
          filter(Year == ye) %>%
          filter(Rep == re) %>%
          filter(GID == GI) %>%
          filter(h2 == h) 
        for (sc in Scenarios){
          partial_obs = full_obs[c(1,sc+1),]
          #for AUDPC, Angle and GDD50
          xi = partial_obs$GDD0
          yi = partial_obs$Value*100 
          #FHB_1
          argmax_var = max_var_obs(sc)
          
          phenos = c(phenos, AUDPC(xi,yi,add0 = F),
                     GDD_angle.mid.point(xi,yi),
                     GDD_interpol_a(xi,yi,ymean = 50),
                     subset(partial_obs,Obs==argmax_var)$Value)
        }
      }
    }
  }
}

phenotypes = c("AUDPC","Angle","GDD50","FHB_1")
phenotype.df = data.frame(
  h2 = rep(h2s,each=length(Year)*length(Rep)*length(GIDs)*length(Scenarios)*length(phenotypes)),
  Year = rep(rep(Year,each=length(Rep)*length(GIDs)*length(Scenarios)*length(phenotypes)),
             length(h2s)),
  Rep = rep(rep(Rep,each = length(GIDs)*length(Scenarios)*length(phenotypes)),
            length(h2s)*length(Year)),
  GID = rep(rep(GIDs,each=length(Scenarios)*length(phenotypes)),
            length(h2s)*length(Year)*length(Rep)),
  Scenario = rep(rep(names(Scenarios),each=length(phenotypes)),
                 length(h2s)*length(Year)*length(Rep)*length(GIDs)),
  Phenotype = rep(phenotypes,
                  length(h2s)*length(Year)*length(Rep)*length(GIDs)*length(Scenarios)),
  Value = phenos)

write_csv(phenotype.df, file = '../data/clean/simulated_metrics.csv')
