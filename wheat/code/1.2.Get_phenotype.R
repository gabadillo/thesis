rm(list = ls())

library(readxl)
library(tidyverse)

#Calculate phenotypes:
#AUDPC
#Angle
#GDD50

#eval and sigmoid functions
source("1.0.Phenotype_functions.R")
source("1.1.Sigmoid_functions.R")

#GDD data
full.temp = read_csv('../data/clean/full.temp.csv')

# original observations
df.FHB = read_csv('../data/field/raw.FHB.Rdata.csv')


full.temp$str.Date = as.character(full.temp$Date)
audpcs = c()
angles = c()
gdd50s = c()

for (i in 1:nrow(df.FHB)){
  y = unlist(df.FHB[i,c(45,25:30)])
  x.date = df.FHB[i,c(20,31:36)]
  remove = unique(c(which(is.na(x.date)),which(is.na(y))))
  if (!length(remove)==0){
    x.date = x.date[-remove]
    y = y[-remove]
  }
  x.Date = data.frame(str.Date = t(x.date))
  
  
  if (nrow(x.Date)>1){
    tmp = full.temp %>% 
      filter(Env == df.FHB[i,]$WS.partner) %>% 
      right_join(x.Date,by = "str.Date") %>%
      na.omit() 
    x = tmp$GDD-min(tmp$GDD)
    if (length(x) != length(y)){
      audpcs = c(audpcs,NA)
      angles = c(angles,NA)
      gdd50s = c(gdd50s,NA)
    }
    else{
      audpcs = c(audpcs,AUDPC(x,y,add0=F))
      angles = c(angles,GDD_angle.mid.point(x,y))
      gdd50s = c(gdd50s,GDD_interpol_a(x,y,ymean = 50))
    }
  }
  else{
    audpcs = c(audpcs,NA)
    angles = c(angles,NA)
    gdd50s = c(gdd50s,NA)
  }
}

# max var per location
tmp = df.FHB %>% 
  group_by(Location) %>% 
  summarise(
    FHB1 = var(FHB1, na.rm = TRUE),
    FHB2 = var(FHB2, na.rm = TRUE),
    FHB3 = var(FHB3, na.rm = TRUE),
    FHB4 = var(FHB4, na.rm = TRUE),
    FHB5 = var(FHB5, na.rm = TRUE),
    FHB6 = var(FHB6, na.rm = TRUE),
    )

maxvars = rep(NA, nrow(df.FHB))
for (loc in unique(tmp$Location)){
  print(loc)
  indices = which(df.FHB$Location %in% loc)
  subdf = df.FHB[indices,]
  subtmp = tmp[tmp$Location %in% loc,]
  fhbmax = names(which.max(subtmp[-1]))
  maxvars[indices] = unlist(df.FHB[indices, fhbmax])
  print(maxvars[indices][1:5])
}

# feldkirchen 2020 does not have FHB5 so FHB3
indices = which(df.FHB$Location %in% 'Feldkirchen' & df.FHB$Year %in% 2020)
maxvars[indices] = df.FHB$FHB3[indices]

df.FHB$Trial = paste(df.FHB$Year, df.FHB$Location, sep = '.') 
FHB.GDD.Traits = df.FHB %>% dplyr::select(Year, WS.partner, Location, Trial, Rep, Row, Col, GID)
colnames(FHB.GDD.Traits)[2] = 'Env'

#fix 18 april 2022
#FHB.GDD.Traits 
FHB.GDD.Traits$trait.AUDPC = audpcs
FHB.GDD.Traits$trait.Angle = angles
FHB.GDD.Traits$trait.GDD50 = gdd50s
FHB.GDD.Traits$trait.FHB_1 = maxvars

FHB.GDD.Traits$cov.Height = df.FHB$Height
FHB.GDD.Traits$cov.Anthesis = df.FHB$Anthesis.Date
FHB.GDD.Traits$cov.AEAR = df.FHB$Anther.extrusion.anther.retention

write_csv(FHB.GDD.Traits, file = "../data/clean/FHB.GDD.Traits.csv")




