# Cross Validation with real data, covariates and Growing Degree Days
# raw data from original partner's datasets
# 21-03-2022

library(readxl)
library(tidyverse)


GDDs = function(temps, base = 5){
  output = c()
  acc = 0
  for (i in 1:length(temps)){
    if (is.na(temps[i])){
      temps[i] = 0
    } 
    GDD = temps[i]-base
    if (GDD > 0){
      acc = acc + GDD
    }
    output = c(output,acc)
  }
  return (output)
}

#### TASK 1: PREPARE DATASET
### BOKU
BOKU.temp <- read_excel("../data/weather/BOKU_Weather_Template.xlsx", 
                                    col_types = c("text", "numeric", "text", 
                                                  "date", "numeric", "text", "numeric", 
                                                  "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", 
                                                  "text", "text"))
BOKU.temp = BOKU.temp[,c(4,12)]
BOKU.temp$Env = "BOKU"
BOKU.temp$GDD = GDDs(BOKU.temp$`Tª mean`)

### NMBU
NMBU.temp <- read_excel("../data/weather/NMBU_Weather_2020_2021.xlsx", 
                                     col_types = c("text", "numeric", "text", 
                                                   "date", "numeric", "text", "numeric", 
                                                   "text", "numeric", "text", "numeric", 
                                                   "numeric", "numeric", "numeric", 
                                                   "text", "text", "numeric", "text", 
                                                   "numeric", "text", "numeric", "numeric", 
                                                   "numeric", "text", "text"))
NMBU.temp = NMBU.temp[,c(4,12)]
NMBU.temp$Env = "NMBU"
NMBU.temp$GDD = GDDs(NMBU.temp$`Tª mean`)

### Secobra
Secobra.temp <- read_excel("../data/weather/Secobra_Weather_Template.xlsx", 
                                       col_types = c("text", "text", "text", 
                                                     "date", "numeric", "text", "numeric", 
                                                     "text", "numeric", "text", "numeric", 
                                                     "numeric", "numeric", "numeric", 
                                                     "numeric", "text", "numeric", "numeric", 
                                                     "numeric", "text", "numeric", "numeric", 
                                                     "text", "text"))
Secobra.temp = Secobra.temp[,c(4,12)]
Secobra.temp$Env = "Secobra"
Secobra.temp$GDD = GDDs(Secobra.temp$`Tª mean`)

######
full.temp = rbind(BOKU.temp, NMBU.temp, Secobra.temp)
full.temp = full.temp[,c(1,3,4)]
full.temp = full.temp[order(full.temp$Date),]

write_csv(full.temp, file = "../data/clean/full.temp.csv")
