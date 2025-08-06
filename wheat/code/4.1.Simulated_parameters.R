rm(list = ls())

# Simulation study
# Each curve can be characterized by two parameters a,b if we assume a sigmoid behavior with saturation parameter 100%
#3-feb-2022
library(readxl)
library(tidyverse)
library(sommer)

set.seed(123)

BOKU.temp = read_csv('../data/clean/full.temp.csv') %>% filter(Env == 'BOKU')
source('1.0.Phenotype_functions.R')
source('1.1.Sigmoid_functions.R')


tmp = read_csv('../data/field/raw.FHB.Rdata.csv') %>% filter(WS.partner == "BOKU")
full.long.BOKU = data.frame(
  n.ob = rep(paste("n",0:6,sep=""),nrow(tmp)),
  Year = rep(tmp$Year, each = 7),
  Rep = rep(tmp$Rep, each = 7),
  GID = rep(tmp$GID, each = 7 ),
  DSA = c(t(as.matrix(tmp[,c(46,37:42)]))),
  Date = c(t(as.matrix(tmp[,c(20,31:36)]))),
  FHB = c(t(as.matrix(tmp[,c(45,25:30)])))
)

BOKU.temp$Date = as.character(BOKU.temp$Date)
tmp = left_join(full.long.BOKU,BOKU.temp[c(1,3)])
BOKU.FHB.GDD = na.omit(tmp[order(tmp$Year,tmp$Rep,tmp$GID,tmp$n.ob),c(2:6,8,7)])
rownames(BOKU.FHB.GDD) = 1:nrow(BOKU.FHB.GDD)

#GDD-GDD0
tmp = BOKU.FHB.GDD %>% 
  filter(DSA == 0) %>%
  .[,c(1:3,6)]
colnames(tmp)[4] = "GDD0"
tmp2 = left_join(tmp,BOKU.FHB.GDD)

####################################################
Year = c(2020,2021)
Rep = c(1,2)
GIDs = unique(BOKU.FHB.GDD$GID)

# estimate best a (position paramater) and b (growth rate)
as = c()
bs = c()

for (ye in Year){
  print(ye)
  for (re in Rep){
    print(re)
    for (GI in GIDs){
      tmp = BOKU.FHB.GDD %>%
        filter(Year == ye) %>%
        filter(Rep == re) %>%
        filter(GID == GI)
      
      X = (tmp$GDD - tmp$GDD[1])/100
      Y = tmp$FHB/100

      params = c(a = 0, b = 1)
      loss_function = function(args,y = Y, x = X){
        a = args[1]
        b = args[2]
        pred = sigmoid(x,a,b)
        sse = sum((pred-y)^2)
        return(sse)
      }
      fit = optim(par = params, fn = loss_function)
      as = c(as, fit$par[1])
      bs = c(bs, fit$par[2])
    }
  }
}

# merge parameters a and b with traits
tmp = read_csv('../data/clean/FHB.GDD.Traits.csv') %>% filter(Env %in% 'BOKU') %>% 
  dplyr::select(Year, Rep, GID, trait.AUDPC:trait.FHB_1)

WBC.df = data.frame(
  Year = rep(Year,each = length(Rep)*length(GIDs)),
  Rep = rep(rep(Rep,each=length(GIDs)),length(Year)),
  GID = rep(GIDs,length(Rep)*length(Year)),
  a = as,
  b = bs) %>% right_join(tmp)

# pintar cualquier obs y ver ajuste de las observaciones a las curvas "ab"
ind = 918
row = WBC.df[ind,]
ye = row$Year
re = row$Rep
GI = row$GID

tmp = BOKU.FHB.GDD %>%
  filter(Year == ye) %>%
  filter(Rep == re) %>%
  filter(GID == GI)

tmp2 = WBC.df %>%
  filter(Year == ye) %>%
  filter(Rep == re) %>%
  filter(GID == GI)

X = (tmp$GDD - tmp[1,]$GDD)/100 
Y = tmp$FHB/100
xg = seq(0,10,0.25)
yg = sigmoid(xg,tmp2$a, tmp2$b)
plot(xg,yg,"l",ylim = c(0,1))
points(X,Y)

#no correlation - we can simulate two independent distributions.
plot(WBC.df$a, WBC.df$b); cor(WBC.df$a, WBC.df$b)

#a may fit a normal
fit.a.20 = MASS::fitdistr(as[1:460],"normal")
fit.a.21 = MASS::fitdistr(as[461:920],"normal")
#b fits a gamma
fit.b.20 = MASS::fitdistr(bs[1:460],"gamma")
fit.b.21 = MASS::fitdistr(bs[461:920],"gamma")

# how GDD are distributed
BOKU.FHB.GDD %>%
  ggplot(aes(x=GDD,fill=DSA))+
  geom_histogram()+
  facet_grid(as.factor(DSA) ~ Year,scales = "free")

###################################################################
# simulation - marker effects for trait a and b (independent)
Geno230 = read_csv('../data/field/Geno.WheatSustain.csv')
GIDs = factor(Geno230$GID)
Geno230 = as.matrix(Geno230[,-1])
rownames(Geno230) = GIDs
n.marker = ncol(Geno230)
n.GID = nrow(Geno230)

a_marker.effects = matrix(rnorm(n.marker),ncol=1)
b_marker.effects = matrix(rnorm(n.marker),ncol=1)

a_genetic.values = scale(Geno230 %*% a_marker.effects)
b_genetic.values = scale(Geno230 %*% b_marker.effects)

#mapping genetic values to the true distributions before adding error
sim.as = c()
sim.bs = c()

for (ye in Year){
  for (re in Rep){
    print(paste(ye,re,sep="."))
    for (GI in GIDs){
      if (ye == 2020){
        new.a = qnorm(pnorm(a_genetic.values[GI,]),
                      mean = fit.a.20$estimate[["mean"]],
                      sd = fit.a.20$estimate[["sd"]])
        new.b = qgamma(pnorm(b_genetic.values[GI,]),
                       shape = fit.b.20$estimate[["shape"]],
                       rate = fit.b.20$estimate[["rate"]])
      }
      else{
        new.a = qnorm(pnorm(a_genetic.values[GI,]),
                      mean = fit.a.21$estimate[["mean"]],
                      sd = fit.a.21$estimate[["sd"]])
        new.b = qgamma(pnorm(b_genetic.values[GI,]),
                       shape = fit.b.21$estimate[["shape"]],
                       rate = fit.b.21$estimate[["rate"]])
      }
      sim.as = c(sim.as,new.a)
      sim.bs = c(sim.bs,new.b)
    }
  }
}

# #### Simulate anthesis variation (to set the position of day 0)
tmp = subset(BOKU.FHB.GDD, DSA == 0)
origin = as.Date("2020-01-01")
tmp$Date.num = as.numeric(as.Date(tmp$Date)-origin)
BOKU.temp$Date.num = as.numeric(as.Date(BOKU.temp$Date)-origin)

Anth.20.1 = MASS::fitdistr(tmp$Date.num[1:230],"normal")
Anth.20.2 = MASS::fitdistr(tmp$Date.num[231:460],"normal")
Anth.21.1 = MASS::fitdistr(tmp$Date.num[461:690],"normal")
Anth.21.2 = MASS::fitdistr(tmp$Date.num[691:920],"normal")

# #### residual variance vare depends on h2 and varg (var.a or var.b))
h2s = c(1,0.8, 0.5, 0.2)
var.a = var(sim.as)
var.b = var(sim.bs)
date.shift = c(0,seq(10,30,4))

##
set.seed(124)
simulate.df = WBC.df[,c(1:3)]
obs = c()
gdds = c()
for (h2 in h2s){
  #calculate ve from vg and h2
  var.ea = var.a*(1-h2)/h2
  var.eb = var.b*(1-h2)/h2
  #split into two distributions and add a different residual
  a_distr = sim.as + rnorm(4*n.GID,sd=sqrt(var.ea))
  b_distr = sim.bs + rnorm(4*n.GID,sd=sqrt(var.eb))
  integs = c()
  if (min(b_distr) < 0){
    b_distr = b_distr - min(b_distr)
  }
  i = 1 #auto inc
  for (ye in Year){
    for (re in Rep){
      print(paste(h2,ye,re,sep="."))
      for (GI in GIDs){
        if (ye == 2020){
          if(re == 1){
            anthesis = round(rnorm(1,Anth.20.1$estimate[["mean"]],Anth.20.1[["sd"]]+1))
          }
          else{
            anthesis = round(rnorm(1,Anth.20.2$estimate[["mean"]],Anth.20.2[["sd"]]+1))
          }
          gdd = suppressMessages(right_join(BOKU.temp,data.frame(Date.num = anthesis + date.shift)))$GDD
          obs = c(obs,sigmoid((gdd-min(gdd))/100,a_distr[i],b_distr[i]))
          gdds = c(gdds,gdd)
        }
        else{
          if(re == 1){
            anthesis = round(rnorm(1,Anth.21.1$estimate[["mean"]],Anth.21.1[["sd"]]+1))
            
          }
          else{
            anthesis = round(rnorm(1,Anth.21.2$estimate[["mean"]],Anth.21.2[["sd"]]+1))
            
          }
          gdd = suppressMessages(right_join(BOKU.temp,data.frame(Date.num = anthesis + date.shift)))$GDD
          obs = c(obs,sigmoid((gdd-min(gdd))/100,a_distr[i],b_distr[i]))
          gdds = c(gdds,gdd)
        }
        i = i + 1
      }
      simulate.df[,paste("a",h2,sep="_")] = a_distr
      simulate.df[,paste("b",h2,sep="_")] = b_distr
      #simulate.df[,paste("integral",h2,sep="_")] = integs
    }
  }
}

n.Obs = paste("n",0:6,sep="")
obs.df = data.frame(
  h2 = rep(h2s,each=length(Year)*length(Rep)*n.GID*length(n.Obs)),
  Year = rep(rep(Year,each=length(Rep)*n.GID*length(n.Obs)),length(h2s)),
  Rep = rep(rep(Rep,each=n.GID*length(n.Obs)),length(h2s)*length(Year)),
  GID = rep(rep(GIDs,each=length(n.Obs)),length(h2s)*length(Year)*length(Rep)),
  Obs = rep(n.Obs,length(h2s)*length(Year)*length(Rep)*n.GID),
  GDD = gdds,
  Value = obs)
tmp = subset(obs.df, Obs == "n0")[,c(1:4,6)]
colnames(tmp)[5] = "GDD0" 
tmp2 = left_join(obs.df, tmp)
obs.df$GDD0 = tmp2$GDD-tmp2$GDD0

#get true integrals
#combination of true a and b parameters
inte0 = obs.df %>% filter(Obs == "n0" & h2 == 1) %>% .[,"GDD"]
inte30 = obs.df %>% filter(Obs == "n6" & h2 == 1) %>% .[,"GDD"]
integs = c()
for (i in 1:nrow(simulate.df)){
  row = simulate.df[i,]
  B = primitive_sigmoid((inte30[i]-inte0[i])/100,row$a_1,row$b_1)
  A = primitive_sigmoid(0,row$a_1,row$b_1)
  integs = c(integs,B-A)
}

simulate.df$Integral = integs
simulate.df$GDD0 = inte0
simulate.df$GDD30 = inte30

obs.df %>%
  filter(h2 == 0.5) %>%
  ggplot(aes(x=GDD0,y=Value,group = interaction(Rep,GID)))+
  geom_line(alpha=0.05)+
  facet_grid(Rep ~ Year,scales="free")+
  theme_classic()

write_csv(simulate.df ,'../data/clean/simulated_parameters.csv')
write_csv(obs.df ,'../data/clean/simulated_observations.csv')


