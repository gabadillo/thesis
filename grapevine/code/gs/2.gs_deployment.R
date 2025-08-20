# GP models deployment
# - kernel methods: ABLUP, GBLUP, HBLUP, RKHS
# - bayesian and ML: BayesB, Bayesian LASSO, XGBoost

library(tidyverse)
library(sommer)
library(BGLR)
library(xgboost)

pheno = read_csv('../../data/gs/pheno.csv') %>% as.data.frame()
pheno$GID = factor(pheno$GID)
pheno = arrange(pheno)
all(pheno$GID == levels(pheno$GID)) # should be TRUE

A.matrix = read_csv('../../data/clean/A.csv')
A.matrix = as.matrix(A.matrix[,-1]); rownames(A.matrix) = colnames(A.matrix)
A.matrix = A.matrix[levels(pheno$GID), levels(pheno$GID)]

G.matrix = read_csv('../../data/clean/K.csv')
G.matrix = as.matrix(G.matrix[,-1]); rownames(G.matrix) = colnames(G.matrix)
G.matrix = G.matrix[pheno$GID, pheno$GID]

M.full = read_csv('../../data/clean/X.csv') %>% as.data.frame()
genonames = M.full[,1]
M.full = as.matrix(M.full[,-1])
rownames(M.full) = genonames
M.full = M.full[levels(pheno$GID),]

# load splits
load('../../data/gs/splits.RData')

# parameter box - this is controlled by an i index when it is parallelized in the cluster
scenario = 'D'
subset = 'all'
r = 1
trait = "B_weight"
# end parameter box 

# k-fold iteration
split = split_list[[scenario]][[subset]][[r]]
GEBVs.prematrix = list()
times = c(Sys.time())

for (k in 1:length(split$folds)){
  times = c(Sys.time())
  test.index = split$folds[[k]]
  train.index = sort(setdiff(split$full.index, test.index))
  
  GEBVs.submatrix = matrix(nrow=length(test.index),ncol=0)
  #Model 1 - ABLUP
  print("ABLUP")
  Z.train = model.matrix(~pheno[train.index,"GID"]-1)
  Z.test = model.matrix(~pheno[test.index,"GID"]-1)
  m1 = mmer(as.formula(paste(trait,1,sep="~")),
            random= ~ vsr(GID, Gu=A.matrix),
            rcov= ~ units,
            data=pheno[train.index,], verbose = F,dateWarning = F)
  GEBVs.submatrix = cbind(GEBVs.submatrix,
                          Z.test %*% matrix(m1$U$`u:GID`[[trait]],ncol=1))
  rm(m1)
  times = c(times,Sys.time())
  
  #Model 2 - GBLUP
  print("GBLUP")
  Z.train = model.matrix(~pheno[train.index,"GID"]-1)
  Z.test = model.matrix(~pheno[test.index,"GID"]-1)
  m2 = mmer(as.formula(paste(trait,1,sep="~")),
            random= ~ vsr(GID, Gu=G.matrix),
            rcov= ~ units,
            data=pheno[train.index,], verbose = F,dateWarning = F)
  GEBVs.submatrix = cbind(GEBVs.submatrix,
                          Z.test %*% matrix(m2$U$`u:GID`[[trait]],ncol=1))
  rm(m2)
  times = c(times,Sys.time())
  
  #Model 3 - HBLUP
  print("HBLUP")
  omega = c(0.1, 0.5, 1)
  
  HBLUP.hfun<-function(o){
    H.matrix = AGHmatrix::Hmatrix(A.matrix, G.matrix, omega = o)
    gblup = suppressMessages(mmer(as.formula(paste(trait,1,sep="~")),
                                  random= ~ vsr(GID, Gu=H.matrix),
                                  rcov= ~ units,
                                  data=pheno[train.index,], verbose = F,dateWarning=F))
    return(gblup$BIC)
  }
  
  HBLUP.models<-sapply(omega,HBLUP.hfun,simplify=TRUE)
  HBLUP.models
  o.best<-omega[which.min(HBLUP.models)]
  H.matrix = AGHmatrix::Hmatrix(A.matrix, G.matrix, omega = o.best)
  
  Z.train = model.matrix(~pheno[train.index,"GID"]-1)
  Z.test = model.matrix(~pheno[test.index,"GID"]-1)
  m3 = mmer(as.formula(paste(trait,1,sep="~")),
            random= ~ vsr(GID, Gu=H.matrix),
            rcov= ~ units,
            data=pheno[train.index,], verbose = F,dateWarning = F)
  GEBVs.submatrix = cbind(GEBVs.submatrix,
                          Z.test %*% matrix(m3$U$`u:GID`[[trait]],ncol=1))
  rm(m3)
  times = c(times,Sys.time())
  
  #M4 - RKHS
  print("RKHS")
  Eucl.distance = as.matrix(dist(M.full))
  nmarkers<-ncol(M.full)
  hvec<-rep(1/nmarkers,9)*c(1/50,1/10,1/3,1/2,1,2,3,10,50)
  
  RHKS.hfun<-function(h){
    Amat.RKHS<-exp(-h*(Eucl.distance)^2)
    gblup = suppressMessages(mmer(as.formula(paste(trait,1,sep="~")),
                                  random= ~ vsr(GID, Gu=Amat.RKHS),
                                  rcov= ~ units,
                                  data=pheno[train.index,], verbose = F,dateWarning=F))
    return(gblup$BIC)
  }
  times = c(times,Sys.time())
  
  RKHS.models<-sapply(hvec,RHKS.hfun,simplify=TRUE)
  RKHS.models
  h.best<-hvec[which.min(RKHS.models)]
  
  Amat.RKHS<-exp(-h.best*(Eucl.distance)^2)
  
  Z.train = model.matrix(~pheno[train.index,"GID"]-1)
  Z.test = model.matrix(~pheno[test.index,"GID"]-1)
  m4 = mmer(as.formula(paste(trait,1,sep="~")),
            random= ~ vsr(GID, Gu=Amat.RKHS),
            rcov= ~ units,
            data=pheno[train.index,], verbose = F,dateWarning = F)
  GEBVs.submatrix = cbind(GEBVs.submatrix,
                          Z.test %*% matrix(m4$U$`u:GID`[[trait]],ncol=1))
  rm(m4)
  times = c(times,Sys.time())
  
  #BGLR models
  dir.create('tmp')
  Bayes_nIter = 5000
  Bayes_burnIn = 1000
  Bayes_thin = 5
  
  # M5 - BayesB
  print("bayesB")
  m5<-BGLR(pheno[train.index,trait],
           ETA=list(list(X=M.full[pheno[train.index,"GID"],],
                         model="BayesB")),
           saveAt = "tmp/",
           nIter=Bayes_nIter,
           burnIn=Bayes_burnIn,
           thin = Bayes_thin,
           verbose=TRUE)
  GEBVs.submatrix <- cbind(GEBVs.submatrix,
                           M.full[pheno[test.index,"GID"],]%*%m5$ETA[[1]]$b)
  rm(m5)
  times = c(times,Sys.time())
  
  # M6 - Bayesian LASSO
  print("bayesL")
  m6<-BGLR(pheno[train.index,trait],
           ETA=list(list(X=M.full[pheno[train.index,"GID"],],
                         model="BL")),
           saveAt = "tmp/",
           nIter=Bayes_nIter,
           burnIn=Bayes_burnIn,
           thin = Bayes_thin,
           verbose=TRUE)
  GEBVs.submatrix <- cbind(GEBVs.submatrix,
                           M.full[pheno[test.index,"GID"],]%*%m6$ETA[[1]]$b)
  rm(m6)
  times = c(times,Sys.time())
  unlink('tmp', recursive = TRUE)
  
  # #Model 7 - XGBoost
  print("XGBoost")
  # xgboost does not allow missing values - impute means
  y_xgb =  pheno[train.index,trait]
  y_xgb[is.na(y_xgb)] = mean(y_xgb, na.rm = TRUE)
  m7 =  xgboost(
    data = M.full[pheno[train.index,"GID"],],
    label = y_xgb,
    nrounds = 500,
    early_stopping_rounds = 20,
    objective = "reg:squarederror",
    verbose = 0,
    nthread = 1
  )
  #end models
  
  GEBVs.submatrix = cbind(GEBVs.submatrix,
                          predict(m7, M.full[pheno[test.index,"GID"],]))
  
  GEBVs.prematrix[[k]] = GEBVs.submatrix
}

GEBVs.matrix = do.call('rbind', GEBVs.prematrix)

cors = apply(GEBVs.matrix, 2, cor, y = pheno[split$full.index,trait], use = 'complete.obs')
models = c('ABLUP', 'GBLUP', 'HBLUP', 'RKHS', 'BayesB', 'BayesL', 'XGBoost')

output = data.frame(Trait = trait, Scenario = scenario, Subset = subset, Rep = r,
                    Model = models, Cor = cors)
#output
#write_csv(output, sprintf('../../output/example_%s.csv', trait))

