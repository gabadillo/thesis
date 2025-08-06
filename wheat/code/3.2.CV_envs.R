rm(list = ls())

#WITHIN ENVIRONMENT AND COVARIATES
#USING PARAMS AND COVARIATES BY TRIALS
library(tidyverse)
library(sommer)
library(lme4)

#Phenotype data (we need covariates)
trials.df = read_csv('../data/clean/FHB.GDD.Traits.csv')  #FIXED 18 APRIL
g.BLUEs = read_csv('../data/clean/g.BLUES.csv')

#CV functions
source('CV_functions.R')

mean.no.nas = function(x){
  return(mean(na.omit(x)))
}


############### COVARIATE PROBLEMS
## Env     Height   Anthesis   AEAR    TOTAL
## -----  -------- ----------  ------  -------
## BOKU     Yes       Yes        Yes      8
## NMBU     Yes       Yes        No       4
## Secobra  Yes       Yes        No       4
## SZD      No        Yes        No       2

#Combinations (AEAR, Height, Anthesis)
# m0: 000 -> no covariates   (ALL)
# m1: 001 -> Anthesis        (ALL)
# m2: 010 -> Height          (BOKU,NMBU,Secobra)
# m3: 011 -> Height Anthesis (BOKU,NMBU,Secobra)
# m4: 100 -> AEAR            (BOKU)
# m5: 101 -> AEAR Anthesis   (BOKU)
# m6: 110 -> AEAR Height     (BOKU)
# m7: 111 -> All covariates  (BOKU)

#binary function
int2bin = function(x,digits=3,radix=2){
  out = c()
  while (x!=0){
    out = c(x%%radix,out)
    x = x%/%2
  }
  if (digits > length(out)){
    out = c(rep(0,digits-length(out)),out)
  }
  return(which(out %in% 1))
}

#Genomic relationship matrix
Geno230 = read_csv('../data/field/Geno.WheatSustain.csv')
genonames = Geno230$GID
Geno230 = as.matrix(Geno230[,-1])
GENO <- rrBLUP::A.mat(Geno230,impute.method = "EM",return.imputed = T)
Kmatrix <- GENO$A
rownames(Kmatrix) = colnames(Kmatrix) = genonames
W <- GENO$imputed
rownames(W) = genonames
intersectnames <- intersect(trials.df$GID,rownames(W))
W <- W[rownames(W)%in%intersectnames,]

########################
#Cross Validation within environments
#using sommer:mmer
tmp = trials.df %>%
.[grep("TS",trials.df$GID),]
tmp = subset(tmp,Env != "SZD")

cors = c()
Reps = 10
Trials = unique(tmp$Trial)
Ys = c("AUDPC","Angle","GDD50","FHB_1")
colnames(tmp)[9:12] = Ys 
remove = which(is.na(tmp$AUDPC))
tmp = tmp[-remove,]
#tmp = na.omit(tmp)
n.models = list("BOKU"=4,"NMBU"=4,"Secobra"=4,"SZD"=2)
covariates = c("cov.Height","cov.Anthesis")

for (t in Trials){
  Pheno = tmp %>% filter(Trial == t)
  Pheno = Pheno[Pheno$GID%in%intersectnames,]
  Pheno$GID = factor(Pheno$GID, levels = rownames((W)))
  cv.object = cv1.index(Pheno, GID.name = 'GID', Reps = Reps)
  for (ys in Ys){
    for (r in 1:length(cv.object)){
      print(paste(t,
                  ys,
                  r,
                  sep="."))
      GEBVs.matrix = matrix(ncol = 4)
      for (k in 1:length(cv.object[[r]]$folds)){
        sp = split.set(cv.object,"cv1",r,k)
        train.set = Pheno[sp$TRS.index,]
        test.set  = Pheno[sp$TS.index,]
        #Ztrain = model.matrix(~train.set$GID-1)
        Ztest = model.matrix(~test.set$GID-1)
        
        model = mmer(as.formula(paste(ys,1,sep="~")),
                     random= ~ vs(GID, Gu=Kmatrix),
                     rcov= ~ units,
                     data=train.set, verbose = F,
                     date.warning = FALSE)
        GEBVs.submatrix = Ztest %*% matrix(model$U$`u:GID`[[ys]],ncol=1)
        for (m in 1:(4-1)){
          hend.form = as.formula(paste(ys,
                                       paste(c(covariates[int2bin(m,digits = 2)]),
                                             collapse="+"),
                                       sep="~"))
          model = mmer(hend.form,
                       random= ~ vs(GID, Gu=Kmatrix),
                       rcov= ~ units,
                       data=train.set, verbose = F,
                       date.warning = FALSE)
          GEBVs.submatrix = cbind(GEBVs.submatrix,
                                  Ztest %*% matrix(model$U$`u:GID`[[ys]],ncol=1))
        }
        GEBVs.matrix = rbind(GEBVs.matrix,GEBVs.submatrix)
      }
      
      GEBVs.matrix = GEBVs.matrix[-1,]
      cor.with.phen = function(x,y=Pheno[cv.object[[r]]$full.index,ys]){
        return(cor(x,y))
      }
      rep_cors = apply(matrix(GEBVs.matrix,nrow=nrow(GEBVs.matrix)),2,cor.with.phen)
      #rep_cors = c(rep_cors,rep(NA,max(unlist(n.models))-length(rep_cors)))
      cors = c(cors,rep_cors)
    }
  }
}

Models = paste("m",c(0:max(unlist(n.models)-1)),sep="")
#4 envs x 3 ys x 20 reps x 8 models = 1920
output.df = data.frame(Trial = rep(Trials,each=length(Ys)*Reps*length(Models)),
                       Trait = rep(rep(Ys,each=Reps*length(Models)),length(Trials)),
                       Rep = rep(rep(1:Reps,each = length(Models)),length(Trials)),
                       Model = rep(Models,length(Ys)*Reps*length(Trials)),
                       Value = cors
)

output.df$Value = c(cors, rep(NA,nrow(output.df)-length(cors)))
##################################################
#write_csv(output.df, file = '../data/output/scenario2.csv')



