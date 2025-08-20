rm(list = ls())

#First and second step with trial and common first step
#COVARIATES IN 2 STEP
library(tidyverse)
library(sommer)
library(lme4)

#Phenotype data (we need covariates)
trials.df = read_csv('../data/clean/FHB.GDD.Traits.csv')  #FIXED 18 APRIL
g.BLUEs = read_csv('../data/clean/g.BLUES.csv')
colnames(g.BLUEs)[2:5] = gsub('.*\\.','',colnames(g.BLUEs)[2:5])

#CV functions
source('CV_functions.R')

mean.no.nas = function(x){
  return(mean(na.omit(x)))
}


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
intersectnames <- intersect(g.BLUEs$GID,rownames(W))
W <- W[rownames(W)%in%intersectnames,]
 
g.covs = trials.df[,c("GID","cov.Height","cov.Anthesis")] %>%
  group_by(GID) %>% 
  summarise_each(funs(mean.no.nas))
g.input.BLUPs = as.data.frame(right_join(g.covs,g.BLUEs))

cors = c()
Reps = 25
Ys = c("AUDPC","Angle","GDD50","FHB_1")
covariates = c("cov.Height","cov.Anthesis")

Pheno = g.input.BLUPs 
Pheno = Pheno[Pheno$GID%in%intersectnames,]
Pheno$GID = factor(Pheno$GID, levels = rownames((W)))
#Artificial Env column to fit cv2.index() function
Pheno$Env = "World"
cv.object = cv2.index(Pheno,Reps = Reps)

for (ys in Ys){
  for (r in 1:length(cv.object)){
    print(paste(
      ys,
      r,
      sep="."))
    GEBVs.matrix = matrix(ncol=2**length(covariates))
    for (k in 1:length(cv.object[[r]]$folds)){
      sp = split.set(cv.object,"cv2",r,k)
      train.set = Pheno[sp$TRS.index,]
      test.set  = Pheno[sp$TS.index,]
      Ztest = model.matrix(~test.set$GID-1)
      
      model = mmer(as.formula(paste(ys,1,sep="~")),
                   random= ~ vs(GID, Gu=Kmatrix),
                   rcov= ~ units,
                   data=train.set, verbose = F,
                   date.warning = FALSE)
      GEBVs.submatrix = Ztest %*% matrix(model$U$`u:GID`[[ys]],ncol=1)
      for (m in 1:(2**length(covariates)-1)){
        hend.form = as.formula(paste(ys,
                                     paste(c(covariates[int2bin(m,
                                                                digits = length(covariates))]),
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
    rep_cors = apply(matrix(GEBVs.matrix,nrow=230),2,cor.with.phen)
    cors = c(cors,rep_cors)
  }
}

Models = paste("m",seq(length(covariates)^2)-1,sep="")
output.df = data.frame(
  Trait = rep(Ys,each=Reps*length(Models)),
  Rep = rep(rep(1:Reps,each = length(Models)),length(Ys)),
  Model = rep(Models,length(Ys)*Reps),
  Value = cors
)

#write_csv(output.df, file = "../output/scenario1.csv")

  
  
  
  

  