library(tidyverse)

pheno = as.data.frame(read_csv('../data/clean/pheno.csv'))
geno = as.data.frame(read_csv('../data/clean/geno.csv'))

ETA = list()

# Environment (E) -> K = I -> can be done with BRR
ZE = as.matrix(model.matrix(~as.factor(pheno$EID)-1))
ETA[['E']] = list(X = ZE, model = 'BRR')

# Line (L) -> K = I -> can be done with BRR
ZG = as.matrix(model.matrix(~as.factor(pheno$GID)-1))
ETA[['L']] = list(X = ZG, model = 'BRR')

# Genotype (g) -> K = MM'/c -> needs RKHS with EVD input (V and d instead of K)
X = as.matrix(geno[,-1])
S=0
for(i in 1:ncol(X)){
   meanXi <- mean(X[,i],na.rm=TRUE)
   X[,i]<-X[,i]-meanXi # centering
   X[,i]<-X[,i]/sd(X[,i]) # standarizing
   S<-S+var(X[,i])
}
G<-tcrossprod(X)/S
ZGZ = ZG%*%G%*%t(ZG)
EVD.G<-eigen(ZGZ)
ETA[['g']] = list(V = EVD.G$vectors, d = EVD.G$values, model = 'RKHS')

# Genotype x Environment (gE) -> K = ZGZ x ZEZ -> same RKHS rules
ZEZ = tcrossprod(ZE)
ZGEZ = ZGZ*ZEZ
EVD.GE<-eigen(ZGEZ)
ETA[['gE']] = list(V = EVD.GE$vectors, d = EVD.GE$values, model = 'RKHS')

save(ETA, file = '../data/clean/ETA.RData')
