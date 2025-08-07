library(tidyverse)
library(sommer)

# compute metrics (AUDPC, ANGLE, GDD50, MAXVAR) from masked simulated values.

#############################################
#

# 30 replicates using slurm
i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(i)

simulate.df = read_csv('../data/clean/simulated_parameters.csv')
phenotype.df = read_csv('../data/clean/simulated_metrics.csv')

# metrics functions
source('1.0.Phenotype_functions.R')
# CV functions
source('CV_functions.R')

##### Cross Validation - CLUSTER
#Task 4: Genomic predictions using GBLUP and Cross Validation scheme.
#Genomic relationship matrix
Geno230 = read_csv('../data/field/Geno.WheatSustain.csv')
genonames = Geno230$GID
Geno230 = as.matrix(Geno230[,-1])
GENO <- rrBLUP::A.mat(Geno230,impute.method = "EM",return.imputed = T)
Kmatrix <- GENO$A
rownames(Kmatrix) = colnames(Kmatrix) = genonames
W <- GENO$imputed
rownames(W) = genonames
intersectnames <- intersect(simulate.df$GID,rownames(W))
W <- W[rownames(W)%in%intersectnames,]

#Set number of repetitions - 1 because parallelization happens outside R
Reps = 1
phen.cors = c() #cor(GEBV, phenotpype)
inte.cors = c() #cor(GEBV, integral)
#integral is the value obtained from substract primitive(GDD30, aikr, bikr) - primitive(GDD0, aikr, bikr)
#where GDD30 and GDD0 are the growing degree days for last measurement day and anthesis, respectively,
#aikr is the true breeding value for parameter a for the i-th line, k-th year and r-th rep and so on
#for b.
#The integral is a way to collapse both independent parameters a and b into one.
a.cors = c() #cor (GEBV, a_1)  #these a_1 and b_1 are the 2 true breeding values!
b.cors = c() #cor (GEBV, b_1)

# factors
h2s = c(1,0.8, 0.5, 0.2)
phenotypes = unique(phenotype.df$Phenotype)

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

for (sc in names(Scenarios)){
  for (h in 1:length(h2s)){
    for (ys in phenotypes){
      Pheno = phenotype.df %>%
        filter(Scenario == sc) %>%
        filter(Phenotype == ys) %>%
        filter(h2 == h2s[h]) %>%
        .[c("GID","Value","Year","Rep")]
      Pheno = Pheno[Pheno$GID%in%intersectnames,]
      Pheno$GID = factor(Pheno$GID, levels = rownames((W)))
      
      #########
      Pheno$Env = "World"
      cv.object = cv2.index(data = Pheno, Reps = Reps)
      for (r in 1:length(cv.object)){
        GEBVs.vector = c()
        for (k in 1:length(cv.object[[r]]$folds)){
          print(paste(sc,h2s[h],ys,r,k,sep="."))
          sp = split.set(cv.object,"cv2",r,k)
          train.set = Pheno[sp$TRS.index,]
          test.set  = Pheno[sp$TS.index,]
          Ztest = model.matrix(~test.set$GID-1)
          
          model = mmer(Value ~ Year + Year:Rep,
                       random= ~ vs(GID, Gu=Kmatrix),
                       rcov= ~ units,
                       data=train.set, verbose = F)
          
          GEBVs.subvector = as.numeric(Ztest %*% matrix(model$U$`u:GID`$Value,ncol=1))
          GEBVs.vector = c(GEBVs.vector, GEBVs.subvector)
        }
        cv.index = cv.object[[r]]$full.index
        
        phen.cor = cor(GEBVs.vector,Pheno[cv.index,]$Value)
        inte.cor = cor(GEBVs.vector,simulate.df[cv.index,]$Integral)
        a.cor    = cor(GEBVs.vector,simulate.df[cv.index,]$a_1)
        b.cor    = cor(GEBVs.vector,simulate.df[cv.index,]$b_1)
        phen.cors = c(phen.cors,phen.cor)
        inte.cors = c(inte.cors,inte.cor)
        a.cors    = c(a.cors,a.cor)
        b.cors    = c(b.cors,b.cor)
      }
    }
  }
}


output.df = data.frame(
  Available = rep(names(Scenarios),each = length(h2s)*length(phenotypes)*Reps),
  h2 = rep(
    rep(h2s, each = length(phenotypes)*Reps),
    length(Scenarios)),
  Y.approach = rep(
    rep(phenotypes, each = Reps),
    length(Scenarios)*length(h2s)),
  Reps = i,
  phen.Value = phen.cors,
  inte.Value = inte.cors,
  a.Value = a.cors,
  b.Value = b.cors)

write_csv(output.df, file = sprintf('../data/output/cluster/rep_%s.csv', i))
