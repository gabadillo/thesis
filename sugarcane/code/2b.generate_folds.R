library(tidyverse)

pheno = read_csv('../data/clean/pheno.csv')
pheno$EID = factor(pheno$EID)
t = length(levels(pheno$EID))
pheno$GID = factor(pheno$GID)
n = length(levels(pheno$GID))

# sanity check - should be TRUE
nrow(pheno) == n*t
SS = n/t

# splitting parameters - magic numbers
SSmin = 12 # should be <= than n/t
A2B = c(3, rep(4, 8)) # vector of genotypes transitioning from non-overlapping (A) to overlapping (B) groups 
seed = 1

# start algorithm
set.seed(seed)

# assign random initial genotypes' environment belonging - geno_shuffle and geno2env are IMMUTABLE
genos = levels(pheno$GID)
geno_shuffle = sample(genos)
geno2env = list()
for (i in 1:length(geno_shuffle)){
  imodt = (i-1)%%t+1
  if (!(i-1)%/%t){
    geno2env[[levels(pheno$EID)[imodt]]] = vector()
  }
  geno2env[[levels(pheno$EID)[imodt]]] = c(geno2env[[levels(pheno$EID)[imodt]]], geno_shuffle[i])
}

# sanity check - intersection between groups is 0
sum(do.call(c, geno2env) %in% geno_shuffle) == length(genos)

# start loop
j = 1 # outer index, controls sample size and transition from NOG to OG
# initial number of NOGs and OGs
A = SS
B = 0
# initial vectors for information storage 
sc = c()
SSlist = list()
output = list()
# pointers: cnt points to geno_shuffle for their removal from NOGs within a fixed sample size (inner loop)
cnt = 1
# pnt does something similar but across different sample sizes (outer loop)
pnt = 0
# dynamic storage of genotypes. 
NOGs = geno2env
OGs = c()
removed = c()

while (!SS < SSmin){
  for (i in j:length(A2B)){
    folds = rep(1, n*t)
    # update OGs and removed
    folds[pheno$GID %in% OGs] = 2
    folds[pheno$GID %in% removed] = 0
    # update NOGs
    for (env in levels(pheno$EID)){
      folds[pheno$EID %in% env & pheno$GID %in% NOGs[[env]]] = 2
    }
    if (length(OGs) > SS){break}
    # update A before saving
    A = SS - B
    newsc = sprintf('%s/%s',A,B)
    sc = c(sc, newsc)
    # save folds
    output[[newsc]] = folds
    # prepare new with more OGs
    B = B + A2B[i]
    for (k in 1:A2B[i]){ # for as many times as a genotype goes from NOG to OG...
      env = levels(pheno$EID)[(cnt-1)%%t+1] # round-robin selection
      OGs = c(OGs, NOGs[[env]][1]) # add the first NOG from the selected environment to OG
      NOGs[[env]] =  NOGs[[env]][-1] # remove that new OG from NOG
      for (notenv in levels(pheno$EID)[!levels(pheno$EID) %in% env]){
        # remove a NOG in each environment, as they are replaced by the new OG
        NOGs[[notenv]] = NOGs[[notenv]][-length(NOGs[[notenv]])] 
      }
      cnt = cnt + 1 # increase counter to keep the environment's round-robin selection
    }
  }
  # save scheme under 
  SSlist[[as.character(SS)]] = sc 
  # update SS by dropping OGs in the first scheme with non-zero B
  # the number of lines dropped is given by A2B[j], j updates every outer iteration
  SS = SS - A2B[j]
  
  # update NOGs by removing current lines in OGs (NOT CURRENT OGs - it contains everything now, but those in the first non-zero B)
  NOGs = geno2env
  
  # update pointer
  pnt = pnt + A2B[j]
  
  for (k in 1:pnt){
    kmodt = (k-1)%%t + 1
    check = unique(c(geno_shuffle[k], NOGs[[levels(pheno$EID)[kmodt]]][1]))
    if (length(check) != 1) stop('Error during OGs drop')
    NOGs[[levels(pheno$EID)[kmodt]]] = NOGs[[levels(pheno$EID)[kmodt]]][-1]
    removed = c(removed, check)
    for (notenv in levels(pheno$EID)[!levels(pheno$EID) %in% levels(pheno$EID)[kmodt]]){
      # remove a NOG in each environment, as they are replaced by the new OG
      NOGs[[notenv]] = NOGs[[notenv]][-length(NOGs[[notenv]])] 
    }
  }

  # reset vector of schemes
  sc = c()
  # resect OGs
  B = 0
  OGs = c()
  # update j
  j = j + 1
}

for (sc in names(output)){
  print(table(table(pheno$GID[output[[sc]] == 2])))
}

pheno = as.data.frame(pheno)
output = as.data.frame(do.call('cbind', output))
output = cbind(EID = pheno$EID, GID = pheno$GID, output)

# constant test?
all(apply(output[-c(1:2)], 2, function(x) sum(x==1)) == (t-1)*n)

write_csv(output, file = sprintf('../data/folds/rep_%s.csv', seed))
