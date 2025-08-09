
## libraries
require(tidyverse)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
#require(GAPIT3)

# LOAD DATA
big.BLUEs = as.data.frame(read_csv('../../data/clean/BLUEs.csv'))
G_grape = as.data.frame(read_csv('../../data/clean/K.csv'))
M_grape = as.data.frame(read_csv('../../data/clean/X.csv')) # tibble kills gapit
mapping = as.data.frame(read_csv('../../data/clean/mapping.csv'))


models = c("Blink") # 6/3 check berry weight
gapit_output <- list() # overwritten take a look!

for (j in 1:length(colnames(big.BLUEs)[-1])){ #
   s <- getwd()
   f <- paste0(getwd(), paste0("/../../output/raw_gwas/"))
   y2 <- colnames(data)[1+j]
   print(y2)

   y <- strsplit(y2,":")[[1]][1]
   dir.create(paste0(f,'/GAPIT_',models[i],'_',y,"/"), recursive = TRUE)
   setwd(paste0(f,'/GAPIT_',models[i],'_',y,"/"))
   # Y
   Y <- data.frame(data[c(1,(j+1))])
   Y = as.data.frame(na.omit(Y))
   # model
   mod <- GAPIT(Y=Y,
            GD=M_grape,
            GM=mapping,
            KI=G_grape,
            #CV=Q,
            PCA.total=6,
            model=models[i])

  gapit_output[[y2]] <- mod$GWAS
  setwd(s)
}

for (trait in names(gapit_output)){
   write_csv(gapit_output[[trait]], file = sprintf('../../data/gwas/%s.csv', trait))
}

