
library(tidyverse)

FDR<-function(pvals, FDR){
  pvalss<-sort(pvals, decreasing=F)
  m=length(pvalss)
  cutoffs<-((1:m)/m)*FDR
  logicvec<-pvalss<=cutoffs
  postrue<-which(logicvec)
  print(postrue)
  k<-max(c(postrue,0))
  cutoff<-(((0:m)/m)*FDR)[k+1]
  return(cutoff)}

files = dir('../../data/gwas/', full.names = TRUE)
clean.GWAS = lapply(files, read_csv) %>% lapply(as.data.frame)
names(clean.GWAS) = str_match(files, '^.*/(.+).csv')[,2]

SNP_table = list()
for (i in 1:length(clean.GWAS)){
  trait = names(clean.GWAS)[i]
  mod <- clean.GWAS[[i]][,1:5]
  colnames(mod)<-c("SNP","Chromosome","Position","P.value", "MAF")
  mod$LOD <- -log10(mod$P.value)
  fdr_value<-FDR(mod$P.value, 0.05)


  ## List of significant SNPs
  snpdf <- mod[mod$LOD>=-log10(fdr_value), ]
  if(dim(snpdf)[1]==0){
    snpdf[1,] <- NA
    snpdf$FDR <- -log10(fdr_value)
    snpdf$Bonferroni <- -log10(0.05/nrow(mod))
    snpdf$trait<-trait
  } else {
    snpdf$FDR <- -log10(fdr_value)
    snpdf$Bonferroni <- -log10(0.05/nrow(mod))
    snpdf$trait<-trait
  }
  #
  SNP_table[[trait]]<-snpdf
}

#save(SNP_table, file = paste0("./GWAS_output/associations.RData"))
assos = do.call("rbind",SNP_table)
assos = assos[order(assos$Chromosome,assos$Position),]
assos$bonfe =  assos$LOD > assos$Bonferroni
assos$Bonferroni = NULL
assos$P.value = NULL
assos = assos[,c('trait','SNP','Chromosome','Position','MAF','FDR','LOD','bonfe')]
colnames(assos) = c("Trait","Marker","Chr","Pos","MAF","FDR","LOD","Bonferroni")
assos = na.omit(assos)
assos = filter(assos, Chr != 20) %>%  arrange(Trait, Chr, Pos)
rownames(assos) = NULL
write_csv(assos, '../../output/gwas_hits.csv')
