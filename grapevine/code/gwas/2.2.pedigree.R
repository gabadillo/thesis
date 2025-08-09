library(tidyverse)

Mmat = as.data.frame(read_csv('../../data/clean/Mmatrix.csv'))
M_grape = as.matrix(Mmat[,-1])
rownames(M_grape) = Mmat[,1]

Gmat = as.data.frame(read_csv('../../data/clean/Gmatrix.csv'))
G_grape = as.matrix(Gmat[,-1])
rownames(G_grape) = Gmat[,1]

rm(Mmat, Gmat)

table_ID = read_csv("../../data/raw/Grape_id_groups.csv")

big.BLUEs = read_csv('../../data/clean/preBLUEs.csv')
colnames(big.BLUEs)[1] = 'Geno'
big.BLUEs$Geno = gsub('gid', '', big.BLUEs$Geno)

colnames(table_ID)[1] = "Geno"
table_ID$Geno = paste0("G_",table_ID$Geno)

common_genos = intersect(big.BLUEs$Geno,colnames(G_grape))
sum(common_genos %in% table_ID$Geno) == length(common_genos)

geno_info = data.frame(Geno = common_genos) %>%
  right_join(table_ID) %>% filter(Geno %in% common_genos) %>% unique(.)

G_grape = G_grape[common_genos,common_genos]
M_grape = M_grape[common_genos,]

#PCA as EVD
dim(G_grape)
EVD = eigen(G_grape)
EVD.varexp = round(100*(EVD$values^2)[1:3]/sum(EVD$values^2),2)
df.EVD = cbind(geno_info,PC1=EVD$vectors[,1],PC2=EVD$vectors[,2])
df.EVD %>%
  ggplot(aes(x=PC1,y=PC2,color=group))+geom_point()+
  xlab(paste0("PC1 (",EVD.varexp[1],"%)"))+ylab(paste0("PC2 (",EVD.varexp[2],"%)"))

#PCA as SVD -> should be the same as prcomp
dim(M_grape)
X = scale(t(M_grape), center = T, scale=F)
SVD = svd(X)
SVD.varexp = round(100*(SVD$d^2)[1:3]/sum(SVD$d^2),2)

res.pca <- prcomp(X)
PRC.varexp = round(100*(res.pca$sdev^2)[1:3]/sum(res.pca$sdev^2),2)
PRC.varexp == SVD.varexp

SVD.G = svd(G_grape)
SVD.G.varexp = round(100*(SVD.G$d^2)[1:3]/sum(SVD.G$d^2),2)

#check that loadings are the same
cor(SVD$v[,1],res.pca$rotation[,1]);cor(SVD$v[,2],res.pca$rotation[,2])

df.SVD = cbind(geno_info,PC1=SVD$v[,1],PC2=SVD$v[,2],PC1E = SVD.G$v[,1],PC2E = SVD.G$v[,2])
df.SVD %>%
  ggplot(aes(x=PC1,y=PC2,color=group))+geom_point()+
  xlab(paste0("PC1 (",SVD.varexp[1],"%)"))+ylab(paste0("PC2 (",SVD.varexp[2],"%)"))

df.SVD %>%
  ggplot(aes(x=PC1E,y=PC2E,color=group))+geom_point()+
  xlab(paste0("PC1 (",SVD.G.varexp[1],"%)"))+ylab(paste0("PC2 (",SVD.G.varexp[2],"%)"))

# curate lines
# first adjustment: values below PC1< -.05 are self-pollinated (G_23 x G_23)
# G_23 itself will belong to jardin anyway

df.SVD$group[df.SVD$PC1E < - 0.05 & df.SVD$group!="G_23"] = "self"

centroids = df.SVD[,c("group","PC1E","PC2E")] %>%
  group_by(group) %>%
  summarise(across(everything(), mean))

dist2cluster = function(pc1,pc2,cx,cy){
  return(sqrt((pc2-cy)^2+(pc1-cx)^2))
}

new_groups = c()
dists = c()

for (i in 1:nrow(df.SVD)){
  subdists = c()
  for (j in 1:nrow(centroids)){
    subdists = c(subdists,dist2cluster(pc1=df.SVD[i,"PC1E"],pc2=df.SVD[i,"PC2E"],
                                       cx=centroids[j,"PC1E"],cy=centroids[j,"PC2E"]))
  }
  new_groups = c(new_groups, unname(which.min(subdists)))
  dists = c(dists,min(unlist(subdists)))
}

df.SVD$new_group = centroids$group[new_groups]
df.SVD$new_group[df.SVD$group == "jardin"] = "jardin"
df.SVD$dist = dists

#Figura 1
A = df.SVD %>%
  #filter(dist<1/50) %>%
  ggplot(aes(x=PC1E*EVD$values[1],y=PC2E*EVD$values[2],fill=new_group))+geom_point(alpha=1,shape=21,size=2)+
  xlab(paste0("PC1 (",SVD.G.varexp[1],"%)"))+ylab(paste0("PC2 (",SVD.G.varexp[2],"%)"))+
  theme_classic()+
  scale_fill_manual(values = c("#D64E12","#F9A52C","#EFDF48",
                               "#8BD346","#60DBE8","#16A4D8",
                               "#9B5FE0","pink","lightgray"))+
  theme(legend.position = "right",
        axis.title = element_text(face="bold"))+
  labs(fill=NULL)

A

# add family to preBLUEs to get BLUEs
big.BLUEs = right_join(big.BLUEs, df.SVD[,c('Geno', 'new_group')])
# keep jardin for G_23
big.BLUEs[1,'new_group'] = 'jardin'
big.BLUEs = big.BLUEs[,c(1,15,2:14)]
colnames(big.BLUEs)[2] = 'family'

write_csv(big.BLUEs, '../../data/clean/BLUEs.csv')

# save subset 588 x 49210 matrix (X) and subset G as kinship (K)
X = cbind(data.frame(Geno = rownames(M_grape)), M_grape)
rownames(X) = NULL
write_csv(X, file = '../../data/clean/X.csv')

K = cbind(data.frame(Geno = rownames(G_grape)), G_grape)
rownames(K) = NULL
write_csv(K, file = '../../data/clean/K.csv')

# save pedigree
pedigree = right_join(geno_info, df.SVD[,c('Geno', 'new_group')])
pedigree[1,'new_group'] = 'jardin'
pedigree$family = pedigree$new_group
pedigree$group = pedigree$new_group = NULL
write_csv(geno_info, file = '../../data/clean/pedigree.csv')

