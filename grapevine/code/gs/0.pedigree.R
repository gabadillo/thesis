# Generate pedigree matrix for ABLUP
# compute A and H matrices

library(tidyverse)
library(AGHmatrix)

pedigree = read_csv('../../data/clean/pedigree.csv')

# adapt parents
pedigree$female_parent[!is.na(pedigree$female_parent)] =
  paste('G', pedigree$female_parent[!is.na(pedigree$female_parent)], sep = '_')
pedigree$male_parent[!is.na(pedigree$male_parent)] =
  paste('G', pedigree$male_parent[!is.na(pedigree$male_parent)], sep = '_')

# add known parents
new_pedigree = pedigree %>% dplyr::select(Geno, family, female_parent, male_parent) %>% 
  as.data.frame() # tibble strikes again
new_pedigree[new_pedigree$Geno=="G_Seleccion.5",3:4]   = c("G_Red","G_Dawn")
new_pedigree[new_pedigree$Geno=="G_23",3:4]            = c("G_Ruby","G_Centennial")
new_pedigree[new_pedigree$Geno=="G_Iniagrape.one",3:4] = c("G_Flame","G_Black")

# names to indices
input_Amatrix = new_pedigree
input_Amatrix[is.na(input_Amatrix)] = "0"
input_Amatrix$family = NULL

# metafounders
metafounders = c('G_Flame', 'G_Red', 'G_18V', 'G_Ruby', 'G_3V', 'G_5LL')
aux = data.frame(Geno = metafounders, female_parent = '0', male_parent = '0')

input_Amatrix2 = rbind(input_Amatrix, aux) %>% 
  arrange(female_parent, male_parent, Geno)

A.matrix = AGHmatrix::Amatrix(input_Amatrix2)
dim(A.matrix)

A.matrix = A.matrix[!rownames(A.matrix) %in% metafounders, !colnames(A.matrix) %in% metafounders]

K = read_csv('../../data/clean/K.csv')
G.matrix = as.matrix(K[,-1])
rownames(G.matrix) = K$Geno

H.matrix = Hmatrix(A=A.matrix, G=G.matrix, markers = marker.matrix, omega = 0)
dim(H.matrix)

genolevels = pedigree %>% arrange(family) %>% .$Geno
vizMatrix = function(m){
  m = reshape2::melt(m)
  m$Var1 = factor(m$Var1, levels = genolevels)
  m$Var2 = factor(m$Var2, levels = genolevels)
  ggplot(data = m, aes(x=Var1, y=Var2, fill=value))+
    geom_tile()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = 'none')+
    scale_fill_viridis_c()
}

plot_list = list()
plot_list[['A']] = vizMatrix(A.matrix)
plot_list[['G']] = vizMatrix(G.matrix)
plot_list[['H']] = vizMatrix(H.matrix)
ggpubr::ggarrange(plotlist = plot_list, nrow = 1, labels = names(plot_list))
boxplot(c(G.matrix) ~ c(A.matrix))

# mat to df
A = cbind(data.frame(Geno = rownames(A.matrix)), A.matrix)
H = cbind(data.frame(Geno = rownames(H.matrix)), H.matrix)

write_csv(new_pedigree, '../../data/clean/new_pedigree.csv')
write_csv(A, '../../data/clean/A.csv')
write_csv(H, '../../data/clean/H.csv')

