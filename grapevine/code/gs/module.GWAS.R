# Module GWAS
# 11-01-23

# INPUT: phenotype, markers, mapping, kinship, trait, train-index
# OUTPUT: GWAS-LODs

library(GAPIT3)

GWAS = function(Phenotype.matrix,
                Marker.matrix,
                Kinship.matrix,
                mapping,
                trait,
                train.index,
                tmp_path){

  path = getwd()
  dir.create(tmp_path)
  setwd(tmp_path)

  mod = GAPIT(
    Y = Phenotype.matrix[train.index,c("GID",trait)],
    GM = mapping,
    GD = cbind(data.frame(Taxa = rownames(Marker.matrix[train.index,])),
               data.frame(Marker.matrix[train.index,])),
    KI = cbind(data.frame(Taxa = rownames(Kinship.matrix[train.index,train.index])),
               data.frame(Kinship.matrix[train.index,train.index])),
    PCA.total = 6,
    model = "Blink"
  )

  setwd(path)
  unlink(tmp_path, recursive=T)

  output = mod$GWAS %>%
    mutate(LOD = -log10(P.value)) %>%
    dplyr::select(SNP, Chromosome, Position, LOD)
  return(output)
}

