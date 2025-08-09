library(httr)
library(jsonlite)
library(xml2)
library(tidyverse)

assos = read_csv('../../output/gwas_hits.csv')
raw_genes = read_csv('../../output/raw_genes.csv') %>% filter(!dist2snp > 10000)

genes = c()
entrys = c()
types = c()
values = c()

# i remove these because unknown erros, just by trial and error
#remove_by_exp = c("Vitvi02g00330","Vitvi04g01447","Vitvi08g01519","Vitvi05g00043")
remove_by_exp = c("Vitvi02g00330")

for (url_gene in raw_genes$gene){
  if(!url_gene %in% remove_by_exp){
    url_init = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cft_intramem%2Ccc_subcellular_location%2Cft_topo_dom%2Cft_transmem%2Ccc_allergen%2Ccc_biotechnology%2Ccc_disruption_phenotype%2Ccc_disease%2Cft_mutagen%2Ccc_pharmaceutical%2Ccc_toxic_dose%2Cgo_p%2Cgo_c%2Cgo_f%2Ccc_tissue_specificity%2Ccc_developmental_stage%2Ccc_induction%2Ccc_interaction%2Cft_binding%2Cft_dna_bind%2Ccc_cofactor%2Cft_site%2Cft_act_site&format=tsv&query=%28"
    url_fin = "%29"
    url = paste0(url_init,url_gene,url_fin)

    r = GET(url)
    raw = read.table(text=content(r),header=TRUE, sep="\t")

    if (nrow(raw)==0){
      genes = c(genes,url_gene)
      entrys = c(entrys,NA)
      types = c(types, NA)
      values = c(values,NA)
    }
    else{
      for (i in 1:nrow(raw)){
        entry = raw[i,]$Entry
        for (j in 2:ncol(raw)){
          if (!is.na(raw[i,j])){
            genes = c(genes,url_gene)
            entrys = c(entrys,entry)
            types = c(types, colnames(raw)[j])
            values = c(values,raw[i,j])
          }
        }
      }
    }
  }
}



gene2info = data.frame(
  Gene = genes,
  Entry = entrys,
  Type = types,
  Value = values
)

entrys = c()
types = c()
IDs = c()
names = c()

# process info
unique(gene2info$Type)
# 10 categories
# 1. subcellular -> remove SUBCELLULAR LOCATION:
for (i in 1:nrow(gene2info)){
  if (gene2info[i,]$Type %in% "Subcellular.location..CC."){
    #Remove
    clean = strsplit(gene2info[i,]$Value, "SUBCELLULAR LOCATION: ")[[1]][2]
    split_point = strsplit(clean, "\\.")[[1]]
    for (segp in split_point){
      split_semicolon = strsplit(segp,";")[[1]]
      for (segs in split_semicolon){
        split_bracket = strsplit(segs, "\\{")[[1]]
        entrys = c(entrys, gene2info[i,]$Entry)
        types = c(types, gene2info[i,]$Type)
        IDs = c(IDs,gsub("\\}","",split_bracket[2]))
        names = c(names,split_bracket[1])
      }
    }
  }
  #2. DNA.binding
  else if (gene2info[i,]$Type %in% "DNA.binding"){
    split_evidence = strsplit(gene2info[i,]$Value,"; /evidence=")[[1]]
    entrys = c(entrys, gene2info[i,]$Entry)
    types = c(types, gene2info[i,]$Type)
    IDs = c(IDs,split_evidence[2])
    names = c(names,strsplit(split_evidence[1],"/note=")[[1]][2])
  }
  #3. Cofactor  -> remove COFACTOR:
  else if (gene2info[i,]$Type %in% "Cofactor"){
    clean = strsplit(gene2info[i,]$Value, "COFACTOR: ")[[1]][2]
    split_semicolon = strsplit(clean,";")[[1]]
    for (segs in split_semicolon){
      if (grepl("Name=",segs)){
        clean_name = strsplit(segs,"=")[[1]][2]
      }
      else if (grepl("Evidence=",segs)){
        clean_id = gsub("\\}","",strsplit(segs,"=\\{")[[1]])[2]
      }
    }
    entrys = c(entrys, gene2info[i,]$Entry)
    types = c(types, gene2info[i,]$Type)
    IDs = c(IDs,clean_id)
    names = c(names,clean_name)
  }
  #4. Active.site
  else if (gene2info[i,]$Type %in% "Active.site"){
    split_semicolon = strsplit(gene2info[i,]$Value,";")[[1]]
    clean_id = NA
    clean_name = NA
    for (segs in split_semicolon){
      if (grepl("evidence",segs)){
        clean_id = strsplit(segs,"=")[[1]][2]
      }
      else if (grepl("note",segs)){
        clean_name = strsplit(segs,"=")[[1]][2]
      }
    }
    entrys = c(entrys, gene2info[i,]$Entry)
    types = c(types, gene2info[i,]$Type)
    IDs = c(IDs,clean_id)
    names = c(names,clean_name)
  }
  #5. Binding.site -> EACH binding is a different info
  else if (gene2info[i,]$Type %in% "Binding.site"){
    split_BINDING = strsplit(gene2info[i,]$Value,"BINDING")[[1]]
    for (segB in split_BINDING){
      split_semicolon = strsplit(segB,";")[[1]]
      clean_id = NA
      clean_name = NA
      for (segs in split_semicolon){
        if (grepl("ligand=",segs)){
          clean_name = strsplit(segs,"=")[[1]][2]
        }
        else if (grepl("evidence=",segs)){
          clean_id = strsplit(segs,"=")[[1]][2]
        }
      }
      entrys = c(entrys, gene2info[i,]$Entry)
      types = c(types, gene2info[i,]$Type)
      IDs = c(IDs,clean_id)
      names = c(names,clean_name)
    }
  }
  #6. Site
  else if (gene2info[i,]$Type %in% "Site"){
    split_complex = strsplit(gene2info[i,]$Value,"; /")[[1]]
    for (segc in split_complex){
      if (grepl("note=",segs)){
        clean_name = strsplit(segs,"=")[[1]][2]
      }
      else if (grepl("evidence=",segs)){
        clean_id = strsplit(segs,"=")[[1]][2]
      }
    }
    entrys = c(entrys, gene2info[i,]$Entry)
    types = c(types, gene2info[i,]$Type)
    IDs = c(IDs,clean_id)
    names = c(names,clean_name)
  }
  #7. Transmembrane -> transform to binary (YES/NO)
  else if (gene2info[i,]$Type %in% "Transmembrane"){
    entrys = c(entrys, gene2info[i,]$Entry)
    types = c(types, gene2info[i,]$Type)
    IDs = c(IDs,"TM")
    names = c(names,"Transmembrane")
  }
  #8:10. Gene Ontologies
  else if (gene2info[i,]$Type %in% c("Gene.Ontology..cellular.component.",
                                     "Gene.Ontology..molecular.function.",
                                     "Gene.Ontology..biological.process.")){
    if (nchar(gene2info[i,]$Value) == 0){
      entrys = c(entrys, gene2info[i,]$Entry)
      types = c(types, gene2info[i,]$Type)
      IDs = c(IDs,NA)
      names = c(names,NA)
    }
    else{
      split_semicolon = strsplit(gene2info[i,]$Value,"; ")[[1]]
      for (segs in split_semicolon){
        split_GO = strsplit(segs," \\[")[[1]]
        entrys = c(entrys, gene2info[i,]$Entry)
        types = c(types, gene2info[i,]$Type)
        IDs = c(IDs,gsub("\\]","",split_GO[2]))
        names = c(names,split_GO[1])
      }
    }
  }
}

join.df = data.frame(
  Entry = entrys,
  Type = types,
  ID = IDs,
  Name = names
)

join.df = na.omit(unique(join.df))
mega.df = unique(right_join(gene2info, join.df))
mega.df$Value = NULL

GO.df = mega.df[grep("Gene.Ontology",mega.df$Type),]
#GO.df$Entry = NULL
GO.db.df = unique(GO.df[,c(3:5)])

# structure of GOs
GO.db.df %>%
  group_by(Type) %>%
  count() %>%
  mutate(rel = n/nrow(GO.df))

#
#### here we start
trait2gene = raw_genes %>%
  select(Trait,gene)
colnames(trait2gene)[2] = "Gene"

pre_trait2go = full_join(trait2gene,GO.df) #!!

f = na.omit(pre_trait2go[,c("Trait","Gene","ID")])
nrow(unique(f[,c("Trait","Gene")]))
length(unique(f$ID))
length(unique(pre_trait2go$ID))

write_csv(pre_trait2go, file = '../../output/pre_trait2go.csv')
