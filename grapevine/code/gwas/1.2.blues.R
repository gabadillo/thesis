##### generate dataframes for each trait with appropiate covariables
#Rasgos de semilla: Sin covariable
#Rasgos de baya:  seed dry weight + SS si sirve (para ratio)
#Rasgos de racimo cosecha: SS + seed dry weight
#Rasgos de racimo a postcosecha: SS

library(tidyverse)

files = dir('../../data/clean/phenotype', full.names = TRUE)
phenoT = lapply(files, read_csv) %>% lapply(as.data.frame)
names(phenoT) = str_match(files, '^.+/(.+).csv$')[,2]

formulas = list(
  "S" = "season + gid",
  "B" = "season + gid + S_dry + soluble_solids",
  "H" = "season + gid + S_dry + soluble_solids",
  "P" = "season + gid + soluble_solids"
)

# get all gids
unique.gids = c()
for (trait in phenoT){
  unique.gids = unique(c(unique.gids,trait$gid))
}

# generate big data.frame
# we need to add gid, removing later for joins
big.BLUEs = data.frame(gid = paste("gid",sort(unique.gids),sep=""))

for (trait in names(phenoT)){
  data = phenoT[[trait]]
  group = strsplit(trait,"_")[[1]][1]
  form = as.formula(paste(trait,"~",formulas[group]))
  mod = lm(form,data)
  hist(mod$residuals,main = trait, breaks = 50)
  mean = unname(mod$coefficients[1])
  blues = mod$coefficients[grep("gidG_",names(mod$coefficients))]
  blues = c(blues,"gidG_23" = 0)
  blues = blues + mean
  df.blues = data.frame(
    gid = names(blues),
    trait = blues)
  colnames(df.blues)[2] = trait
  big.BLUEs = left_join(big.BLUEs,df.blues)
}

big.BLUEs$gid = str_sub(big.BLUEs$gid,4)

#write_csv(big.BLUEs, file = '../../data/clean/preBLUEs.csv') # family has to be added after processing genomics

# load final BLUEs
big.BLUEs = read_csv('../../data/clean/BLUEs.csv')

#
levs = c("B_height","B_shape","B_weight","B_width",
         "H_cluster","H_rachis",
         "P_cluster_loss", "P_cluster_weight","P_rachis_loss", "P_rachis_weight",
         "S_dry","S_fresh","S_number")

labs = levs

# families
color_names <- c("111", "406", "411", "900", "902", "912", "929", "jardin", "self")
color_values <- c("#D64E12", "#F9A52C", "#EFDF48", "#8BD346", "#60DBE8", "#16A4D8", "#9B5FE0", "pink", "gray")
named_colors <- setNames(color_values, color_names)

family_levels = c('self', '406', '411', '900', '912', '111', '902', '929', 'jardin')

# figure plot
B = big.BLUEs %>%
  na.omit() %>%
  gather(trait,value,3:15) %>%
  mutate(trait = factor(trait,levels = levs,labels=labs)) %>%
  mutate(family = factor(family, levels = family_levels)) %>% 
  ggplot(aes(x=trait,y=value,fill=family))+
  geom_boxplot()+
  theme_classic()+
  theme(legend.position = "top")+
  facet_wrap(~trait,scales="free",ncol=3)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+xlab("")+ylab("")+
  scale_fill_manual(values = named_colors,name="")+
  theme(
    strip.background = element_rect(
      fill="black", color="black", size=1, linetype="solid"
    ))+
  theme(
    strip.text.x = element_text(
      size = 10, color = "white", face="bold"
    ),
    axis.text.y = element_text(
      size = 10, face="bold"
    )
  )

B

