library(tidyverse)

pre_trait2go = read_csv('../../output/pre_trait2go.csv')

GO2name = na.omit(unique(pre_trait2go[,5:6]))
colnames(GO2name)[1] = "ID"

na.omit((pre_trait2go[,4:6])) %>%
  group_by(Type) %>%
  count()

tmp = na.omit(pre_trait2go[,c("Trait","Gene","ID")])
#430 ontologies -> 430/172 = 2.5 go/gene
length(unique(tmp$ID))
#217 unique gos
nrow(unique(tmp[,c("Trait","Gene")]))
#83 trait-gene with gos -> 83/172 = 48.25%
length(unique(tmp[,c("Gene")]))
#80 if we see only the gene but we need to take into account that the same gene
#can be related with two or more traits

count_go = pre_trait2go %>%
  unique() %>%
  count(Type,ID)

highest = count_go %>%
  filter(n>5) %>%
  count(Type)

p = 10 # number of ontologies per class
 
component.gos = count_go %>%
  filter(Type %in% "Gene.Ontology..cellular.component.") %>%
  arrange(-n) %>%
  .[1:p,]
function.gos = count_go %>%
  filter(Type %in% "Gene.Ontology..molecular.function.") %>%
  arrange(-n) %>%
  .[1:p,]
process.gos = count_go %>%
  filter(Type %in% "Gene.Ontology..biological.process.") %>%
  arrange(-n) %>%
  .[1:p,]


selected.gos = rbind(component.gos,function.gos,process.gos)
selected.gos$Type = rep(c("Cellular Component","Molecular Function","Biological Process"),
                        each = p)
selected.names = left_join(selected.gos,GO2name)
### step one finished!

# now, extract average by type
#component
tmp = pre_trait2go %>%
  unique() %>%
  filter(Type %in% "Gene.Ontology..cellular.component.") %>%
  count(Trait)
tmp$n = tmp$n/sum(tmp$n)
tmp$Type = "Cellular Component"
# function
tmp1 = pre_trait2go %>%
  unique() %>%
  filter(Type %in% "Gene.Ontology..molecular.function.") %>%
  count(Trait)
tmp1$n = tmp1$n/sum(tmp1$n)
tmp1$Type = "Molecular Function"
# process
tmp2 = pre_trait2go %>%
  unique() %>%
  filter(Type %in% "Gene.Ontology..biological.process.") %>%
  count(Trait)
tmp2$n = tmp2$n/sum(tmp2$n)
tmp2$Type = "Biological Process"

average.freq = rbind(tmp,tmp1,tmp2)
average.freq$Name = "average"
average.freq$freq = average.freq$n
average.freq$n = NULL

######
# now, freqs of top5 for each go category
component.freq = selected.names %>%
  filter(Type %in% "Cellular Component") %>%
  dplyr::select(ID) %>%
  left_join(pre_trait2go) %>%
  count(Trait,Name) %>%
  group_by(Name) %>%
  mutate(freq = n / sum(n))
component.freq$n = NULL
component.freq$Type = "Cellular Component"

function.freq = selected.names %>%
  filter(Type %in% "Molecular Function") %>%
  dplyr::select(ID) %>%
  left_join(pre_trait2go) %>%
  count(Trait,Name) %>%
  group_by(Name) %>%
  mutate(freq = n / sum(n))
function.freq$n = NULL
function.freq$Type = "Molecular Function"

process.freq = selected.names %>%
  filter(Type %in% "Biological Process") %>%
  dplyr::select(ID) %>%
  left_join(pre_trait2go) %>%
  count(Trait,Name) %>%
  group_by(Name) %>%
  mutate(freq = n / sum(n))
process.freq$n = NULL
process.freq$Type = "Biological Process"

df = rbind(average.freq[,c(1,3,4,2)],component.freq,function.freq,process.freq) %>%
  arrange(Trait)
# done!

component.gos = left_join(component.gos,GO2name)
function.gos = left_join(function.gos,GO2name)
process.gos = left_join(process.gos,GO2name)

name_levels = rev(c("average",component.gos$Name, function.gos$Name, process.gos$Name))

text_labels = rev(c("all",component.gos$Name, function.gos$Name, process.gos$Name))
number_labels = rev(c("NN",component.gos$n,function.gos$n,process.gos$n))


name_labels = paste0(text_labels, " (",number_labels,")")


# name_labels= rev(c("All",
#  "nucleus","membrane","cytoplasm","extracellular\nregion","cytosol",
#  "Zn2+\nbinding","ATP\nbinding", "RNA\nbinding", "sequence-specific\nDNA binding","nucleic acid\nbinding",
#  "proteolysis","cell redox\nhomeostasis", "protein\nphosphorylation","cellular response\nto oxidative stress",
#   "regulation of\nDNA-templated\ntranscription"
# ))

levs = c("S_fresh","S_dry",
  "B_height","B_width","B_shape","B_weight",
  "H_cluster","H_rachis",
  "P_cluster_weight","P_cluster_loss","P_rachis_weight")

trait_colors = c("#F7D08A","#F79F79",
                 "#D3E75F","#87D733","#68AA22","#2E4E2D",
                 "#E28DC4","#D65CAB",
                 "#99D2EB","#66BCE1","#258EBB")

labs = levs
labs[c(9,11)] = c("P_cluster","P_rachis")

df %>%
  mutate(Type = factor(Type,levels = sort(unique(Type))[c(2,3,1)])) %>%
  mutate(Trait = factor(Trait, levels = levs, labels = labs)) %>%
  mutate(Name = factor(Name,levels = name_levels,labels=name_labels)) %>%
  ggplot(aes(y=Name,x=-freq,fill=Trait,alpha=Name%in%"All"))+
  geom_bar(position="fill",stat="identity")+
  facet_wrap( ~ Type,scales="free")+
  theme_classic()+
  scale_fill_manual(values = trait_colors)+
  scale_alpha_manual(values = c(1,1))+
  theme(legend.position = "top")+
  theme(axis.ticks.x = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold",size=10),
        strip.background = element_rect(fill="black"),
        strip.text = element_text(color="white",face="bold",size=20))+
  xlab("")+ylab("")+
  guides(alpha = "none")


