library(tidyverse)

pre_trait2go = read_csv('../../output/pre_trait2go.csv')
GO.table = unique(pre_trait2go[,3:5])
pre_trait2go = pre_trait2go %>% na.omit() %>% unique()

# full table
n = length(unique(pre_trait2go$Trait))
p = length(unique(pre_trait2go$ID))

trait2go = matrix(0,nrow=n,ncol=p)
rownames(trait2go) = c("S_fresh","S_dry",
                                     "B_height","B_width","B_shape","B_weight",
                                     "H_cluster","H_rachis",
                                     "P_cluster_weight","P_cluster_loss","P_rachis_weight")
colnames(trait2go) = sort(unique(pre_trait2go$ID))

for (i in 1:nrow(pre_trait2go)){
  tr = pre_trait2go[i,]$Trait
  go = pre_trait2go[i,]$ID
  trait2go[tr,go] = trait2go[tr,go] + 1
}

colnames(trait2go) = sort(unique(pre_trait2go$ID))

#### flags #########################
white = 0.3
gray = 0.7

type_list = list(c("Gene.Ontology..cellular.component."),
               c("Gene.Ontology..molecular.function."),
               c("Gene.Ontology..biological.process."),
               c( "Gene.Ontology..cellular.component.",
                  "Gene.Ontology..molecular.function.",
                  "Gene.Ontology..biological.process.")
               )
titles = c("Cellular component","Molecular function","Biological process","All Gene Ontologies")
plot_list = list()
## filter
# type
trait.order = c("S_fresh","S_dry",
                "B_height","B_width","B_shape","B_weight",
                "H_cluster","H_rachis",
                "P_cluster_weight","P_cluster_loss","P_rachis_weight")

for (cnt in 1:length(type_list)){
  types = type_list[[cnt]]
  go.types = GO.table$ID[GO.table$Type %in% types]

  column.counts = colSums(trait2go)
  min_count  = 2
  go.counts = names(which(column.counts>min_count))

  go.filtered = unique(c(intersect(go.types,go.counts)))
  #print(length(go.filtered))
  filt.trait2go = trait2go[,go.filtered]

  #ones = matrix(1,ncol=1,nrow=nrow(filt.trait2go))
  #mnorm = kronecker(ones,matrix(column.counts[go.filtered],nrow=1))
  #norm.filt.trait2go = filt.trait2go/mnorm

  #filt.trait2go[filt.trait2go!=0] = 1
  a = cor(t(filt.trait2go))
  a[is.na(a)] = 0
  # b = cbind(a,rep(0,nrow(a)))
  # c = rbind(b,c(rep(0,ncol(a)),1))
  # rownames(c)[11] = colnames(c)[11] = "B_weight"
  # d = c[trait.order, trait.order]

  full = list(cors = a, PCA = prcomp(a)$rotation[,1:3])
  ####
  ####

  ##  plot network ###################
  # nodes
  df.nodes = data.frame(full$PCA)
  df.nodes$Trait = rownames(df.nodes)
  rownames(df.nodes) = 1:nrow(df.nodes)

  # edges
  m = full$cors
  df.edges = data.frame(row=rownames(m)[row(m)], col=colnames(m)[col(m)], corr=c(m))

  #using arbitrary circle points
  n = nrow(df.nodes)
  r = 3
  angles = ((1:n)-1)/((n)/360)
  circ.x = c()
  circ.y = c()
  for (i in 1:n){
    circ.x = c(circ.x, r*cos(angles[i]*pi/180))
    circ.y = c(circ.y, r*sin(angles[i]*pi/180))
  }
  df.nodes$X = circ.x
  df.nodes$Y = circ.y

  m = full$cors
  df.edges = data.frame(row=rownames(m)[row(m)], col=colnames(m)[col(m)], corr=c(m))

  # edge attributes
  x1 = c()
  x2 = c()
  y1 = c()
  y2 = c()
  show = c()

  for (i in 1:nrow(df.edges)){
    n1 = df.edges[i,]$row
    n2 = df.edges[i,]$col

    i1 = which(df.nodes$Trait==n1)
    i2 = which(df.nodes$Trait==n2)

    x1 = c(x1, df.nodes[i1,]$X)
    y1 = c(y1, df.nodes[i1,]$Y)
    x2 = c(x2, df.nodes[i2,]$X)
    y2 = c(y2, df.nodes[i2,]$Y)

    show = c(show, i1>i2 )
  }

  df.edges$x1 = x1
  df.edges$x2 = x2
  df.edges$y1 = y1
  df.edges$y2 = y2
  df.edges$show = show

  df.edges$highlight = df.edges$corr > gray
  df.nodes$Trait = factor(df.nodes$Trait, levels = trait.order)
  trait_colors = c("#F7D08A","#F79F79",
                   "#D3E75F","#87D733","#68AA22","#2E4E2D",
                   "#E28DC4","#D65CAB",
                   "#99D2EB","#66BCE1","#258EBB")


  final_plot = ggplot()+
    geom_segment(data=filter(df.edges,show,corr>white),
                 aes(x=x1,xend=x2,y=y1,yend=y2,
                     size=abs(corr),alpha=highlight))+
    # geom_label(data=df.nodes,
    #                  aes(x=X,y=Y,label=Trait))+
    geom_point(data=df.nodes,aes(x=X,y=Y,fill=Trait),
               size=12,shape=21)+
    theme_classic()+
    scale_size("corr", range = c(0,5),limits = c(0,1))+
    scale_alpha_manual(values = c(0.1,0.7))+
    scale_fill_manual(values = trait_colors)+
    ggtitle(titles[cnt])+
    theme(legend.position = "none")+
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank())+xlab("")+ylab("")
  plot_list[[cnt]] = final_plot
}

ggpubr::ggarrange(plotlist = plot_list)


# 12-04-2023



