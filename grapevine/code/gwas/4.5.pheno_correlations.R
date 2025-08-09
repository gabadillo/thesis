library(tidyverse)

big.BLUEs = read_csv('../../data/clean/BLUEs.csv')

trait.order = c("S_fresh","S_dry",
                "B_height","B_width","B_shape","B_weight",
                "H_cluster","H_rachis",
                "P_cluster_weight","P_cluster_loss","P_rachis_weight")

tmp = cor(big.BLUEs[,3:15],use = "pairwise.complete.obs")

#### flags 3#########################
white = 0.3
gray = 0.7

type_list = list(c("jardin"),c("406","411","929","900",
                               "902","912","111"))

titles = c("Cultivars","Breeding lines")
plot_list = list()

## filter
# type
for (cnt in 1:length(type_list)){
  types = type_list[[cnt]]
  lines = big.BLUEs$Geno[big.BLUEs$family %in% types]

  blues = filter(big.BLUEs, Geno %in% lines) %>%
    .[,-c(1:2,3,15)]
  blues = blues[,sort(colnames(blues))]

  full = list(cors = cor(blues, use = "pairwise.complete.obs"),
  PCA = prcomp(na.omit(blues))$rotation[,1:3])

  ####
  ####

  ##  plot network ###################
  # nodes
  df.nodes = data.frame(full$PCA)
  df.nodes$Trait = trait.order
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
  trait.order = c("S_fresh","S_dry",
                  "B_height","B_width","B_shape","B_weight",
                  "H_cluster","H_rachis",
                  "P_cluster_weight","P_cluster_loss","P_rachis_weight")

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

ggpubr::ggarrange(plotlist = plot_list, ncol = 1)

