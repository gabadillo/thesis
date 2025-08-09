
library(tidyverse)

mapping = read_csv('../../data/clean/mapping.csv')
assos = read_csv('../../output/gwas_hits.csv')
raw_genes = read_csv('../../output/raw_genes.csv')

raw_genes = filter(raw_genes, !dist2snp > 10000)


trait.order = c("S_fresh","S_dry",
                "B_height","B_width","B_shape","B_weight",
                "H_cluster","H_rachis",
                "P_cluster_weight","P_cluster_loss","P_rachis_weight")

assos.chrom = raw_genes %>%
  dplyr::select(1:8) %>%
  unique() %>%
  count(Chr,Trait)


chromo_maxpos =
  mapping %>%
  dplyr::select(CHROM,POS) %>%
  group_by(CHROM) %>%
  summarise(max=max(POS))
chromo_maxpos$Chr = as.numeric(chromo_maxpos$CHROM)
chromo_maxpos = na.omit(chromo_maxpos)

raw_genes$delta = abs(rnorm(nrow(raw_genes),0,0.1))+0.15
raw_genes$Trait = factor(raw_genes$Trait,levels=trait.order)

mbp = 1e6
trait_colors = c("#F7D08A","#F79F79",
                 "#D3E75F","#87D733","#68AA22","#2E4E2D",
                 "#E28DC4","#D65CAB",
                 "#99D2EB","#66BCE1","#258EBB")

ggplot()+
  geom_segment(data=chromo_maxpos,
               aes(x=0, xend=max/mbp,y=Chr,yend=Chr),
               size=0.2)+
  geom_point(data=raw_genes,
             aes(x=(Start+End)/2/mbp,y=Chr-0.2,fill=Trait),
             shape=24,color="black",size=3,alpha=1)+
  geom_point(data=raw_genes,
             aes(x=(gene_start+gene_end)/2/mbp,y=Chr+delta,color=Trait),
             size=0.1,alpha=0.8,shape=8)+
  scale_fill_manual(values=trait_colors)+
  scale_color_manual(values=trait_colors)+
  # geom_segment(data=raw_genes,
  #              aes(x=gene_start/mbp,xend=gene_end/mbp,y=Chr+delta,yend=Chr+delta,
  #                  color = Trait))+
  theme_classic()+
  scale_x_continuous(breaks = seq(0,30,5))+
  scale_y_continuous(breaks = seq(19,1,-1))+
  xlab("Mbp")

### all together
mapping.chr = filter(mapping, CHROM != "Un") %>%
  arrange(CHROM,POS)

curr.chrom = "01"
curr.pos = 0

new.pos = c(mapping.chr[1,]$POS)
for (i in 2:nrow(mapping.chr)){
  if (!curr.chrom %in% mapping.chr[i+1,]$CHROM){
    curr.pos = curr.pos + mapping.chr[i,]$POS
    new.pos = c(new.pos, curr.pos)
    curr.chrom = mapping.chr[i,]$CHROM
  }
  else{
    new.pos = c(new.pos,curr.pos + mapping.chr[i,]$POS)
  }
}

mapping.chr$new.pos = new.pos
delta = 0
chromo_maxpos$init = c(-delta,cumsum(chromo_maxpos$max))[1:19]+delta
chromo_maxpos$fin = c(cumsum(chromo_maxpos$max))

tmp = assos[,c("Trait","Marker","LOD")]
colnames(tmp)[2] = "ID"
df = unique(left_join(tmp,mapping.chr))
df$delta = rnorm(nrow(df),0.1,0.5)
trait_colors = c("#F7D08A","#F79F79",
                 "#D3E75F","#87D733","#68AA22","#2E4E2D",
                 "#E28DC4","#D65CAB",
                 "#99D2EB","#66BCE1","#258EBB")
mbp = 0.5*1e6

window = data.frame(new.pos = mbp*mapping.chr$new.pos%/%(mbp))
window = window %>% count(new.pos)

#### check hits with no genes !!! 11-4-23
hits_with_genes = unique(raw_genes[,c("Trait","Marker")])
addregions = assos[which(!assos$Marker %in% hits_with_genes$Marker),1:2]

#####
delta = 5e5
regions.counts = raw_genes[,c(4:5)] %>%
  rbind(addregions) %>%
  left_join(mapping.chr,by=c("Marker"="ID")) %>%
  mutate(new.init = new.pos-delta,
         new.fin  = new.pos+delta)

regions.genes = raw_genes[,c(4:5)] %>%
  #rbind(addregions) %>%
  left_join(mapping.chr,by=c("Marker"="ID")) %>%
  mutate(new.init = new.pos-delta,
         new.fin  = new.pos+delta)

regions = regions.counts %>%
  count(Trait,Marker,new.init,new.fin)

###### also new approach aaaaaa
# regionswithoutgenes = c()
# for (i in 1:nrow(assos)){
#   indextraits = which(regions$Trait %in% assos$Trait[i])
#   indexmarkers = which(regions$Marker %in% assos$Marker[i])
#   indexintersect = intersect(indextraits,indexmarkers)
#   if (length(indexintersect)==0){
#     regionswithoutgenes = c(regionswithoutgenes,i)
#   }
# }
# addregions = assos[regionswithoutgenes,1:2]
# save(addregions, file = "postGWAS/hits_no_genes.RData")
####### dont touch hereee


eps = 1e7
tmp = c()
curr.pos = regions[1,]$new.init
cnt = 0
for (i in 1:nrow(regions)){
  if ((regions[i,]$new.init-curr.pos)>eps){
    curr.pos = regions[i,]$new.init
    tmp = c(tmp,0)
    cnt = 0
  }
  else{
    cnt = cnt + 1
    tmp = c(tmp,cnt%%3)
  }
}
regions$jitter = tmp

#regions = regions[1:nrow(raw_genes),]
trait.order = c("S_fresh","S_dry",
                "B_height","B_width","B_shape","B_weight",
                "H_cluster","H_rachis",
                "P_cluster_weight","P_cluster_loss","P_rachis_weight")

# tmp comes from snps with no genes 6.7

df$Trait = factor(df$Trait,levels = trait.order)
levels(regions$Trait) = trait.order
levels(regions.genes$Trait) = trait.order


ggplot()+
  geom_rect(data=chromo_maxpos,
            aes(xmin=init,xmax=fin,alpha=as.factor(as.numeric(CHROM)%%4)),
            ymin=log(6),ymax=3)+
  geom_segment(data=df,
    aes(y=1, yend=log(LOD),color=Trait,x=new.pos,xend=new.pos),
    size=1,alpha=0.7)+
  geom_point(data=df,
             aes(x=new.pos,y=log(LOD),color=Trait,size=log(LOD)),alpha=0.7)+
  theme_classic()+
  scale_alpha_manual(values = c(0.05,0.15,0.25,0.35))+
  coord_polar(theta="x")+
  scale_color_manual(values = trait_colors)+
  scale_fill_manual(values = trait_colors)+
  guides(alpha="none")+
  theme(legend.position = "none")+
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_blank())+
  xlab("")+ylab("")+
  geom_text(data=chromo_maxpos,
            aes(label=paste0("Chr",CHROM),x=(init+fin)/2,y=2.75))+
  geom_line(data=window,
            aes(x=new.pos,y=3.15+(n/500)),alpha=0.6,color="#65334D")+
  geom_segment(data = regions,
               aes(x=new.init,xend=new.fin,color=Trait),
               y=3.8,yend=3.8,size=6.5,alpha=0.8)+
  geom_hline(yintercept=3.7)+
  geom_hline(yintercept=3.9)+
  geom_jitter(data = regions.genes,
               aes(x=new.pos,y=4.5,fill=Trait),shape=21)




