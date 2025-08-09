#15-3-23 - boostrapping
# input: markers, blues, hits

M_grape = read_csv('../../data/clean/X.csv')
big.BLUEs = read_csv('../../data/clean/BLUEs.csv')
assos = read_csv('../../output/gwas_hits.csv')

output = data.frame(matrix(nrow=0,ncol=8))
colnames(output) = c("Taxa","Value","SNP","rank","i","Trait","Marker","rep")

samp = 250
rep = 200
for (r in 1:rep){
  print(r)
  for (i in 1:nrow(assos)){
    asso = assos[i,]
    tr = asso$Trait
    pheno = na.omit(big.BLUEs[,c("Geno",tr)]) %>%
      .[sample(1:nrow(.),samp,T),]
    #pheno$order = rank(pheno[,tr])
    marker = asso$Marker
    snp = M_grape[,c("Geno",marker)]
    all = na.omit(suppressMessages(left_join(pheno,snp)))
    colnames(all)[2:3] = c("Value","SNP")
    all$rank = rank(-all$Value,ties.method = "random")
    all$i = i
    all$Trait = tr
    all$Marker = marker
    all$rep = r
    output = rbind(output,all)
  }
}


### dont run more than once
tmp = output %>%
  group_by(rank,i,Trait,Marker) %>%
  summarise_at(vars(SNP), list(name = mean)) %>%
  arrange(i)

tmp %>%
  ggplot(aes(x=rank,y=i,alpha=name,fill=Trait))+
  geom_tile(color="transparent")+
  theme_classic()+
  facet_wrap(Trait~.,scales="free")


best.assos = assos %>% 
  group_by(Trait) %>% 
  mutate(best = LOD == max(LOD)) %>% 
  ungroup() %>% filter(best) %>% dplyr::select(Trait, Marker)

trait_order = c("S_fresh","S_dry",
  "B_height","B_width","B_shape","B_weight",
  "H_cluster","H_rachis",
  "P_cluster_weight","P_cluster_loss","P_rachis_weight")

trait_colors = c("#F7D08A","#F79F79",
                 "#D3E75F","#87D733","#68AA22","#2E4E2D",
                 "#E28DC4","#D65CAB",
                 "#99D2EB","#66BCE1","#258EBB")

tmp = tmp %>%
  right_join(best.assos) %>%
  mutate(Trait = factor(Trait,levels = trait_order)) #

rank.adj = c()

###3!!!!!!!!!!!!!!
### ...

change = c("S_fresh","S_dry","P_cluster_weight")
#change = c()
#
for (i in 1:nrow(tmp)){
  if (tmp[i,]$Trait %in% change){
    rank.adj = c(rank.adj, samp+1-tmp[i,]$rank)
  }
  else{
    rank.adj = c(rank.adj, tmp[i,]$rank)
  }
}
tmp$rank.adj = rank.adj
tmp$bar = "Ranking"

tmp %>%
  ggplot(aes(x=rank.adj,y=Trait,alpha=name,fill=Trait))+
  geom_tile(color="black")+
  theme_classic()+
  facet_grid(Trait~.,scales="free")+
  theme(legend.position = "top")+
  theme(axis.ticks.y = element_blank(),
        axis.text.y= element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_blank())+
  xlab("")+ylab("")+
  scale_alpha(range=c(0,1), limits=c(0,2), na.value = 1)+
  scale_fill_manual(values = trait_colors)+
  labs(alpha="Average SNP dosage\n per raking position")+
  scale_x_continuous(breaks=c(1,25,50,75,100))


# add freqs

freqs = c()
nmarks = nrow(M_grape)
for (i in 1:nrow(assos)){
  freqs = c(freqs, samp*cumsum(table(M_grape[[assos[i,]$Marker]])/nmarks))
}

df.freqs = data.frame(
  Trait = rep(assos$Trait, each = 3),
  Marker = rep(assos$Marker, each = 3),
  i = rep(1:nrow(assos), each = 3),
  Dose = rep(0:2,nrow(assos)),
  Freq = freqs)
df.freqs$Trait = factor(df.freqs$Trait, levels = trait_order)

# frequency of alleles
# df.freqs %>%
#   filter(i %in% bests) %>%
#   ggplot(aes(x=-Freq,y=as.factor(i),fill=Trait, alpha = factor(Dose)))+
#   theme_classic()+
#   geom_bar(position="stack", stat="identity",
#            color="black")+
#   scale_alpha_manual(values = c(0,0.5,1))+
#   facet_grid(Trait~.,scales="free")+
#   scale_fill_manual(values = trait_colors)+
#   theme(legend.position = "top")+
#   theme(axis.ticks.y = element_blank(),
#         axis.text.y= element_blank(),
#         axis.line.y = element_blank(),
#         strip.text = element_blank())+
#   xlab("")+ylab("")+
#   labs(alpha="Average SNP dosage\n per raking position")


tmp2 = left_join(best.assos,df.freqs)
allele_fill = c()
for (i in 1:nrow(best.assos)){
  m = best.assos[i,]$Marker
  t = best.assos[i,]$Trait
  dose.dist = tmp2 %>%
    filter(Marker == m & Trait == t) %>%
    dplyr::select(Freq) %>% cumsum()
  for (j in 1:samp){
    if (j < dose.dist$Freq[1]){
      allele_fill = c(allele_fill, 0)
    }
    else if (j < dose.dist$Freq[2]){
      allele_fill = c(allele_fill, 1)
    }
    else{
      allele_fill = c(allele_fill, 2)
    }
  }
}

tmp$allele = NULL

tmp2 = tmp
tmp2$name = allele_fill
tmp2$rank.adj = tmp2$rank
tmp2$bar = "Distribution"

full = rbind(tmp,tmp2)

full %>%
  mutate(bar = factor(bar,labels = c("SNP\nDistribution","Ranking"))) %>%
  ggplot(aes(x=rank.adj,y=bar,fill=Trait,alpha=name))+
  geom_tile(aes(color=bar))+
  facet_grid(Trait~.,scales="free")+
  scale_color_manual(values = c("transparent","black"))+
  theme_classic()+
  facet_grid(Trait~.,scales="free")+
  scale_fill_manual(values = trait_colors)+
  theme(legend.position = "top")+
  theme(#axis.ticks.y = element_blank(),
        #axis.text.y= element_blank(),
        #axis.line.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        strip.text = element_blank())+
  xlab("")+ylab("")+
  labs(alpha="Average SNP dosage\n per raking position")+
  guides(color="none")+
  scale_x_continuous(breaks=c(1,25,50,75,100))+
  scale_alpha(range=c(0.1,1), limits=c(0,2), na.value = 1)+
  #theme(panel.background = element_rect(color="black"))+
  geom_hline(yintercept = 1.5, linetype = "dashed")


avg.df = data.frame(
  Trait = trait_order,
  Value = colMeans(M_grape[,best.assos$Marker[match(trait_order,best.assos$Trait)]])
)


trait_label = c("S_fresh", "S_dry","B_height","B_width","B_shape","B_weight",
                "H_cluster","H_rachis","P_cluster","P_cluster_loss","P_rachis")
full$Trait = factor(full$Trait, levels = trait_order, labels = trait_label)
avg.df$Trait = factor(avg.df$Trait, levels = trait_order, labels = trait_label)

full$TM = paste(full$Trait,full$Marker,sep=" | ")
avg.df$TM = paste(avg.df$Trait, rownames(avg.df), sep=" | ")

full$TM = factor(full$TM, levels = sort(unique(full$TM))[c(11,10,1,4,2,3,5:9)])
avg.df$TM = factor(avg.df$TM, levels = sort(unique(avg.df$TM))[c(11,10,1,4,2,3,5:9)])

full %>%
  filter(bar %in% "Ranking") %>%
  ggplot(aes(x=rank,y=name,color=TM,fill=TM,group=i))+
  geom_area(alpha=0.7)+
  geom_line(alpha=0.8,size=1.1,show.legend = F)+
  facet_wrap(TM~.,nrow=4)+
  theme_classic()+
  theme(strip.background = element_rect(fill="black"),
        strip.text = element_text(size=12,color="white",face="bold"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))+
  scale_fill_manual(values = trait_colors)+
  scale_color_manual(values = trait_colors)+
  geom_hline(data = avg.df, aes(yintercept = Value),
             linetype = "dashed",show.legend = F)+
  scale_y_continuous(breaks = c(0,1,2))+
  scale_x_continuous(breaks = c(1,50,100,150,200,250))+
  ylab("Average genotypic value")+xlab("Phenotypic ranking")+
  theme(legend.position = "none")+
  theme(legend.key.size = unit(11,"mm"))

for (tr in trait_order){
  a = filter(full, Trait == tr)
  print(cor(a$rank, a$name, method = "spearman"))
}

full %>%
  filter(bar == "Ranking") %>%
  group_by(i, Trait, Marker) %>%
  summarise(spe = cor(rank,name,method="spearman"))

RR = 200
SS = 250
MM = 10
output = c()
# same but random markers
for (tr in trait_order){
  for (m in 1:MM){
    marker = colnames(M_grape)[sample(2:ncol(M_grape),1)]
    print(c(tr))
    for (rr in 1:RR){
      pheno = na.omit(big.BLUEs[,c("Geno",tr)]) %>%
        .[sample(1:nrow(.),SS,T),]
      snp = M_grape[,c("Geno",marker)]
      all = na.omit(suppressMessages(left_join(pheno,snp)))
      colnames(all)[2:3] = c("Value","SNP")
      all$rank = rank(-all$Value,ties.method = "random")
      all$Trait = tr
      all$Marker = marker
      all$rep = rr
      output = rbind(output,all)
    }
  }
}

a = output %>%
  filter(rep<10) %>%
  group_by(Trait,Marker,rank) %>%
  summarise(name = mean(SNP)) %>%
  group_by(Trait, Marker) %>%
  summarise(cor = cor(rank,name,method="spearman"))

a %>%
  mutate(Trait = factor(Trait,levels = trait_order, labels = trait_label)) %>%
  ggplot(aes(x=Trait, y=cor))+
  geom_hline(yintercept=0, color = "black", linetype="dashed", size=1.1)+
  geom_boxplot(aes(fill=Trait,color=Trait),alpha=0.6)+
  xlab("")+ylab("Average Spearman's correlation")+
  theme_classic()+
  scale_fill_manual(values = trait_colors)+
  scale_color_manual(values = trait_colors)+
  theme(legend.position = "none")

#+  ggtitle("Average Spearman's correlation (R = 100, N = 20)")

b = a %>%
  group_by(Trait) %>%
  summarise(meancor = mean(abs(cor)))


hist(a$cor)
# check
# a = full %>%
#   group_by(Trait,Marker,bar) %>%
#   summarise(mean = mean(name)) %>%
#   arrange(Trait,Marker)
#
# sig.markers = unique(full$Marker)
# colMeans(M_grape[,sig.markers])
#
# tmp.lines = M_grape %>%
#   filter(M_11_13989722 %in% 0) %>%
#   select(Taxa) %>% unlist()
#
# big.BLUEs %>% filter(Taxa %in% tmp.lines) %>%
#   select(S_fresh) %>% unlist () %>% mean(na.rm = T)

### check 929
pheno2 = big.BLUEs[,c("Geno","S_fresh")]
marker2 = M_grape[,c("Geno","M_11_13989722")]

famili2 = big.BLUEs[,c('Geno', 'family')]

tmp3 = pheno2 %>% right_join(marker2) %>% right_join(famili2) %>%
  arrange(-S_fresh) %>%
  mutate(rank = rank(-S_fresh)) #%>% .[1:50,]

tmp7 = tmp3[1:50,c(1,4,2,5,3)]



table(tmp3$family)
table(tmp3$M_11_13989722)

tmp4 = pheno2 %>% right_join(marker2) %>% right_join(famili2) %>%
  arrange(-S_fresh) %>% filter(family == "929")
table(tmp4$M_11_13989722)

tmp5 = pheno2 %>% right_join(marker2) %>% right_join(famili2) %>%
  arrange(-S_fresh) %>% filter(family == "jardin")
table(tmp5$M_11_13989722)
