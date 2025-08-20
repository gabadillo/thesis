library(tidyverse)
#sparse testing representation

# INPUT
n.g = 12  # no of genotyes per "group"
n.e = 3 # no of environments

# total budget is n.g*n.e plots
# total combinations is n.g*(n.e)^2
checks = c(3,7,11,15,19,23,27,31) # how many checks we want to have - omit 0

checks = 1:12

# this should be a huge function
design = data.frame(
  env = rep(1:n.e,each=n.g),
  gen = 1:(n.g*n.e),
  ind = rep(1:n.g,n.e),
  status = factor("line",levels=c("test","check","line")))
  
m = kronecker(diag(n.e),matrix(1,ncol=1,nrow=n.g))
m_list = list()
m_list[[paste(n.g,0,sep="/")]] = list(m=m,design=design,uniq=n.e*n.g,avg.reps=1)

for (r in 1:n.g){
  #remove n.e-1 lines 
  for (r2 in 1:(n.e-1)){
    genomax = (max(design$ind[which(apply(m,1,sum)==1)]))
    remove = max(intersect(which(design$status=="line"),which(design$ind==genomax)))
    m[remove,] = 0
    design$status[remove] = "test"
  }
  
  #add check
  genomin = (min(design$ind[which(apply(m,1,sum)==1)]))
  add = min(intersect(which(design$status=="line"),which(design$ind==genomin)))
  m[add,] = 1
  design$status[add] = "check"
  
  if (r %in% checks){
    print(table(design$status))
    print(paste0("total plots: ",sum(m)))
    m_list[[paste(n.g-r,r,sep="/")]] = list(design=design,m=m,
                                            uniq = sum(table(design$status)*c(0,1,1)),
                                            avg.reps = sum(table(design$status)*c(0,n.e^2,1)/c(n.g*n.e)))
  }
}

plot_design = function(item,title="",subtitle=""){
  m = item$m
  design = item$design
  n = nrow(m); p = ncol(m)
  aux = data.frame(xintercept=0:p+0.5,
                   yintercept=0:-p/(p/n))
  
  colnames(m) = paste0("E",1:p)
  
  mdesign = cbind(design,m)
  mmelt = mdesign %>% gather(env,value,(1:n.e)+4)
  mmelt$status = factor(mmelt$status,levels = c("test","check","line"), 
                        labels = c('test', 'OG', 'NOG'))
  mmelt$env = factor(mmelt$env, levels = colnames(m))
  
  mmelt %>%
    ggplot()+
    geom_tile(aes(y=-gen,x=env,alpha=as.character(value),
                  fill=status), color = 'black')+
    geom_vline(data=aux,aes(xintercept=xintercept))+
    geom_hline(data=aux,aes(yintercept=yintercept))+
    theme_classic()+
    theme(
      legend.position="none",
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(face="bold"))+
    xlab("")+ylab("")+
    scale_alpha_manual(values = c(0,1))+
    scale_fill_manual(values = c(test="white",
                                 OG="#FA4616",
                                 NOG="#0021A5"),name="")+
    labs(title = title, subtitle = subtitle)
}

plot_list = list()
for (st in names(m_list)){
  title = st
  subtitle = paste0("unique genotypes: ",m_list[[st]]$uniq,
                 "\naverage phenotypes per genotype: ",round(m_list[[st]]$avg.reps,2))
  plot_list[[st]] = plot_design(m_list[[st]],title,subtitle)
}

ggpubr::ggarrange(plotlist = plot_list)

table(m_list$`20/11`$design$status)

#get the legend
# add color="black" in line 63
plot_list[[2]]+theme(legend.position = "top")+guides(alpha="none")

# gif
library(magick)
# --- 2. Define Output Directory and File Prefix ---
output_dir <- "temp_plots" # A temporary directory for your individual image files
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
file_prefix <- "plot_"

# --- 3. Save Each ggplot Object as a Separate Image ---
image_files <- character(length(plot_list))

for (i in seq_along(plot_list)) {
  file_name <- file.path(output_dir, paste0(file_prefix, sprintf("%03d", i), ".png"))
  ggsave(filename = file_name, plot = plot_list[[i]],
         width = 6, height = 4, units = "in", dpi = 150) # Adjust dimensions and DPI as needed
  image_files[i] <- file_name
}

# --- 4. Create the GIF using magick ---
# Read the images
images <- image_read(image_files)
#'delay' is the time in 1/100ths of a second between frames (e.g., 100 = 1 second)
# 'loop' sets the number of times the GIF should loop (0 for infinite loop)
gif_name <- "../output/sparseTesting.gif"
image_animate(images, fps = 1, loop = 0) %>% # fps is frames per second
  image_write(path = gif_name)

unlink('temp_plots/', recursive = TRUE)
