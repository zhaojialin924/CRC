###cor##
###plot network####
setwd("D:/Postgraduation/microbe_crc")

library(readxl)
library(ggthemes)
library(RColorBrewer)
library(igraph)
library(dplyr)
library(stringr)
micro_pathway=read.delim("IV_res_bacteria_new_pathway_80%_wilcox_fdr_sig.txt",sep = "\t",header = T)
micro_pathway=micro_pathway %>% filter(FC>1.25)
data1=micro_pathway %>% select(bacteria,pathway,FC)
data1[which(data1$pathway=="peptidoglycan biosynthesis V (&beta;-lactam resistance)"),]$pathway="peptidoglycan biosynthesis V (beta-lactam resistance)"
data1$pathway=str_wrap(data1$pathway,width = 20)
data1$bacteria=str_wrap(data1$bacteria,width = 10)
graph1=graph_from_data_frame(data1)
plot(graph1)
node=as.data.frame(V(graph1))
V(graph1)$type=c(rep("Bacteria",6),rep("Pathway-Kegg",2),rep("Pathway-MetaCYC",6))
vcolor<-c("#f8d4b3","lightblue","tomato")
V(graph1)$color<-vcolor[factor(V(graph1)$type)]
E(graph1)$width=E(graph1)$FC*2.3
plot(graph1,
     layout=layout_nicely,
     edge.color="gray70",edge.arrow.size=.4, edge.curved=.1,
     vertex.label.font=1,
     vertex.label.color="black",
)

tkplot(graph1,
       layout=layout_nicely,
       edge.color="gray70",edge.arrow.size=.4, edge.curved=.1,
       vertex.label.font=1,
       vertex.label.color="black",
)

ll2=tkplot.getcoords(2)
load("network_plot_type.Rdata")
pdf("bacteria_pathway_network1.pdf",width = 9,height = 9)
plot(graph1,
     layout = ll2,
     edge.color="gray70",
     edge.arrow.size=1, 
     edge.curved=.1,
     
     vertex.label.font=1,
     vertex.label.dist=-2,
     vertex.label.color="black",)
legend(x=.8,y=1,
       levels(factor(V(graph1)$type)),pch=21,col="#777777",pt.bg=vcolor,cex = .8,
       pt.cex=2, bty="n", ncol=1)
dev.off()



######add spearman#####
micro_pathway=read.delim("IV_res_bacteria_new_pathway_80%_wilcox_fdr_sig.txt",sep = "\t",header = T)
micro_pathway=micro_pathway %>% filter(FC>1.25)
data1=micro_pathway %>% select(bacteria,pathway,FC)
data1[which(data1$pathway=="peptidoglycan biosynthesis V (&beta;-lactam resistance)"),]$pathway="peptidoglycan biosynthesis V (beta-lactam resistance)"
data1$pathway=str_wrap(data1$pathway,width = 20)
data1$bacteria=str_wrap(data1$bacteria,width = 10)
graph1=graph_from_data_frame(data1)
plot(graph1)
node=as.data.frame(V(graph1))
V(graph1)$type=c(rep("Bacteria",6),rep("Pathway-Kegg",2),rep("Pathway-MetaCYC",6))
vcolor<-c("#f8d4b3","lightblue","tomato")
V(graph1)$color<-vcolor[factor(V(graph1)$type)]
E(graph1)$width=E(graph1)$FC*2.3
plot(graph1,
     layout=layout_nicely,
     edge.color="gray70",edge.arrow.size=.4, edge.curved=.1,
     vertex.label.font=1,
     vertex.label.color="black",
)
legend(x=.8,y=1.5,
       levels(factor(V(graph1)$type)),pch=21,col="#777777",pt.bg=vcolor,cex = .8,
       pt.cex=2, bty="n", ncol=1)

###add eges###
cor1=read.csv("ivres_ba_cor_spearman_sig_result.csv")
cor1=cor1 %>% filter(cor!=1)
cor1=cor1[-6,]
cor1=cor1[,c(1:3)]
cor1$ba1=str_wrap(data1$bacteria,width = 10)
graph2=add.edges(graph1,c("Dialister\npneumosintes","Gemella\nmorbillorum",
                          "Dialister\npneumosintes","Parvimonas\nmicra",
                          "Dialister\npneumosintes","Peptostreptococcus\nporci",
                          "Dialister\npneumosintes","Peptostreptococcus\nstomatis",
                          "Dialister\npneumosintes","Solobacterium\nmoorei"),
                 color="#D6604D",width=cor1$cor*5,directed=F)

E(graph2)[is.na(E(graph2)$color)]$color="gray70"
tkplot(graph2,
       layout=ll2,edge.arrow.size=.4, edge.curved=.1,
       vertex.label.font=1,
       vertex.label.color="black",
)

ll2=tkplot.getcoords(2)
# save(ll2,file = "network_pathway_plot_type_1012.Rdata")
pdf("bacteria_pathway_network.pdf",width = 7,height = 7)
plot(graph2,
     layout=ll2,
     vertex.label.dist=-2,
     edge.arrow.size=.5, edge.curved=.1,
     vertex.label.font=1,
     vertex.label.color="black",
)
legend(x=1.2,y=.8,
       levels(factor(V(graph2)$type)),pch=21,col="#777777",pt.bg=vcolor,cex = .8,
       pt.cex=2, bty="n", ncol=1)
legend(
  x=1.2,y=1.,
  title = "Correlation (±)",
  legend = c("Positive"),
  col = c("#D6604D"),
  lty=1,
  lwd=2,
  bty="n",
  cex = .8,
  pt.cex=2
)
dev.off()
