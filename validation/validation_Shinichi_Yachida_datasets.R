####Shinichi Yachida et al. dataset #####
library(igraph)
library(readxl)
library(ggplot2)
library(Hmisc)
setwd("D:/Postgraduation/microbe_crc")

clin=readxl::read_xlsx("01data/SY_dataset/41591_2019_458_MOESM3_ESM (1).xlsx",sheet = "Table_S2-1",
                       skip = 2)
data1=readxl::read_xlsx("01data/SY_dataset/41591_2019_458_MOESM3_ESM (1).xlsx",sheet = "Table_S7-1",
                        skip = 3)

IVres=read.csv("02res/ours/IV_res_sig.csv")
res_sig=read.csv("ivres_ba_cor_spearman_all_result.csv")
res_sig=res_sig[which(res_sig$fdr<0.05),]
res_sig=res_sig[which(res_sig$ba1=="Dialister pneumosintes" | res_sig$ba1=="Fusobacterium animalis"),]

ba2=unique(c(res_sig$ba2))

ba2[4]="Fusobacterium nucleatum subsp. animalis"
ba2[6]="Fusobacterium nucleatum subsp. polymorphum"
ba2[5]="Fusobacterium nucleatum subsp. nucleatum"
data2=data1[which(data1$...1 %in% ba2),]
data2=as.data.frame(data2)
rownames(data2)=data2$...1
data2=data2[,-1]
normal=clin[which(clin$Group=="Healthy"),]
case=clin[which(clin$Group %in% c("Stage_I_II","Stage_III_IV")),]

normal_data=data2[,which(colnames(data2) %in% normal$Subject_ID)]
case_data=data2[,which(colnames(data2) %in% case$Subject_ID)]

normal=rcorr(as.matrix(t(normal_data)),type = "spearman")

r=normal$r
p=normal$P

edges <- as.data.frame(as.table(r))
p1=as.data.frame(as.table(p))
colnames(p1)=c("from","to","pvalue")
colnames(edges) <- c("from", "to", "correlation")
edges=merge(edges,p1,by=c("from","to"))
edges=edges[which(edges$from=="Dialister pneumosintes"| edges$from=="Fusobacterium nucleatum subsp. animalis"),]
edges=edges[which(edges$correlation<1),]
edges=edges[-2,]

edges$from=str_wrap(edges$from,width = 10)
edges$to=str_wrap(edges$to,width = 10)
edges=edges[which(edges$pvalue<0.05),]
edges$color=ifelse(edges$correlation>0,"#D6604D","#4393C3")
edges$weight <- edges$correlation  
vertices <- data.frame(name = str_wrap(colnames(r),width = 10))
g <- graph_from_data_frame(edges, vertices = vertices, directed = FALSE)
plot(g, 
     layout=layout.circle,
     vertex.label.cex=1,
     edge.width = abs(E(g)$weight*10), 
     vertex.color = "skyblue", 
     vertex.size = 20,
     vertex.label.dist=-2,
     vertex.label.color="black",)
tkplot(g, 
       layout=layout.circle,
       vertex.label.cex=1,
       edge.width = abs(E(g)$weight*10),
       vertex.color = "skyblue", 
       vertex.size = 20,
       vertex.label.dist=-2,
       vertex.label.color="black",)
ll4=tkplot.getcoords(1)

pdf("sydataset_normal_dp_fa_cor_network2.pdf",width = 6,height = 6)
plot(g, 
     layout=ll4,
     vertex.label.cex=1.5,
     # vertex.label.col="black",
     edge.width = abs(E(g)$weight*12), 
     vertex.color = c("#f8d4b3","#f8d4b3","#D6604D",rep("#f8d4b3",4),"#D6604D","#f8d4b3","#f8d4b3"),  
     vertex.size = 20,
     vertex.label.dist=-2,
     vertex.label.color="black",)

legend(
  x=.9,y=1.,
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




case=rcorr(as.matrix(t(case_data)),type = "spearman")
r=case$r
p=case$P

edges <- as.data.frame(as.table(r))
p1=as.data.frame(as.table(p))
colnames(p1)=c("from","to","pvalue")
colnames(edges) <- c("from", "to", "correlation")
edges=merge(edges,p1,by=c("from","to"))
edges=edges[which(edges$from=="Dialister pneumosintes"| edges$from=="Fusobacterium nucleatum subsp. animalis"),]
edges=edges[which(edges$correlation<1),]
edges=edges[-2,]
write.table(edges,"sydataset_case_DP_Cluster_microbe_cor_1101.txt",sep = "\t",
            row.names = F,quote = F)
edges$from=str_wrap(edges$from,width = 10)
edges$to=str_wrap(edges$to,width = 10)
edges=edges[which(edges$pvalue<0.05),]
edges$lty=ifelse(edges$pvalue<0.05,1,2)
edges$color=ifelse(edges$correlation>0,"#D6604D","#4393C3")
edges$weight <- edges$correlation  
vertices <- data.frame(name = str_wrap(colnames(r),width = 10))
g <- graph_from_data_frame(edges, vertices = vertices, directed = FALSE)

pdf("sydataset_case_dp_fa_cor_network2.pdf",width = 6,height = 6)
plot(g, 
     layout=ll4,
     vertex.label.cex=1.5,
     # vertex.label.col="black",
     edge.width = abs(E(g)$weight*12), 
     vertex.color = c("#f8d4b3","#f8d4b3","#D6604D",rep("#f8d4b3",4),"#D6604D","#f8d4b3","#f8d4b3"), 
     vertex.size = 20,
     vertex.label.dist=-2,
     vertex.label.color="black",)

legend(
  x=.9,y=1.,
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

