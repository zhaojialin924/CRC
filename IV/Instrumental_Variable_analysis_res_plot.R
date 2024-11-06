#####Instrumental Variable analysis####
####res plot####
####iv_res_small_plot###
setwd("D:/Postgraduation/microbe_crc")
library(ggplot2)
library(aplot)
library(readxl)
library(tidyverse)
library(patchwork)
library(cowplot)
library(grid)
library(ggplotify)
library(enrichplot)
clin=readxl::read_xlsx("01data/ours/clinical_infomation_20240413.xlsx")
activity=read.csv("02res/ours/WES_signature.csv",row.names = 1)
activity$sample=rownames(activity)
clin1=merge(clin1,activity)
###new###
data11=read.delim("count/species.filter.tsv",sep = "\t")
data11=as.data.frame(t(data11))
data11=data11/rowSums(data11)
data11=as.data.frame(t(data11))
data11$var=apply(data11,1,var)
data11=data11[order(-data11$var),]
top50=head(data11,50)
top50=top50 %>% select(-var)
top50hvb=as.data.frame(t(top50))
top50hvb$sample=rownames(top50hvb)
load("WES_mutation_gene_no_artifacts.Rdata")
mat=mat[,c(1:50)]
mat1 <- data.frame(sapply(mat, function(x) as.factor(as.character(x))))
mat1$sample=rownames(mat)
data1=merge(clin1,mat1,by="sample")


IV_res=read.csv("02res/ours/IV_res_sig.csv")
col=c("#E58579","#8AB1D2")
for(i in 1:nrow(IV_res)){
  gene=IV_res$Gene[i]
  microbe=IV_res$Bacteria[i]
  sig=IV_res$Sig[i]
  plotdf=data.frame(Siganture=data1[[sig]],Gene=data1[[gene]],Bacteria=data1[[microbe]])
  plotdf[which(plotdf$Bacteria==0),]$Bacteria=1e-4
  
  p1=ggplot(plotdf,aes(x=Siganture,y=Bacteria,color=as.factor(Gene),
                       shape=as.factor(Gene),fill=as.factor(Gene)))+
    geom_point(size=3,na.rm = T)+
    scale_color_manual(values = rev(col[-2]),
                       breaks=c("0", "1"),
                       labels=c("WT", "Mut"))+
    scale_fill_manual(values = rev(col[-2]),
                      breaks=c("0", "1"),
                      labels=c("WT", "Mut"))+
    scale_shape_manual(values = c(2,5),
                       breaks=c("0", "1"),
                       labels=c("WT", "Mut"))+
    scale_y_log10()+
    ylab("")+
    xlab("")+
    theme_bw()+
    labs(color=gene,shape=gene,fill=gene)
  
  p2=ggplot(plotdf,aes(x=as.factor(Gene),y=Bacteria,shape=as.factor(Gene),
                       color=as.factor(Gene)))+
    geom_boxplot(outlier.shape =NA)+
    ylab(paste0(microbe))+
    xlab("")+
    scale_color_manual(values = rev(col[-2]))+
    scale_shape_manual(values = c(2,5))+
    scale_y_log10() +
    scale_x_discrete(breaks=c("0","1"),
                     labels= c("WT","Mut"))+
    theme_bw()+
    geom_jitter(width = 0.2, size = 2) +
    theme(axis.text.y = element_text(size = 14, colour = "black"),
          axis.text.x = element_text(size = 14,colour="black"),
          axis.title.x = element_text(size = 16, colour = "black"),
          axis.title.y = element_text(size = 16, colour = "black"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5,size = 12))

  p3=ggplot(plotdf,aes(x=as.factor(Gene),y=Siganture,shape=as.factor(Gene),
                       color=as.factor(Gene)))+
    geom_boxplot(outlier.shape = NA)+
    ylab(ifelse(sig=="Sig1","SBS5","SBS6"))+
    xlab("")+
    scale_color_manual(values = rev(col[-2]))+
    scale_shape_manual(values = c(2,5))+
    scale_x_discrete(breaks=c("0","1"),
                     labels= c("WT","Mut"))+
    theme_bw()+
    geom_jitter(width = 0.2, size = 2) +
    theme(axis.text.y = element_text(size = 14, colour = "black",angle = 270,
                                     hjust = 0.5,vjust = 0.5),
          axis.text.x = element_text(size = 14,colour="black"),
          axis.title.x = element_text(size = 16, colour = "black"),
          axis.title.y = element_text(size = 16, colour = "black"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5,size = 12))+
    coord_flip()

  
  
  px=p1+theme(axis.ticks.y=element_blank(),
              axis.text = element_blank(),
              axis.ticks.x=element_blank(),
              legend.position = "top",
              # axis.text.x = element_text(size = 14,colour="black"),
              # axis.title.x = element_text(size = 16, colour = "black")
  )
  p <- px %>% 
    insert_bottom(p3,height = 0.4) %>%
    insert_left(p2,width=0.4) %>%
    as.ggplot()
  p
  ggsave(paste0("03picture/",i,"_",microbe,"_",gene,"_",sig,"_com_22.pdf"),
         width =7.2,height = 4.8)
  }
  
  




