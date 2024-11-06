####eppk1 mutation status####
setwd("D:/Postgraduation/microbe_crc")
library(covTestR)
library(ggstatsplot)
library(ggplot2)
library(egg)
library(ggpubr)
library(Hmisc)
library(corrplot)
library(RColorBrewer)
col=c("#E58579","#f3a361","#8AB1D2")
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582","#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE","#4393C3", "#2166AC",
                           "#053061"))(100)
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
top50hvb=top50hvb[,c(51,1:50)]
res_sig=read.csv("1011_ivres_ba_cor_spearman_all_result.csv")
res_sig=res_sig[which(res_sig$fdr<0.05),]
res_sig=res_sig[which(res_sig$ba1=="Dialister pneumosintes" | res_sig$ba1=="Fusobacterium animalis"),]
top50hvb=top50hvb[,which(colnames(top50hvb) %in% c(unique(res_sig$ba2),"sample"))]
####EPPK1 MUT WT####

load("WES_mutation_gene_no_artifacts.Rdata")
mat=mat[,c(1:50)]
mat1 <- data.frame(sapply(mat, function(x) as.factor(as.character(x))))
mat1$sample=rownames(mat)
data1=merge(top50hvb,mat1[,c("sample","EPPK1")],by="sample")
WT=data1[which(data1$EPPK1==0),]
rownames(WT)=WT$sample
WT=WT[,-c(1,2,13)]
WT=WT[,sort(colnames(WT))]
mut=data1[which(data1$EPPK1==1),]
rownames(mut)=mut$sample
mut=mut[,-c(1,2,13)]
mut=mut[,sort(colnames(mut))]
mut=mut[-3,]
xx=rcorr(as.matrix(mut),type = "spearman")
mut_r=xx$r
p=xx$P
pdf("EPPK1_MUT_cor_heatmap.pdf",width = 6,height = 6)
corrplot(corr=mut_r,
         insig="label_sig",
         sig.level = c(.01, .05),
         pch.cex = 0.8,
         col = rev(col2),
         tl.col = "black")
dev.off()


xx=rcorr(as.matrix(WT),type = "spearman")
WT_r=xx$r
p=xx$P
pdf("eppk1_wt_cor_heatmap1.pdf",width = 6,height = 6)
corrplot(corr=WT_r,
         insig="label_sig",
         sig.level = c(.01, .05),
         pch.cex = 0.8,
         col = rev(col2),
         tl.col = "black")
dev.off()

data2=data1
for(ba1 in c("Dialister pneumosintes","Fusobacterium animalis","Parvimonas micra")){
  for (ba2 in unique(colnames(data2)[c(3:12)])){
    plotdf=data.frame(EPPK1=data2[["EPPK1"]],BA1=data2[[ba1]],BA2=data2[[ba2]])
    plotdf[which(plotdf$BA1==0),]$BA1=1e-4
    plotdf[which(plotdf$BA2==0),]$BA2=1e-4
    ggplot(plotdf,aes(BA1,BA2,
                      color=EPPK1))+
      geom_point(size=3,aes(color=EPPK1,shape=EPPK1))+ #"#a42c36"
      scale_x_log10()+
      scale_y_log10()+
      xlab(ba1)+
      ylab(ba2)+
      scale_color_manual(values = rev(col[-2]))+
      scale_shape_manual(values = c(2,5))+
      scale_color_manual(values = rev(col[-2]),
                         breaks=c("0", "1"),
                         labels=c("WT", "Mut"))+
      scale_fill_manual(values = rev(col[-2]),
                        breaks=c("0", "1"),
                        labels=c("WT", "Mut"))+
      scale_shape_manual(values = c(2,5),
                         breaks=c("0", "1"),
                         labels=c("WT", "Mut"))+
      theme_article()+
      theme(axis.text = element_text(size = 14,colour = "black"),
            axis.title.x = element_text(size = 16,colour = "black"),
            axis.title.y=element_text(size = 16,color = "black"))+
      geom_smooth(method = lm,fill="lightgray")+
      stat_cor(method = "spearman",cor.coef.name = "rho",)
    ggsave(paste0("eppk1_dotplot/1_eppk1_",ba1,"_",ba2,"_cor11111.pdf"),width = 4.6,height =3.6 )
    ggsave(paste0("eppk1_dotplot/1_eppk1_",ba1,"_",ba2,"_cor11111.png"),width = 4.6,height =3.6 )
  }
}

