#####beta diversityZ####
setwd("D:\\Postgraduation\\微生物")
library(tidyverse)
library(magrittr)
library(ggh4x)
library(ggsci)
library(ggplot2)
library(data.table)
library(vegan)
library(ggalt)
library(pairwiseAdonis)
library("ggplotify")
library(aplot)

data1=read.delim("count/species.tsv",sep = "\t",row.names = 1)
data1 <- data.frame(apply(data1,2,function(x) x/sum(x)),stringsAsFactors = FALSE)
otu.distance <- vegdist(t(data1),method = "bray", binary=F)
data1=as.data.frame(t(data1))

pcoa <- cmdscale (otu.distance,k=3,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)

pc12 <- as.data.frame(pc12)
pc12$sample <- row.names(pc12)
head(pc12)

sample1=readxl::read_xlsx("01data/ours/clinical_infomation_20240413.xlsx")
colnames(sample1)[1]="sample"
df1=merge(pc12,sample1)
act=read.csv("02res/ours/WES_signature.csv",row.names = 1)
act$sample=rownames(act)
df1=merge(df1,act)
df1$group=factor(df1$group,levels = c("SBS5","SBS10a","SBS6"))
rownames(df1)=df1$sample


df1=df1
dune.div <- adonis2(data1 ~ group,data = df1, permutations = 50000, method="bray")
dune.div

dune.pairwise.adonis <- pairwise.adonis(x=data1, factors=df1$group, sim.function = "vegdist",
                                        sim.method = "bray",
                                        p.adjust.m = "BH",
                                        reduce = NULL,
                                        perm = 999)

dune.pairwise.adonis
tab2 <- ggtexttable(dune.pairwise.adonis[,c("pairs","R2","p.value")], rows = NULL,
                    theme = ttheme("blank")) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 1)  %>%
  tab_add_hline(at.row = nrow(dune.pairwise.adonis)+1, row.side = "bottom", linewidth = 1)
tab2
ggsave("pcoa_signature_group_paire.png",width = 3,height = 2)
ggsave("pcoa_signature_group_paire.pdf",width = 3,height = 2)


dune_adonis <- paste0("PERMANOVA F-value: ",round(dune.div$F,2),"\n", 
                      "R2: ",round(dune.div$R2,2), "; P-value: ", dune.div$`Pr(>F)`)
col=rev(col)
p2<-ggplot(data=df1,aes(x=V1,y=V2,
                        color=group,shape=group))+
  theme_bw()+
  geom_point(size=3)+
  scale_shape_manual(values = c(17,15,18))+
  theme(panel.grid = element_blank())+
  scale_color_manual(values = col)+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=paste0("PCoA1 (",pc[1],"%)"),
       y=paste0("PCoA2 (",pc[2],"%)"),
  )+
  theme(axis.title.x=element_text(size=16,color="black"),
        axis.title.y=element_text(size=16,angle=90),
        axis.text.y=element_text(size=14,color="black"),
        axis.text.x=element_text(size=14,color="black"),
        panel.grid=element_blank())
p2
p3=p2+stat_ellipse(data=df1,geom = "polygon",
                   level=0.9,linetype = 2,size=0.5,
                   aes(fill=group),alpha=0.2,show.legend = T)+
  scale_fill_manual(values = col)
p3
pc1.density <-
  ggplot(df1) +
  geom_density(aes(x=V1, group=group,fill=group,linetype=group),
               color="black", alpha=0.6,position = 'identity',
               show.legend = F) +
  scale_fill_manual(values = col)+
  scale_y_discrete(expand = c(0,0.001))+
  labs(x=NULL,y=NULL)+
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())

pc2.density <-
  ggplot(df1) +
  geom_density(aes(x=V2, group=group, fill=group, linetype=group),
               color="black", alpha=0.6, position = 'identity',show.legend = F) +
  scale_fill_manual(values=col) +
  theme_classic()+
  scale_y_discrete(expand = c(0,0.001))+
  labs(x=NULL,y=NULL)+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  coord_flip()
p1 <- p3 %>% 
  insert_top(pc1.density,height = 0.3) %>% 
  insert_right(pc2.density,width=0.3) %>% 
  as.ggplot()
p1
ggsave("03picture/pcoa/pcoa_cluster4.pdf",width = 6.5,height = 4.8)


dune.div <- adonis2(data1 ~ Sex, data = df1, permutations = 5000, method="bray")
dune.div
dune_adonis <- paste0("PERMANOVA F-value: ",round(dune.div$F,2),"\n", 
                      "R2: ",round(dune.div$R2,2), "; P-value: ", dune.div$`Pr(>F)`)
p2<-ggplot(data=df1,aes(x=V1,y=V2,
                        color=Sex,shape=Sex))+
  theme_bw()+
  geom_point(size=3)+
  scale_color_manual(values = col[-2])+
  scale_shape_manual(values = c(17,18))+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=paste0("PCoA1 (",pc[1],"%)"),
       y=paste0("PCoA2 (",pc[2],"%)"),
  )+
  theme(axis.title.x=element_text(size=16,color="black"),
        axis.title.y=element_text(size=16,angle=90,color="black"),
        axis.text.y=element_text(size=14,color="black"),
        axis.text.x=element_text(size=14,color="black"),
        panel.grid=element_blank())
p3=p2+stat_ellipse(data=df1,geom = "polygon",
                   level=0.9,linetype = 2,size=0.5,
                   aes(fill=Sex),alpha=0.2,show.legend = T)+
  scale_fill_manual(values = col[-2])


pc1.density <-
  ggplot(df1) +
  geom_density(aes(x=V1, group=Sex,fill=Sex,linetype=Sex),
               color="black", alpha=0.6,position = 'identity',
               show.legend = F) +
  scale_fill_manual(values=col[-2]) +
  #scale_linetype_manual(values = c("solid","dashed"))+
  scale_y_discrete(expand = c(0,0.001))+
  labs(x=NULL,y=NULL)+
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())

pc2.density <-
  ggplot(df1) +
  geom_density(aes(x=V2, group=Sex, fill=Sex, linetype=Sex),
               color="black", alpha=0.6, position = 'identity',show.legend = F) +
  scale_fill_manual(values=col[-2]) +
  theme_classic()+
  scale_y_discrete(expand = c(0,0.001))+
  labs(x=NULL,y=NULL)+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  coord_flip()
p1 <- p3 %>% 
  insert_top(pc1.density,height = 0.3) %>% 
  insert_right(pc2.density,width=0.3) %>% 
  as.ggplot()
p1

ggsave("03picture/pcoa/pcoa_sex3.pdf",width = 6,height = 4.8)


dune.pairwise.adonis <- pairwise.adonis(x=data1, factors=df1$Tumor_Location_combine, sim.function = "vegdist",
                                        sim.method = "bray",
                                        p.adjust.m = "BH",
                                        reduce = NULL,
                                        perm = 999)

dune.pairwise.adonis
write.table(dune.pairwise.adonis,"pcoa_location_paire.txt",row.names = F,quote = F,sep = "\t")
tab2 <- ggtexttable(dune.pairwise.adonis[,c("pairs","p.value")], rows = NULL,
                    theme = ttheme("blank")) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 1)  %>%
  tab_add_hline(at.row = nrow(dune.pairwise.adonis)+1, row.side = "bottom", linewidth = 1)
tab2
ggsave("pcoa_location_paire.png",width = 4,height = 2)


dune.div <- adonis2(data1 ~ Tumor_Location_combine, data = df1, permutations = 5000, method="bray")
dune.div
dune_adonis <- paste0("PERMANOVA F-value: ",round(dune.div$F,2),"\n", 
                      "R2: ",round(dune.div$R2,2), "; P-value: ", dune.div$`Pr(>F)`)
p2<-ggplot(data=df1,aes(x=V1,y=V2,
                        color=Tumor_Location_combine,shape=Tumor_Location_combine))+
  theme_bw()+
  geom_point(size=3)+
  scale_color_manual(values = col[c(1,3,2)])+
  scale_shape_manual(values = c(17,18,15))+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=paste0("PCoA1 (",pc[1],"%)"),
       y=paste0("PCoA2 (",pc[2],"%)")
  )+
  theme(axis.title.x=element_text(size=16,color="black"),
        axis.title.y=element_text(size=16,angle=90,color="black"),
        axis.text.y=element_text(size=14,color="black"),
        axis.text.x=element_text(size=14,color="black"),
        panel.grid=element_blank())

p2+stat_ellipse(data=df1,geom = "polygon",
                level=0.9,linetype = 2,size=0.5,
                aes(fill=Tumor_Location_combine),alpha=0.2,show.legend = T)+
  scale_fill_manual(values = col[c(1,3,2)])
p3=p2+stat_ellipse(data=df1,geom = "polygon",
                   level=0.9,linetype = 2,size=0.5,
                   aes(fill=Tumor_Location_combine),alpha=0.2,show.legend = T)+
  scale_fill_manual(values = col[c(1,3,2)])+
  labs(fill="Tumor Location",
       color="Tumor Location",
       shape="Tumor Location")

pc1.density <-
  ggplot(df1) +
  geom_density(aes(x=V1, group=Tumor_Location_combine,fill=Tumor_Location_combine,linetype=Tumor_Location_combine),
               color="black", alpha=0.6,position = 'identity',
               show.legend = F) +
  scale_fill_manual(values = col[c(1,3,2)])+
  scale_y_discrete(expand = c(0,0.001))+
  labs(x=NULL,y=NULL)+
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())

pc2.density <-
  ggplot(df1) +
  geom_density(aes(x=V2, group=Tumor_Location_combine, fill=Tumor_Location_combine, linetype=Tumor_Location_combine),
               color="black", alpha=0.6, position = 'identity',show.legend = F) +
  scale_fill_manual(values = col[c(1,3,2)])+
  theme_classic()+
  scale_y_discrete(expand = c(0,0.001))+
  labs(x=NULL,y=NULL)+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  coord_flip()
p1 <- p3 %>% 
  insert_top(pc1.density,height = 0.3) %>% 
  insert_right(pc2.density,width=0.3) %>% 
  as.ggplot()
p1
ggsave("03picture/pcoa/pcoa_Tumor_Location_combine33.pdf",width = 6.7,height = 4.8)




dune.pairwise.adonis <- pairwise.adonis(x=data1, factors=df1$stage1, sim.function = "vegdist",
                                        sim.method = "bray",
                                        p.adjust.m = "BH",
                                        reduce = NULL,
                                        perm = 999)

dune.pairwise.adonis
write.table(dune.pairwise.adonis,"pcoa_stage_paire.txt",row.names = F,quote = F,sep = "\t")
tab2 <- ggtexttable(dune.pairwise.adonis[,c("pairs","p.value")], rows = NULL, 
                    theme = ttheme("blank")) %>% 
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 1)  %>% 
  tab_add_hline(at.row = nrow(dune.pairwise.adonis)+1, row.side = "bottom", linewidth = 1) 
tab2
ggsave("pcoa_stage_paire.png",width = 3.2,height = 2)
ggsave("pcoa_stage_paire.pdf",width = 3.2,height = 2)

dune.div <- adonis2(data1 ~ stage1, data = df1, permutations = 999, method="bray")
dune.div
dune_adonis <- paste0("PERMANOVA F-value: ",round(dune.div$F,2),"\n", 
                      "R2: ",round(dune.div$R2,2), "; P-value: ", dune.div$`Pr(>F)`)
p2<-ggplot(data=df1,aes(x=V1,y=V2,
                        color=stage1,shape=stage1))+
  theme_bw()+
  geom_point(size=3)+
  scale_color_manual(values = col)+
  scale_shape_manual(values = c(17,15,18))+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=paste0("PCoA1 (",pc[1],"%)"),
       y=paste0("PCoA2 (",pc[2],"%)")
  )+
  theme(axis.title.x=element_text(size=16,color="black"),
        axis.title.y=element_text(size=16,angle=90,color="black"),
        axis.text.y=element_text(size=14,color="black"),
        axis.text.x=element_text(size=14,color="black"),
        panel.grid=element_blank())
p2+stat_ellipse(data=df1,geom = "polygon",
                level=0.9,linetype = 2,size=0.5,
                aes(fill=stage1),alpha=0.2,show.legend = T)
p3=p2+stat_ellipse(data=df1,geom = "polygon",
                   level=0.9,linetype = 2,size=0.5,
                   aes(fill=stage1),alpha=0.2,show.legend = T)+
  scale_fill_manual(values = col)+
  labs(fill="Stage",
       color="Stage",
       shape="Stage")


pc1.density <-
  ggplot(df1) +
  geom_density(aes(x=V1, group=stage1,fill=stage1,linetype=stage1),
               color="black", alpha=0.6,position = 'identity',
               show.legend = F) +
  scale_fill_manual(values = col)+
  scale_y_discrete(expand = c(0,0.001))+
  labs(x=NULL,y=NULL)+
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())



pc2.density <-
  ggplot(df1) +
  geom_density(aes(x=V2, group=stage1, fill=stage1, linetype=stage1),
               color="black", alpha=0.6, position = 'identity',show.legend = F) +
  scale_fill_manual(values = col)+
  theme_classic()+
  scale_y_discrete(expand = c(0,0.001))+
  labs(x=NULL,y=NULL)+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  coord_flip()
p1 <- p3 %>% 
  insert_top(pc1.density,height = 0.3) %>% 
  insert_right(pc2.density,width=0.3) %>% 
  as.ggplot()
p1
ggsave("03picture/pcoa/pcoa_stage3.png",width = 5.6,height = 4.8)
ggsave("03picture/pcoa/pcoa_stage3.pdf",width = 5.6,height = 4.8)
