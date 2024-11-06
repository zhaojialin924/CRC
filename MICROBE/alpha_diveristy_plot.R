####alpha diveersity####
library(vegan)
library(ggpubr)
library(patchwork)
library(gridExtra)
library(egg)
setwd("D:/Postgraduation/microbe_crc")


tax=c("species")
calculate_diversity <- function(otu) { 
  # Richness (number of species)
  richness <- estimateR(t(otu))[1, ]  
  # ACE (Abundance-based Coverage Estimator)
  ace <- estimateR(t(otu))[4, ] 
  # Chao
  chao <- estimateR(t(otu))[2, ]  
  # Shannon diversity index
  shannon <- diversity(t(otu), index = "shannon")  
  # Simpson diversity index
  simpson <- diversity(t(otu), index = "simpson")  
  # Pielou's evenness index
  pielou <- diversity(t(otu), index = "shannon") / log(specnumber(t(otu)))  
  # Inverse Simpson index
  invsimpson <- 1 / simpson  
  # Combine results into a dataframe
  results_df <- data.frame(Richness = richness, ACE = ace, Chao = chao, 
                           Shannon = shannon, Simpson = simpson, 
                           Pielou = pielou, Invsimpson = invsimpson)  
  return(results_df)
}

lapply(1:length(tax), function(i) {
  j=tax[i]
  data = final %>% select(j,9:92) 
  data = data %>% group_by(across(all_of(j))) %>% summarise(across(everything(), sum), .groups = 'drop')
  write.table(data,paste("count/",paste(j,".tsv",sep=""),sep=""),row.names=F,col.names=T,sep="\t",quote=F)
  data = data.frame(data)
  row.names(data)=data[,1]
  data = data[,-1]
  norm =as.data.frame(t(rrarefy(t(data),min(colSums(data)))))
  write.table(norm,paste("norm/",paste(j,".norm.tsv",sep=""),sep=""),row.names=T,col.names=T,sep="\t",quote=F)
  alpha=calculate_diversity(norm)
  write.table(alpha,paste("alpha/",paste(j,"alpha.tsv",sep="."),sep=""),row.names=T,col.names=T,sep="\t",quote=F)
})


####plot#####
alpha11=read.delim("D:\\Postgraduation\\microbe_crc\\alpha\\species.alpha.tsv",row.names = 1)
alpha11$sample=rownames(alpha11)


####siganture#####
activity=read.csv("02res/ours/WES_signature.csv",row.names = 1)
activity$sample=rownames(activity)
alpha11_act=merge(alpha11,activity)
alpha11_act$signature=paste0("Sig",max.col(alpha11_act[,c(9:11)]))
alpha11_act[which(alpha11_act$signature=="Sig1"),]$group="SBS5"
alpha11_act[which(alpha11_act$signature=="Sig2"),]$group="SBS10a"
alpha11_act[which(alpha11_act$signature=="Sig3"),]$group="SBS6"
alpha11_act$group=factor(alpha11_act$group,levels = c("SBS5","SBS10a","SBS6"))

col=c("#E58579","#ffe5c3","#8AB1D2")
col=rev(col)
p1=ggplot(alpha11_act,aes(x=group,y=Shannon,color=group,shape=group))+
  geom_boxplot()+
  geom_jitter(aes(fill=group),position = position_jitter(width = 0.2))+
  ylab("Shannon")+
  xlab("")+
  theme_bw()+
  scale_color_manual(values = col)+
  scale_fill_manual(values=col)+
  scale_shape_manual(values = c(17,15,18))+
  theme(axis.text = element_text(color = "black"),
        axis.text.x=element_text(angle = 45,hjust = 1,size=14),
        axis.text.y=element_text(size=14),
        axis.title.y = element_text(size=16),
        legend.position = "none"
  )+
  stat_compare_means(#aes(group =  group),
    comparisons = list(c("SBS5","SBS10a"),
                       c("SBS5","SBS6"),
                       c("SBS10a","SBS6")),
    label = "p.signif",
    label.x.npc = "middle",
    method = "wilcox.test",
    #hide.ns = T
  )

p1
p2=ggplot(alpha11_act,aes(x=group,y=Simpson,color=group,shape=group))+
  geom_boxplot()+
  geom_jitter(aes(fill=group),position = position_jitter(width = 0.2))+
  scale_color_manual(values = col)+
  scale_fill_manual(values = col)+
  scale_shape_manual(values = c(17,15,18))+
  ylab("Simpson")+
  xlab("")+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        axis.text.x=element_text(angle = 45,hjust = 1,size=14),
        axis.text.y=element_text(size=14),
        axis.title.y = element_text(size=16),
        legend.position = "none"
        #legend.position='top'
  )+
  stat_compare_means(#aes(group =  group),
    comparisons = list(c("SBS5","SBS10a"),
                       c("SBS5","SBS6"),
                       c("SBS10a","SBS6")),
    label = "p.signif",
    label.x.npc = "middle",
    method = "wilcox.test",
    #hide.ns = T
  )

p2

ggsave("signature_shannon_alpha_boxplot3.pdf",width = 4,height = 4,
       egg::set_panel_size(p1,width=unit(1.8, "in"), height=unit(3, "in")))
ggsave("signature_simpson_alpha_boxplot3.pdf",width = 4,height = 4,
       egg::set_panel_size(p2,width=unit(1.8, "in"), height=unit(3, "in")))



sample1=readxl::read_xlsx("01data/ours/clinical_infomation_20240413.xlsx")
alpha11_act_clinical=merge(alpha11_act,sample1)


###stage###
col=rev(col)

p2=ggplot(alpha11_act_clinical,aes(x=stage1,y=Shannon,color=stage1,shape=stage1))+
  geom_boxplot()+
  geom_jitter(aes(fill=stage1),position = position_jitter(width = 0.2))+
  scale_color_manual(values = col)+
  scale_fill_manual(values = col)+
  scale_shape_manual(values = c(17,15,18))+
  ylab("Shannon")+
  xlab("")+
  theme_bw()+
  labs(fill="Stage",
       color="Stage")+
  theme(axis.text = element_text(color = "black"),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        axis.title.y = element_text(size=16),
        legend.position = "none"
        #legend.position='top'
  )+
  stat_compare_means(#aes(group =  group),
    comparisons = list(c("I","II"),
                       c("I","III"),
                       c("II","III")),
    label = "p.signif",
    label.x.npc = "middle",
    method = "wilcox.test",
    #hide.ns = T
  )


p3=ggplot(alpha11_act_clinical,aes(x=stage1,y=Simpson,color=stage1,shape=stage1))+
  geom_boxplot()+
  geom_jitter(aes(fill=stage1),position = position_jitter(width = 0.2))+
  scale_color_manual(values = col)+
  scale_fill_manual(values = col)+
  scale_shape_manual(values = c(17,15,18))+
  ylab("Simpson")+
  xlab("")+
  theme_bw()+
  labs(fill="Stage",
       color="Stage")+
  theme(axis.text = element_text(color = "black"),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        axis.title.y = element_text(size=16),
        legend.position = "none"
        #legend.position='top'
  )+
  stat_compare_means(#aes(group =  group),
    comparisons = list(c("I","II"),
                       c("I","III"),
                       c("II","III")),
    label = "p.signif",
    label.x.npc = "middle",
    method = "wilcox.test"
  )


ggsave("stage_shannon_alpha_boxplot3.pdf",width = 4,height = 4,
       egg::set_panel_size(p2,width=unit(1.8, "in"), height=unit(3, "in")))
ggsave("stage_simpson_alpha_boxplot3.pdf",width = 4,height = 4,
       egg::set_panel_size(p3,width=unit(1.8, "in"), height=unit(3, "in")))



###Sex###

p2=ggplot(alpha11_act_clinical,aes(x=Sex,y=Shannon,color=Sex,shape=Sex))+
  geom_boxplot()+
  geom_jitter(aes(fill=Sex),position = position_jitter(width = 0.2))+
  scale_color_manual(values = col[-2])+
  scale_fill_manual(values = col[-2])+
  scale_shape_manual(values = c(17,18))+
  ylab("Shannon")+
  xlab("")+
  theme_bw()+
  labs(fill="Sex",
       color="Sex")+
  theme(axis.text = element_text(color = "black"),
        axis.text.x=element_text(angle = 45,hjust = 1,size=14),
        axis.text.y=element_text(size=14),
        axis.title.y = element_text(size=16),
        legend.position = "none"
        #legend.position='top'
  )+
  stat_compare_means(#aes(group =  group),
    comparisons = list(c("Female","Male")),
    label = "p.signif",
    label.x.npc = "middle",
    method = "wilcox.test",
    #hide.ns = T
  )

p3=ggplot(alpha11_act_clinical,aes(x=Sex,y=Simpson,color=Sex,shape=Sex))+
  geom_boxplot()+
  geom_jitter(aes(fill=Sex),position = position_jitter(width = 0.2))+
  scale_color_manual(values = col[-2])+
  scale_fill_manual(values = col[-2])+
  scale_shape_manual(values = c(17,18))+
  ylab("Simpson")+
  xlab("")+
  theme_bw()+
  labs(fill="Sex",
       color="Sex")+
  theme(axis.text = element_text(color = "black"),
        axis.text.x=element_text(angle = 45,hjust = 1,size=14),
        axis.text.y=element_text(size=14),
        axis.title.y = element_text(size=16),
        legend.position = "none"
        #legend.position='top'
  )+
  stat_compare_means(#aes(group =  group),
    comparisons = list(c("Female","Male")),
    label = "p.signif",
    label.x.npc = "middle",
    method = "wilcox.test",
    #hide.ns = T
  )

ggsave("sex_shannon_alpha_boxplot3.pdf",width = 4,height = 4,
       egg::set_panel_size(p2,width=unit(1.2, "in"), height=unit(3, "in")))
ggsave("sex_simpson_alpha_boxplot3.pdf",width = 4,height = 4,
       egg::set_panel_size(p3,width=unit(1.2, "in"), height=unit(3, "in")))


####location####
alpha11_act_clinical$Tumor_Location_combine=factor(alpha11_act_clinical$Tumor_Location_combine,
                                                   levels = c("Left hemicolon","Right hemicolon","Rectum"))

p2=ggplot(alpha11_act_clinical,aes(x=Tumor_Location_combine,y=Shannon,color=Tumor_Location_combine,
                                   shape=Tumor_Location_combine))+
  geom_boxplot()+
  geom_jitter(aes(fill=Tumor_Location_combine),position = position_jitter(width = 0.2))+
  scale_color_manual(values = col[c(1,3,2)])+
  scale_fill_manual(values = col[c(1,3,2)])+
  scale_shape_manual(values = c(17,18,15))+
  ylab("Shannon")+
  xlab("")+
  theme_bw()+
  labs(fill="Tumor Location",
       color="Tumor Location")+
  theme(axis.text = element_text(color = "black"),
        axis.text.x=element_text(angle = 45,hjust = 1,size=14),
        axis.text.y=element_text(size=14),
        axis.title.y = element_text(size=16),
        legend.position = "none"
        #legend.position='top'
  )+
  stat_compare_means(#aes(group =  group),
    comparisons = list(c("Left hemicolon","Rectum"),
                       c("Left hemicolon","Right hemicolon"),
                       c("Rectum","Right hemicolon")),
    label = "p.signif",
    label.x.npc = "middle",
    method = "wilcox.test",
    #hide.ns = T
  )

p3=ggplot(alpha11_act_clinical,aes(x=Tumor_Location_combine,y=Simpson,
                                   color=Tumor_Location_combine,
                                   shape=Tumor_Location_combine))+
  geom_boxplot()+
  geom_jitter(aes(fill=Tumor_Location_combine),
              position = position_jitter(width = 0.2))+
  scale_color_manual(values = col[c(1,3,2)])+
  scale_fill_manual(values = col[c(1,3,2)])+
  scale_shape_manual(values = c(17,18,15))+
  ylab("Simpson")+
  xlab("")+
  theme_bw()+
  labs(fill="Tumor Location",
       color="Tumor Location")+
  theme(axis.text = element_text(color = "black"),
        axis.text.x=element_text(angle = 45,hjust = 1,size=14),
        axis.text.y=element_text(size=14),
        axis.title.y = element_text(size=16),
        legend.position = "none"
        #legend.position='top'
  )+
  stat_compare_means(#aes(group =  group),
    comparisons = list(c("Left hemicolon","Rectum"),
                       c("Left hemicolon","Right hemicolon"),
                       c("Rectum","Right hemicolon")),
    label = "p.signif",
    label.x.npc = "middle",
    method = "wilcox.test",
    #hide.ns = T
  )


ggsave("location_shannon_alpha_boxplot3.pdf",width = 6,height = 6,
       egg::set_panel_size(p2,width=unit(1.8, "in"), height=unit(3, "in")))
ggsave("location_simpson_alpha_boxplot3.pdf",width = 6,height = 6,
       egg::set_panel_size(p3,width=unit(1.8, "in"), height=unit(3, "in")))

