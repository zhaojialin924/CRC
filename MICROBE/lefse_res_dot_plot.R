#####lefse result plot####
###lefse
setwd("D:/Postgraduation/microbe_crc/")
library(dplyr)
library(ggplot2)
col=c("#E58579","#f3a361","#8AB1D2")


####SBS5&SBS6#####
lefres=read.delim("02res/ours/micro_sig1_sig3_signature_group_lefse_count.res",sep = "\t",header = F)
lefres=lefres %>% subset(!is.na(V4))
colnames(lefres)=c("names","Logarith value","Group","LDA_value","P_value")
lefres$names=unlist(lapply(lefres$names,function(x) tail(strsplit(x,"\\.")[[1]],1)))

lefres=lefres[order(lefres$Group,lefres$LDA_value),]
lefres$P_value=as.numeric(lefres$P_value)
lefres$names=factor(lefres$names,levels = lefres$names)

ggplot(lefres,aes(x=-log10(P_value),y=names,fill=Group))+
  geom_point(aes(size=LDA_value,color=Group))+
  theme_bw()+
  labs(y="",size="LDA SCORE (log 10)")+
  theme(axis.text.x = element_text(colour = "black",size = 12),
        axis.text.y=element_text(size=12,color="black"))+
  scale_fill_manual(values = col)+
  scale_color_manual(values = col)

ggsave("lda_dot_count_all_sig1_sig3.png",width = 8,height = 6)
ggsave("lda_dot_count_all_sig1_sig3.pdf",width = 8,height = 6)

s_lefres=lefres[grepl("^s_",lefres$names),]

s_lefres$Group=ifelse(s_lefres$Group=="C1-Sig1","SBS5","SBS6")
s_lefres$Group=factor(s_lefres$Group,levels = c("SBS6","SBS5"))
ggplot(s_lefres,aes(x=-log10(P_value),y=str_wrap(names,width = 10),fill=Group))+
  geom_point(aes(size=LDA_value,color=Group))+
  theme_bw()+
  scale_size_continuous(range=c(3,8))+
  labs(y="",size="LDA SCORE (log 10)")+
  theme(axis.text.x = element_text(colour = "black",size = 14),
        axis.text.y=element_text(size=14,color="black"),
        axis.title = element_text(size = 16,color = "black"))+
  scale_fill_manual(values = rev(col[-2]))+
  scale_color_manual(values = rev(col[-2]))+
  scale_y_discrete(labels=function(x) str_wrap(x,width = 10))

ggsave("lda_dot_count_species_sig1_sig3_lda2_2.png",width = 8,height = 4)
ggsave("lda_dot_count_species_sig1_sig3_lda2_2.pdf",width = 8,height = 4)

