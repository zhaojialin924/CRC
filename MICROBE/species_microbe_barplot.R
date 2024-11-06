library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(ggsci)
library(ggalluvial)
setwd("D:/Postgraduation/microbe_crc")
col=pal_d3("category20")(20) 
col2 = pal_d3("category20b",alpha = 0.7)(20) 
mypal=c(col[c(1:7,9:17,19:20)],col2[c(1:7)],"gray")
data=read.delim("count/species.tsv",sep = "\t")
rownames(data)=data$species
data=data[,-1]
data=as.data.frame(t(data))
data1=data/rowSums(data)
data1=as.data.frame(t(data1))
data1$species=rownames(data1)
data1 <- data1 %>%
  reshape2::melt(id=c('species')) %>%
  mutate(sample=variable) %>%
  dplyr::select(sample,species,abundance=value)

plot_dat=data1 %>%
  mutate(Species=ifelse(species %in% top_species$species,species,'Other')) %>%
  mutate(abundance=100*abundance)

plot_dat$Species=factor(plot_dat$Species,levels = c("Other",rev(top_species$species)))
sample=unique(plot_dat$sample)
top_abund=plot_dat[which(plot_dat$species==top_species$species[1]),]
top_abund=top_abund %>%
  arrange(abundance)
plot_dat$sample=factor(plot_dat$sample,levels = top_abund$sample)

plot_dat1=plot_dat %>%
  select(-species) %>%
  group_by(sample,Species) %>%
  mutate(abundance1=sum(abundance)) %>%
  select(-abundance) %>%
  distinct()
# save(plot_dat1,file = "Species_barplot_flow_plotdat.Rdata")
###plot####
ggplot(plot_dat1,aes(sample,abundance1,fill=Species,
                     stratum=Species,alluvium=Species))+
  geom_stratum(width = 0.5, color='white')+
  geom_alluvium(alpha = 0.5,  
                width = 0.5,
                color='white',
                linewidth = 0.5,
                curve_type = "xspline")+
  theme(axis.text=element_text(size=16,color='black'),
        axis.title=element_text(size=16),legend.position='top',
        axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black",size = 0.5),
        legend.text = element_text(size=16,color="black"),
        legend.title = element_text(size=16,color="black"))+
  xlab("Sample")+ylab("Relative abundance (%)")+
  scale_fill_manual(values = rev(mypal))+
  guides(fill=guide_legend(reverse = T))
ggsave("Species_all_barplot_flow2.pdf",width = 20,height = 10 )

#####add siganture group####
activity=read.csv("02res/ours/WES_signature.csv",row.names = 1)
activity$Sig_group=paste0("Sig",max.col(activity))
activity$sample=rownames(activity)
activity=activity[,c(4,5)]
plot_dat12=merge(plot_dat1,activity,by="sample")
# save(plot_dat12,file = "Species_barplot_flow_siganture_plotdat.Rdata")
ggplot(plot_dat12,aes(sample,abundance1,fill=Species,
                      stratum=Species,alluvium=Species))+
  geom_stratum(width = 0.5, color='white')+
  geom_alluvium(alpha = 0.5,  
                width = 0.5,
                color='white',
                linewidth = 0.5,
                curve_type = "xspline")+
  theme(axis.text=element_text(size=16,color='black'),
        axis.title=element_text(size=16),legend.position='top',
        axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black",size = 0.5),
        legend.text = element_text(size=16,color="black"),
        legend.title = element_text(size=16,color="black"))+
  xlab("Sample")+ylab("Relative abundance (%)")+
  scale_fill_manual(values = rev(mypal))+
  guides(fill=guide_legend(reverse = T))+
  facet_grid(~Sig_group,scales = "free", space = "free")
ggsave("Species_all_barplot_flow_sig_group1.pdf",width = 20,height = 10 )



####

activity=read.csv("02res/ours/WES_signature.csv",row.names = 1)
activity$Sig_group=paste0("C",max.col(activity),"-Sig",max.col(activity))
activity$sample=rownames(activity)
activity$Group=ifelse(activity$Sig_group=="C1-Sig1","SBS5",
                      ifelse(activity$Sig_group=="C2-Sig2","SBS10a","SBS6"))
data=read.delim("count/species.tsv",sep = "\t")
rownames(data)=data$species
data=data[,-1]
data=as.data.frame(t(data))
data1=data/rowSums(data)
data1=as.data.frame(t(data1))
data1$species=rownames(data1)
data1 <- data1 %>%
  reshape2::melt(id=c('species')) %>%
  mutate(sample=variable) %>%
  dplyr::select(sample,species,abundance=value)

top_species <- data1 %>%
  group_by(species) %>%
  summarise(abundance1=mean(abundance)) %>%
  arrange(desc(abundance1)) %>%
  ungroup() %>%
  dplyr::filter(!duplicated(species)) %>%
  head(25)
data2=merge(data1,top_species,by="species")
top_species1 <- data3 %>%
  select(-sample,Sig1,Sig2,Sig3,Sig_group) %>%
  group_by(Group,species) %>%
  summarise(abundance1=mean(abundance)) %>%
  arrange(desc(abundance1)) %>%
  ungroup() %>%
  distinct()

sig_ave=reshape2::dcast(top_species1,species~Group,value.var = "abundance1")
colnames(top_species)[2]="Overall"
data2=merge(top_species,sig_ave,by="species")

data3=reshape2::melt(data2)
data3$species=factor(data3$species,levels = rev(top_species$species))
data3$variable=factor(data3$variable,levels = rev(c("Overall","SBS5",
                                                    "SBS10a","SBS6")))

ggplot(data3,aes(x=variable,y=value,fill=species))+
  geom_bar(stat='identity',position="fill",width = 0.6)+
  theme(axis.text=element_text(size=12,color='black'),
        axis.title=element_text(size=16),
        legend.position='bottom',
        axis.text.x = element_text(size = 12,colour = "black"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black",size = 0.5),
        legend.text = element_text(size=10,color="black"),
        legend.title = element_text(size=12,color="black"))+
  xlab("")+ylab("Percentage")+
  scale_fill_manual(values = rev(mypal[-26]))+
  guides(fill=guide_legend(reverse = T))+
  scale_y_continuous(expand = c(0,0))+
  labs(fill="Species")+
  coord_flip()
ggsave("species_ave_sbs_t.pdf",width = 12,height = 6)