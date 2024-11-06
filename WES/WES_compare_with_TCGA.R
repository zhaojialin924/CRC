####WES#####
####compare with TCGA#####
setwd("D:/Postgraduation/microbe_crc")

library(ggplot2)
library(tidyverse)
library(ggsci)
library(maftools)
library(data.table)
library(readxl)
library(sigminer)
library(NMF)
library(ggpubr)
library(tidyverse)
library(EnvStats)

####TMB####
maf_2=fread("Somatic_mutation.maf")
dim(maf_2)  #[1] 294189     55
maf_2$GnomAD_genome_AF_eas1=as.numeric(maf_2$GnomAD_genome_AF_eas)
maf_2$CADD1=as.numeric(maf_2$CADD)
maf_2$`1000g2015aug_all1`=as.numeric(maf_2$`1000g2015aug_all`)
maf_2$ExAC_ALL1=as.numeric(maf_2$ExAC_ALL)
maf_2$SIFT_score1=as.numeric(maf_2$SIFT_score)
maf_2$Polyphen2_HVAR_score1=as.numeric(maf_2$Polyphen2_HVAR_score)
maf_2$Polyphen2_HDIV_score1=as.numeric(maf_2$Polyphen2_HDIV_score)


maf_2_1=maf_2  %>%  
  filter (GnomAD_genome_AF_eas1 < 0.01 | is.na(GnomAD_genome_AF_eas1)) %>%
  filter (`1000g2015aug_all1` <0.01 | is.na(`1000g2015aug_all1`)) %>%
  filter (ExAC_ALL1 <0.01 | is.na(ExAC_ALL1)) %>%
  filter (CADD1 >=20 | is.na(CADD1)) %>%
  filter (SIFT_score1 <0.05 | is.na(SIFT_score1)) %>%
  filter (Polyphen2_HVAR_score1 >0.908 | is.na(Polyphen2_HVAR_score1)) %>%
  filter (Polyphen2_HDIV_score1 >0.908 | is.na(Polyphen2_HDIV_score1)) 
dim(maf_2_1)  #[1] 209135     55


table(maf_2_1$Variant_Classification)
finalmaf = maf_2_1 %>%
  filter (Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", 
                                        "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation",
                                        "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site"))
table(finalmaf$Variant_Classification)
save(finalmaf,file = "tmb_finalmaf.Rdata")
load("tmb_finalmaf.Rdata")
dim(finalmaf) #[1] 20995    55
library(maftools)
maf <- read.maf(maf=finalma)


pdf("101sample_TCGA_compare_QC1.pdf",width = 12,height =8 )
mutload <- tcgaCompare(maf = maf, cohortName = 'SYSUCC-CRC',tcga_capture_size = 40, 
                       capture_size = 40,
                       bg_col=c("white","#E7DAD2"),
                       medianCol = "#FA7F6F")
dev.off()
x=mutload$mutation_burden_perSample %>%
  as.tibble() %>%
  dplyr::filter(cohort=="SYSUCC-CRC") %>%
  head()

ggplot(mutload$mutation_burden_perSample,aes(x=cohort,y=total_perMB))+
  #geom_point(size=2)+
  geom_violin(aes(fill=cohort),show.legend = F,color="white")+
  geom_jitter(size=0.5,aes(color=cohort),show.legend = F)+
  stat_summary(fun="median",geom = "point",size=2,shape=21,color="black",fill="white")+
  scale_y_log10(breaks=c(0.1,1,10,100,1000),
                labels=c(0.1,1,10,100,1000))+
  theme_bw()+
  stat_n_text(size=2.3,y.pos = 3.2,vjust=0)+
  xlab("")+
  ylab("TMB(per MB)")+
  scale_fill_manual(values = c(rep("gray70",28),"#5f97d2","#ef7a6d",c(rep("gray70",4))))+
  scale_color_manual(values = c(rep("gray70",28),"#5f97d2","#ef7a6d",c(rep("gray70",4))))+
  theme(axis.text.x = element_text(angle = 45,hjust=1,
                                   color=c(rep("black",28),"#5f97d2","#ef7a6d",c(rep("black",4))),
                                   size=14),
        axis.text.y = element_text(color = "black",size=14),
        axis.title.y=element_text(color = "black",size = 16))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line())

ggsave("TMB_violin_SYSUCC22.pdf",width = 12,height = 8)