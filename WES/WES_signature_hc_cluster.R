#####wes###
####signature###
setwd("D:/Postgraduation/microbe_crc")

library(ggplot2)
library(tidyverse)
library(ggsci)
library(maftools)
library(data.table)
library(readxl)
library(sigminer)
library(NMF)
library(data.table)
library(maftools)
library(pheatmap)
library(factoextra)
library(igraph)

####Signature####
maf_2=fread("01data/ours/ying-analysis/Somatic_mutation.maf")
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
  filter (ExAC_ALL1 <0.01 | is.na(ExAC_ALL1)) 

save(maf_2_1,file = "siganture_finalmaf.Rdata")
load("siganture_finalmaf.Rdata")

maf <- read.maf(maf=maf_2_1)

mats <- mt_tally <- sig_tally(
  maf,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  useSyn = TRUE,
  mode = "ALL"
)

est <- sig_estimate(mt_tally$SBS_96, range = 2:5, nrun = 50, seed = 123456,verbose = TRUE)

show_sig_number_survey2(est$survey)

sigs <- sig_extract(mt_tally$SBS_96, n_sig = 3, nrun = 50, seed = 123456)
pic <- show_sig_profile(sigs, mode = "SBS", style = "cosmic")
sim <- get_sig_similarity(sigs, sig_db = "SBS")
add_labels(pic, x = 0.72, y = 0.25, y_end = 0.9, labels = sim, n_label = 3)
show_sig_exposure(sigs, rm_space = TRUE, style = "cosmic")

activity <- get_sig_exposure(sigs,type = "relative")
head(activity)
write.csv(activity,"02res/ours/WES_signature.csv",row.names=F,quote=F)
activity=read.csv("02res/ours/WES_signature.csv")
activity=as.data.frame(activity)
rownames(activity)=activity$sample
activity=activity[,-1]

##caculate distance
d<-dist(activity)
fit1<-hclust(d,method = "average")
plot(fit1,hang = -1,cex=.8)

pdf("03picture/p1/101sample_euclidean_distance_signature_heatmap.pdf",width = 20,height = 20)
heatmap(as.matrix(d))
dev.off()


fviz_dend(fit1,k=3)
fviz_dend(fit1,k=3,rect =T,cex=0.5,
          k_colors = c("#2E9FDF", "#00AFBB","#FC4E07"))
ggsave(filename = "03picture/p1/HC_signature_101sample.png",width = 10,height = 6)
