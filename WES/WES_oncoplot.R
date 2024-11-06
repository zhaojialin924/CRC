###WES######
###oncoplot####
#####filter Gene####
library(maftools)
library(ggplot2)
library(tidyverse)
library(ggsci)
library(data.table)
library(readxl)
library(sigminer)
library(NMF)
library(ggpubr)
library(tidyverse)
library(EnvStats)
library(ComplexHeatmap)
load("tmb_finalmaf.Rdata")
genes = finalmaf[, c("Hugo_Symbol"), drop = FALSE]
artifacts = c("LOC","ENS","FAM","GOL","PRA","NBP","POT","DEF","MUC","KRT","WAS","ANK","TRI","FRG",
              paste0("OR",1:9))
genes = genes[-which(substring(genes$Hugo_Symbol, 1, 3) %in% artifacts), , drop = FALSE]
artifacts2 = c("PLIN","CELA","SRA1")
genes = genes[-which(substring(genes$Hugo_Symbol, 1, 4) %in% artifacts2), , drop = FALSE]
artifacts_genes = c("ATXN1","PBRM1","ZNF814","MSH3","TTN","USH2A")

genes = genes[-which(genes$Hugo_Symbol %in% artifacts_genes),]
length(unique(genes$Hugo_Symbol))

finalmaf= finalmaf[which(finalmaf$Hugo_Symbol %in% genes$Hugo_Symbol ),]
dim(finalmaf) 

###for matrix###
finalmaf$Tumor_Sample_Barcode = as.character(finalmaf$Tumor_Sample_Barcode)
finalmaf$Hugo_Symbol = as.character(finalmaf$Hugo_Symbol)
MAF_gene_Count= finalmaf %>% 
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>% summarise(count=n()) %>% 
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%     
  summarise(count=n()) %>% group_by(Hugo_Symbol) %>% 
  summarise(count=n()) %>% arrange(desc(count))

matrix = matrix(nrow = length(unique(finalmaf$Tumor_Sample_Barcode)), ncol = nrow(MAF_gene_Count))
colnames(matrix) = MAF_gene_Count$Hugo_Symbol
rownames(matrix) = unique(finalmaf$Tumor_Sample_Barcode)

matrix[which(is.na(matrix))] = ""
matrix = matrix[,c(1:100)]

i = 1
for(i in 1:ncol(matrix)){
  gene = colnames(matrix)[i]
  sub = finalmaf[which(finalmaf$Hugo_Symbol == gene),]
  samples = sub$Tumor_Sample_Barcode
  matrix[,gene][which(rownames(matrix) %in% samples)] = 1
  matrix[,gene][which(!(rownames(matrix) %in% samples))] = 0
}

mat=as.data.frame(matrix)
mat$Patient_ID=rownames(mat)
save(mat,file = "WES_mutation_gene_no_artifacts.Rdata")
# load("WES_mutation_gene_no_artifacts.Rdata")

save(finalmaf,file = "oncoplot_finalmaf.Rdata")
load("oncoplot_finalmaf.Rdata")
maf <- read.maf(maf=finalmaf)
plotmafSummary(maf = maf, 
               rmOutlier = TRUE, 
               addStat = 'median', 
               dashboard = TRUE, 
               titvRaw = FALSE)
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)

options(repr.plot.width=16, repr.plot.height=5)
oncoplot(maf = maf,colors = vc_cols,
         top = 50,
         fontSize=0.8,
         sortByAnnotation = T,writeMatrix = T)
mat=read.table("onco_matrix.txt",sep="\t")
mat = mat %>% select(-c("P39"))
col = c("Missense_Mutation" = "#1F78B4", "Splice_Site" = "#FDBF6F", 
        "Nonsense_Mutation" = "#B2DF8A","Frame_Shift_Del" = "#A6CEE3",
        "Frame_Shift_Ins" = "#FB9A99","Multi_Hit" = "#33A02C",
        "In_Frame_Del"="#ff7f00","In_Frame_Ins"="#E31A1C")

alter_fun = list(
  background = alter_graphic("rect", fill = "white"),  ##e9e9e9 
  Missense_Mutation = alter_graphic("rect", fill = col["Missense_Mutation"]),
  Splice_Site = alter_graphic("rect", fill = col["Splice_Site"]),
  Nonsense_Mutation = alter_graphic("rect",fill = col["Nonsense_Mutation"],),
  Frame_Shift_Del = alter_graphic("rect", fill = col["Frame_Shift_Del"]),
  Frame_Shift_Ins = alter_graphic("rect", fill = col["Frame_Shift_Ins"]),
  Multi_Hit = alter_graphic("rect", fill = col["Multi_Hit"]),
  In_Frame_Del = alter_graphic("rect", fill = col["In_Frame_Del"]),
  In_Frame_Ins = alter_graphic("rect", fill = col["In_Frame_Ins"])
)

table(finalmaf$Variant_Classification)
###caculate TMB###
frequency_df = data.frame(table(finalmaf$Tumor_Sample_Barcode))
colnames(frequency_df) = c("sampleID", "Non_silent_Mutation_frequency")
frequency_df$Nonsilent_mutational_burden_per_Mb = frequency_df$Non_silent_Mutation_frequency / 40

clin=readxl::read_xlsx("01data/ours/clinical_infomation_20240413.xlsx")
pdata=clin %>% select(Patient_ID,Sex,Age,Tumor_Location,Stage1)
colnames(pdata)=c("sampleID","Sex","Age","Tumor_Location","Stage")
activity=read.csv("02res/ours/WES_signature.csv",row.names = 1)
activity$sampleID=rownames(activity)
pdata=merge(pdata,activity)
pdata=pdata[which(pdata$sampleID %in% finalmaf$Tumor_Sample_Barcode),]
pdata=merge(pdata,frequency_df)
pdata=pdata[which(pdata$sampleID %in% colnames(mat)),]
pdata = pdata[order(pdata$Non_silent_Mutation_frequency),]
sample_order = pdata$sampleID
pdata$Non_silent_Mutation_frequency1=log(pdata$Non_silent_Mutation_frequency)
mat = mat[,sample_order]
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
ha = HeatmapAnnotation(`log(Mutational load)` = anno_barplot(pdata$Non_silent_Mutation_frequency1, baseline = 0),
                       Sex = pdata$Sex,
                       `Tumor Location`=pdata$Tumor_Location,
                       Stage=pdata$Stage,
                       `SBS5`=pdata$Sig1,
                       `SBS10a`=pdata$Sig2,
                       `SBS6`=pdata$Sig3,
                       col = list(Stage = c("III" = "#EF767B", "II" = "#E9e7a7", "I" = "#43A3EF"),
                                  `Tumor Location` = c("Left hemicolon" = "#fd763f", 
                                                       "Rectum" = "#eeca40", 
                                                       "Right hemicolon" = "#23BAC5"),
                                  Sex = c("Female" = "pink", "Male" = "skyblue"),
                                  `SBS5`=col_fun,
                                  `SBS10a`=col_fun,
                                  `SBS6`=col_fun),
                       annotation_legend_param = list(
                         Stage=list(title_gp = gpar(col = "black", fontsize = 7),
                                    labels_gp = gpar(col = "black", fontsize = 6),
                                    direction = "horizontal",nrow=1),
                         `Tumor Location`=list(direction = "horizontal",nrow=1,
                                               title_gp = gpar(col = "black", fontsize = 7),
                                               labels_gp = gpar(col = "black", fontsize = 6)),
                         Sex=list(direction = "horizontal",nrow=1,
                                  title_gp = gpar(col = "black", fontsize = 7),
                                  labels_gp = gpar(col = "black", fontsize = 6)),
                         `SBS5`=list(direction = "horizontal",
                                     title_gp = gpar(col = "black", fontsize = 7),
                                     labels_gp = gpar(col = "black", fontsize = 6)),
                         `SBS10a`=list(direction = "horizontal",
                                       title_gp = gpar(col = "black", fontsize = 7),
                                       labels_gp = gpar(col = "black", fontsize = 6)),
                         `SBS6`=list(direction = "horizontal",
                                     title_gp = gpar(col = "black", fontsize = 7),
                                     labels_gp = gpar(col = "black", fontsize = 6))),
                       simple_anno_size = unit(.15,"cm"),
                       height = unit(2.1,"cm")
)

heatmap_legend_param = list(title = "Alternations", 
                            at = c('Missense_Mutation', 
                                   'Multi_Hit',                       
                                   'Nonsense_Mutation',
                                   'Frame_Shift_Del', 
                                   'Frame_Shift_Ins', 'Splice_Site',
                                   'In_Frame_Del','In_Frame_Ins'), 
                            labels = c('Missense_Mutation', 'Multi_Hit', 
                                       'Nonsense_Mutation',
                                       'Frame_Shift_Del', 
                                       'Frame_Shift_Ins', 'Splice_Site','In_Frame_Del','In_Frame_Ins'),
                            direction = "horizontal",nrow=2,
                            title_gp = gpar(col = "black", fontsize = 7),
                            labels_gp = gpar(col = "black", fontsize = 6))

pdf('oncoPrint4_top50.pdf',width = 7,height = 8)
ht=oncoPrint(mat,
             top_annotation = ha,
             alter_fun  = alter_fun, 
             alter_fun_is_vectorized = FALSE,
             pct_gp = gpar(fontsize = 8), 
             col = col,
             heatmap_legend_param = heatmap_legend_param,
             row_order = rownames(mat),
             column_order = sample_order) 
draw(ht, merge_legend = TRUE, heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")
dev.off()