#####validation Sidra cohort#####
library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(sjstats)
library(epiR)
setwd("D:/Postgraduation/microbe_crc")

M_gene=read.delim("clinical_info_loca_gene_cramer_sig_or_0.1_0920.txt")
M_gene=unique(M_gene$gene)
M_gene=M_gene[-c(5,26,28)]
load("sidra_mutation_gene_status.Rdata")
matrix=as.data.frame(matrix)
sidra_mgene=matrix[,M_gene]
sidra_mgene$Patient_ID=substr(rownames(sidra_mgene),1,17)

os=readxl::read_xlsx("01data/ACICAM/41591_2023_2324_MOESM4_ESM.xlsx",
                     sheet = "2.1_Clinical_Data_AC_ICAM",
                     skip = 1)
os_sidra_mgene=merge(os,sidra_mgene)
mmr=read_xlsx("01data/ACICAM/41591_2023_2324_MOESM4_ESM.xlsx",
              sheet = "11_MANTIS",
              skip = 1)
mmr=mmr %>% select(Patient_ID,MSI)
os_sidra_mgene=merge(os_sidra_mgene,mmr)

os_sidra_mgene$Tumor_Location=ifelse(os_sidra_mgene$Tumor_anatomic_location %in% 
                                       c("ceceum","colon ascendens","flexura hepatica",
                                         "colon transversum"),"Right","Left")
os_sidra_mgene$Stage=ifelse(os_sidra_mgene$AJCC_path_stage %in% c(1,2),0,1)
colnames(os_sidra_mgene)[2]="Sex"
os_sidra_mgene$MSI=factor(os_sidra_mgene$MSI,levels = c("MSI-H","MSS"))
os_sidra_mgene$Sex=factor(os_sidra_mgene$Sex,levels = c("FEMALE","MALE"))
os_sidra_mgene$Tumor_Location=factor(os_sidra_mgene$Tumor_Location,levels = c("Left","Right"))
os_sidra_mgene$Stage=factor(os_sidra_mgene$Stage,levels = c(1,0))
x1=c("Stage","Sex","MSI","Tumor_Location")

res <- data.frame(clin=character(), gene=character(), Cramers_v=numeric(), 
                  pvalue=numeric(),chi_squard=numeric(),OR=numeric(),
                  PR=numeric(),explose=character(),note=character(),
                  stringsAsFactors=FALSE)

results <- lapply(x1, function(x) {
  lapply(M_gene, function(y) {
    xx <- os_sidra_mgene[, c(x, y)]
    colnames(xx) <- c("clin", "gene")
    corr1 <- crosstable_statistics(xx, x1 = clin, x2 = gene, statistics = "cramer")
    Cramers_v1=corr1$estimate
    pvalue1=corr1$p.value
    chi_squard1=corr1$statistic
    xx$gene=factor(xx$gene,levels = c("MUT","WT"))
    if(all(table(xx$clin,xx$gene)!=0))
    {
      try(re <- epi.2by2(table(xx$clin,xx$gene)),TRUE)
      pr=re$massoc.summary[1,2]
      or=re$massoc.summary[2,2]
      note1=""
    }else{
      or=as.numeric(fisher.test(table(xx$clin,xx$gene))$estimate["odds ratio"])
      pr=""
      note1="Fisher"
    }
    c(clin=x, gene=y, Cramers_v=Cramers_v1,pvalue=pvalue1,chi_squard=chi_squard1,OR=or,
      PR=pr,explose=levels(xx$clin)[1],note=note1)
  })
})


res1 <- do.call(rbind, unlist(results, recursive=FALSE))
res=rbind(res,res1)
res$Cramers_v=as.numeric(res$Cramers_v)
res$pvalue=as.numeric(res$pvalue)
res$`chi_squard.Chi-squared`=as.numeric(res$`chi_squard.Chi-squared`)
res$OR=as.numeric(res$OR)
res$PR=as.numeric(res$PR)
res_fdr <- lapply(unique(x1), function(i){
  df <- subset( res,clin == i)
  df$fdr <- p.adjust(df$pvalue, method = "BH")
  return(df)
})
res_fdr <- do.call(rbind, res_fdr)
res_fdr$fdr=as.numeric(res_fdr$fdr)
res_sig=res_fdr %>% filter(fdr<0.1)
write.table(res_fdr,"sidra_clinical_gene_res_fdr.txt",row.names = F,quote = F,sep = "\t")
write.table(res_sig,"Sidra_clinical_gene_res_0.1_sig_fdr.txt",row.names = F,quote = F,
            sep = "\t")


clinical_info=unique(res_sig$clin)
gene=sort(unique(res_sig$gene))

hm=matrix(NA,nrow=length(clinical_info),ncol = length(gene),
          dimnames = list(clinical_info,gene))

for(i in 1:nrow(res_fdr)) {
  clin_idx <- match(res_fdr$clin[i], clinical_info)
  gene_idx <- match(res_fdr$gene[i], gene)
  hm[clin_idx,gene_idx] <- res_fdr$OR[i]
}

p=matrix(NA,nrow=length(clinical_info),ncol = length(gene),

for(i in 1:nrow(res_fdr)) {
  clin_idx <- match(res_fdr$clin[i], clinical_info)
  gene_idx <- match(res_fdr$gene[i], gene)
  p[clin_idx,gene_idx] <- res_fdr$fdr[i]
}
p=as.matrix(p)

p[p<0.01]="**"
p[p>0.01 & p<0.05]="*"

p[p>0.05 &p<0.1]="·"
p[p>0.1]=""

ht <- Heatmap(hm,
              cluster_rows = F,
              cluster_columns = F,
              name = "OR",
              row_title = NULL,
              col=circlize::colorRamp2(c(0,1,20),c("#619cff", "white", "red")),
              show_row_names = TRUE,
              show_column_names = TRUE,
              width = ncol(hm)*unit(6,"mm"),
              height = nrow(hm)*unit(6,"mm"),
              row_names_gp = gpar(fontsize=10),
              column_names_rot = 45,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(p[i,j], x = x, y = y,gp = gpar(fontsize = 10,col="black"))},
              heatmap_legend_param = list(
                title = "** FDR < 0.01\n*  0.01 < FDR < 0.05\n·  0.05 < FDR < 0.1\n\nOR"))

draw(ht)
pdf("heatmap_sidra_clin_gene_fdr_OR.pdf",width = 13,height = 5)
draw(ht)
dev.off()

IV_res=read.csv("ivres_ba_cor_spearman_all_result.csv")
IV_res=IV_res[which(IV_res$fdr<0.05),]

ba=unique(IV_res$ba2)
anno=read.delim("count/annotation.tsv")
anno1=anno[which(anno$species %in% ba),]
micro=readxl::read_xlsx("01data/ACICAM/41591_2023_2324_MOESM4_ESM.xlsx",
                        sheet = "13_16S_Relative_Abundance",
                        skip = 1)
micro=as.data.frame(micro)
rownames(micro)=micro$Taxa
micro=micro[,grepl("AN",colnames(micro))]
micro$Taxa=c(unlist(lapply(rownames(micro), function(x) regmatches(x,regexpr("(D_5__.*)$", x, perl = TRUE)))),"Unknow")
head(micro$Taxa)
micro$Taxa=gsub("D_5__","",micro$Taxa)

micro1=micro[which(micro$Taxa %in% unique(anno1$genus)),]
rownames(micro1)=micro1$Taxa
micro1=micro1[,-247]
micro1=as.data.frame(t(micro1))
micro1$Patient_ID=substr(rownames(micro1),1,17)

for(mi in colnames(micro1)[c(1:11)]){
  micro1[[mi]]=ifelse(micro1[[mi]]==0,"Absent","Present")
}

os_sidra_micro=merge(os,micro1)
mmr=read_xlsx("01data/ACICAM/41591_2023_2324_MOESM4_ESM.xlsx",
              sheet = "11_MANTIS",
              skip = 1)
mmr=mmr %>% select(Patient_ID,MSI)
os_sidra_micro=merge(os_sidra_micro,mmr)

os_sidra_micro$Tumor_Location=ifelse(os_sidra_micro$Tumor_anatomic_location %in% 
                                       c("ceceum","colon ascendens","flexura hepatica",
                                         "colon transversum"),"Right","Left")

colnames(os_sidra_micro)[2]="Sex"
os_sidra_micro$Stage=ifelse(os_sidra_micro$AJCC_path_stage %in% c(1,2),0,1)
os_sidra_micro$MSI=factor(os_sidra_micro$MSI,levels = c("MSI-H","MSS"))
os_sidra_micro$Sex=factor(os_sidra_micro$Sex,levels = c("FEMALE","MALE"))
os_sidra_micro$Tumor_Location=factor(os_sidra_micro$Tumor_Location,levels = c("Left","Right"))
os_sidra_micro$Stage=factor(os_sidra_micro$Stage,levels = c(1,0))
x1=c("Stage","Sex","MSI","Tumor_Location")


res <- data.frame(clin=character(), bacteria=character(), Cramers_v=numeric(), 
                  pvalue=numeric(),chi_squard=numeric(),OR=numeric(),
                  PR=numeric(),explose=character(),
                  stringsAsFactors=FALSE)

results <- lapply(x1, function(x) {
  lapply(colnames(micro1)[c(1:11)], function(y) {
    xx <- os_sidra_micro[, c(x, y)]
    colnames(xx) <- c("clin", "micro")
    corr1 <- crosstable_statistics(xx, x1 = clin, x2 = micro, statistics = "cramer")
    Cramers_v1=corr1$estimate
    
    pvalue1=corr1$p.value
    chi_squard1=corr1$statistic
    xx$micro=factor(xx$micro,levels = c("Present","Absent"))
    try(re <- epi.2by2(table(xx$clin,xx$micro)),TRUE)
    pr=re$massoc.summary[1,2]
    or=re$massoc.summary[2,2]
    c(clin=x, bacteria=y, Cramers_v=Cramers_v1,pvalue=pvalue1,chi_squard=chi_squard1,OR=or,
      PR=pr,explose=levels(xx$clin)[1])
  })
})


res1 <- do.call(rbind, unlist(results, recursive=FALSE))
res=rbind(res,res1)
res$Cramers_v=as.numeric(res$Cramers_v)
res$pvalue=as.numeric(res$pvalue)
res$`chi_squard.Chi-squared`=as.numeric(res$`chi_squard.Chi-squared`)
res$OR=as.numeric(res$OR)
res$PR=as.numeric(res$PR)
res_fdr <- lapply(unique(x1), function(i){
  df <- subset( res,clin == i)
  df$fdr <- p.adjust(df$pvalue, method = "BH")
  return(df)
})
res_fdr <- do.call(rbind, res_fdr)
res_fdr$fdr=as.numeric(res_fdr$fdr)

res_sig=res_fdr %>% filter(fdr<0.1)
write.table(res_fdr,"sidra_clinical_bacteria_res_fdr.txt",row.names = F,quote = F,sep = "\t")
write.table(res_sig,"Sidra_clinical_bacteria_res_0.1_sig_fdr.txt",row.names = F,quote = F,
            sep = "\t")


clinical_info=unique(res_sig$clin)
micro=unique(res_sig$bacteria)
hm=matrix(NA,nrow=length(clinical_info),ncol = length(micro),
          dimnames = list(clinical_info,micro))

for(i in 1:nrow(res_fdr)) {
  clin_idx <- match(res_fdr$clin[i], clinical_info)
  micro_idx <- match(res_fdr$bacteria[i], micro)
  hm[clin_idx,micro_idx] <- res_fdr$OR[i]
}
p=matrix(NA,nrow=length(clinical_info),ncol = length(micro),
         dimnames = list(clinical_info,micro))
for(i in 1:nrow(res_fdr)) {
  clin_idx <- match(res_fdr$clin[i], clinical_info)
  micro_idx <- match(res_fdr$bacteria[i], micro)
  p[clin_idx,micro_idx] <- res_fdr$pvalue[i]
}
p=as.matrix(p)
p[p<0.01]="**"
p[p>0.01 & p<0.05]="*"
p[p>0.05]=""
ht <- Heatmap(hm,
              cluster_rows = F,
              cluster_columns = F,
              name = "OR",
              row_title = NULL,
              col=circlize::colorRamp2(c(0,1,3),c("#619cff", "white", "red")),
              show_row_names = TRUE,
              show_column_names = TRUE,
              width = ncol(hm)*unit(6,"mm"),
              height = nrow(hm)*unit(6,"mm"),
              row_names_gp = gpar(fontsize=10),
              column_names_rot = 45,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(p[i,j], x = x, y = y,gp = gpar(fontsize = 10,col="black"))},
              heatmap_legend_param = list(
                title = "** pvalue < 0.01\n*  0.01 < pvalue < 0.05\n\n\nOR"))

draw(ht)
pdf("heatmap_sidra_clin_pvalue_bacteria_OR.pdf",width = 5,height = 5)
draw(ht)
dev.off()

