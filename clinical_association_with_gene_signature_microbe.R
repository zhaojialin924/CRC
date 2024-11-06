library(gmodels)
library(sjPlot)
library(dplyr)
library(sjstats)
library(epitools)
library(epiR)
library(ggstatsplot)
library(ComplexHeatmap)
library(ggstatsplot)
setwd("D:/Postgraduation/microbe_crc")
####clinical phenotype-Mutation association####
###clinical###
clin1=readxl::read_xlsx("01data/ours/clinical_infomation_20240413.xlsx")
clin1$stage1=ifelse(grepl("III",clin1$Stage),1,0)
clin1$`Her-2`=ifelse(clin1$`Her-2`==0,0,1)
clin1$Positive_Lymph_Nodes1=ifelse(clin1$Positive_Lymph_Nodes>0,1,0)


load("WES_mutation_gene_no_artifacts.Rdata")
mat=mat[,c(1:50)]
mat1 <- data.frame(sapply(mat, function(x) as.factor(as.character(x))))
mat1$sample=rownames(mat)
data1=merge(clin1,mat1,by="sample")
x1=colnames(data1)[c(3,4,6,15,19,21,27,28)] ###clinical phenotype
x2=colnames(mat) ###gene
data1$`Lymphovascular Invasion(LVI)`=factor(data1$`Lymphovascular Invasion(LVI)`,
                                            levels = c(1,0))
data1$`Perineural Invasion`=factor(data1$`Perineural Invasion`,
                                   levels = c(1,0))
data1$Sex=factor(data1$Sex)
data1$`Lymph_Nodes>12`=factor(data1$`Lymph_Nodes>12`,levels = c(1,0))
data1$MMR=factor(data1$MMR,levels = c("dMMR","pMMR"))
data1$BRAF=factor(data1$BRAF,levels = c(1,0))
data1$stage1=factor(data1$stage1,levels = c(1,0))
data1$Positive_Lymph_Nodes1=factor(data1$Positive_Lymph_Nodes1,levels = c(1,0))

res <- data.frame(clin=character(), gene=character(), Cramers_v=numeric(), 
                  pvalue=numeric(),chi_squard=numeric(),OR=numeric(),PR=numeric(),
                  note=character(),explose=character(),
                  stringsAsFactors=FALSE)

results <- lapply(x1, function(x) {
  lapply(x2, function(y) {
    xx <- data1[, c(x, y)]
    colnames(xx) <- c("clin", "gene")
    xx$gene=factor(xx$gene,levels = c(1,0))
    corr1 <- crosstable_statistics(xx, x1 = clin, x2 = gene, statistics = "cramer")
    Cramers_v1=corr1$estimate
    pvalue1=corr1$p.value
    chi_squard1=as.numeric(corr1$statistic)
    explose1=as.character(levels(xx$clin)[1])
    if(all(table(xx$clin,xx$gene)!=0))
    {
      try(re <- epi.2by2(table(xx$clin,xx$gene)),silent = T)
      if(!inherits(re, "try-error")){
        pr=re$massoc.summary[1,2]
        or=re$massoc.summary[2,2]
        note1="unadjust"}
      else{
        pr=NA
        or=NA
        note1="error"
      }}
    else{
      try(re <- epi.2by2(table(xx$clin,xx$gene)+0.5),silent = T)
      if(!inherits(re,"try-error")){
        pr=re$massoc.summary[1,2]
        or=re$massoc.summary[2,2]
        note1="adjust"
      }
      else{
        pr=NA
        or=NA
        note1="error"
      }}
    c(clin=x, gene=y, Cramers_v=Cramers_v1,pvalue=pvalue1,chi_squard=chi_squard1,OR=or,
      PR=pr,note=note1,explose=explose1)
  })
})


res1 <- do.call(rbind, unlist(results, recursive=FALSE))
res=rbind(res,res1)
res[which(res$clin=="stage1"),]$clin="Stage"
res[which(res$clin=="Positive_Lymph_Nodes1"),]$clin="Positive_Lymph_Nodes"
res[which(res$clin=="MMR"),]$clin="MMR_Status"
res$pvalue=as.numeric(res$pvalue)
res$Cramers_v=as.numeric(res$Cramers_v)
res$chi_squard=as.numeric(res$chi_squard)
res$OR=as.numeric(res$OR)
res$PR=as.numeric(res$PR)

####clinical phenotype_Location#####
da=data1[which(data1$Tumor_Location!="Rectum"),]
da$Tumor_Location=factor(da$Tumor_Location,levels = c("Left hemicolon","Right hemicolon"))

results1 <- lapply(x2, function(y) {
  xx <- da[, c("Tumor_Location", y)]
  colnames(xx) <- c("clin", "gene")
  xx$gene=factor(xx$gene,levels = c(1,0))
  corr1 <- crosstable_statistics(xx, x1 = clin, x2 = gene, statistics = "cramer")
  Cramers_v1=corr1$estimate
  pvalue1=corr1$p.value
  chi_squard1=as.numeric(corr1$statistic)
  explose1=as.character(levels(xx$clin)[1])
  if(all(table(xx$clin,xx$gene)!=0))
  {
    try(re <- epi.2by2(table(xx$clin,xx$gene)),TRUE)
    pr=re$massoc.summary[1,2]
    or=re$massoc.summary[2,2]
    note1=""
  }else{
    try(re <- epi.2by2(table(xx$clin,xx$gene)+0.5),TRUE)
    pr=re$massoc.summary[1,2]
    or=re$massoc.summary[2,2]
    note1="adjust"
  }
  c(clin="Tumor_Location", gene=y, Cramers_v=Cramers_v1,pvalue=pvalue1,chi_squard=chi_squard1,OR=or,
    PR=pr,note=note1,explose=explose1)
})

results1=list(results1)
res1 <- do.call(rbind, unlist(results1, recursive=FALSE))
res=rbind(res,res1)
res$pvalue=as.numeric(res$pvalue)
res$Cramers_v=as.numeric(res$Cramers_v)
res$chi_squard=as.numeric(res$chi_squard)
res$OR=as.numeric(res$OR)
res$PR=as.numeric(res$PR)

res_sig=res %>% filter(pvalue<0.1)
write.table(res,"clinical_info_loca_gene_cramer_or.txt",row.names = F,quote = F,
            sep = "\t")
write.table(res_sig,"clinical_info_loca_gene_cramer_sig_or_0.1.txt",row.names = F,quote = F,
            sep = "\t")


clinical_info=unique(res_sig$clin)
gene=sort(unique(res_sig$gene))
or=matrix(NA,nrow=length(clinical_info),ncol = length(gene),
          dimnames = list(clinical_info,gene))

# for matrix
for(i in 1:nrow(res)) {
  clin_idx <- match(res$clin[i], clinical_info)
  gene_idx <- match(res$gene[i], gene)
  or[clin_idx,gene_idx] <- res$OR[i]
}

p=matrix(NA,nrow=length(clinical_info),ncol = length(gene),
         dimnames = list(clinical_info,gene))
# for p matrix
for(i in 1:nrow(res)) {
  clin_idx <- match(res$clin[i], clinical_info)
  gene_idx <- match(res$gene[i], gene)
  p[clin_idx,gene_idx] <- res$pvalue[i]
}
p=as.matrix(p)

p[p<0.01]="**"
p[p>0.01 & p<0.05]="*"

p[p>0.05 &p<0.1]="·"
p[p>0.1]=""


ht <- Heatmap(or,
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
                title = "** pvalue < 0.01\n*  0.01 < pvalue < 0.05\n·  0.05 < pvalue < 0.1\n\nOR"))

draw(ht)
pdf("heatmap_clin_gene_OR.pdf",width = 13,height = 5)
draw(ht)
dev.off()

####clinical phenotype-siganture association####
activity=read.csv("02res/ours/WES_signature.csv",row.names = 1)
activity$sample=rownames(activity)
data2=merge(clin1,activity)
data2$group=max.col(data2[,c(29:31)])
data2$Sig1_SBS5=0
data2$Sig2_SBS10a=0
data2$Sig3_SBS6=0
data2[which(data2$group==1),]$Sig1_SBS5=1
data2[which(data2$group==2),]$Sig2_SBS10a=1
data2[which(data2$group==3),]$Sig3_SBS6=1

res <- data.frame(clin=character(), signature=character(), Cramers_v=numeric(), 
                  pvalue=numeric(),chi_squard=numeric(),OR=numeric(),PR=numeric(),
                  note=character(),explose=character(),
                  stringsAsFactors=FALSE)
x2=c("Sig1_SBS5","Sig2_SBS10a","Sig3_SBS6") ##signature
results <- lapply(x1, function(x) {
  lapply(x2, function(y) {
    xx <- data2[, c(x, y)]
    colnames(xx) <- c("clin", "sig")
    xx$sig=factor(xx$sig,levels = c(1,0))
    corr1 <- crosstable_statistics(xx, x1 = clin, x2 = sig, statistics = "cramer")
    Cramers_v1=corr1$estimate
    pvalue1=corr1$p.value
    chi_squard1=as.numeric(corr1$statistic)
    explose1=as.character(levels(xx$clin)[1])
    if(all(table(xx$clin,xx$gene)!=0))
    {
      try(re <- epi.2by2(table(xx$clin,xx$gene)),TRUE)
      pr=re$massoc.summary[1,2]
      or=re$massoc.summary[2,2]
      note1=""
    }else{
      try(re <- epi.2by2(table(xx$clin,xx$gene)+0.5),TRUE)
      pr=re$massoc.summary[1,2]
      or=re$massoc.summary[2,2]
      note1="adjust"
    }
    c(clin=x, signature=y, Cramers_v=Cramers_v1,pvalue=pvalue1,chi_squard=chi_squard1,
      OR=or,PR=pr,explose=levels(xx$clin)[1],note=note1)
  })
})


res1 <- do.call(rbind, unlist(results, recursive=FALSE))
res=rbind(res,res1)
res[which(res$clin=="stage1"),]$clin="Stage"
res[which(res$clin=="Positive_Lymph_Nodes1"),]$clin="Positive_Lymph_Nodes"
res[which(res$clin=="MMR"),]$clin="MMR_Status"
res$OR=as.numeric(res$OR)
res$PR=as.numeric(res$PR)
res$Cramers_v=as.numeric(res$Cramers_v)
res$pvalue=as.numeric(res$pvalue)
res$chi_squard=as.numeric(res$chi_squard)
####Location#####
da=data2[which(data2$Tumor_Location!="Rectum"),]
da$Tumor_Location=factor(da$Tumor_Location,levels = c("Left hemicolon","Right hemicolon"))
x3="Tumor_Location"
results <- lapply(x3, function(x) {
  lapply(x2, function(y) {
    xx <- da[, c(x, y)]
    colnames(xx) <- c("clin", "sig")
    xx$sig=factor(xx$sig,levels = c(1,0))
    corr1 <- crosstable_statistics(xx, x1 = clin, x2 = sig, statistics = "cramer")
    Cramers_v1=corr1$estimate
    pvalue1=corr1$p.value
    chi_squard1=as.numeric(corr1$statistic)
    explose1=as.character(levels(xx$clin)[1])
    if(all(table(xx$clin,xx$sig)!=0))
    {
      try(re <- epi.2by2(table(xx$clin,xx$sig)),TRUE)
      pr=re$massoc.summary[1,2]
      or=re$massoc.summary[2,2]
      note1=""
    }else{
      try(re <- epi.2by2(table(xx$clin,xx$sig)+0.5),TRUE)
      pr=re$massoc.summary[1,2]
      or=re$massoc.summary[2,2]
      note1="adjust"
    }
    c(clin=x, signature=y, Cramers_v=Cramers_v1,pvalue=pvalue1,chi_squard=chi_squard1,
      OR=or,PR=pr,explose=levels(xx$clin)[1],note=note1)
  })
})


res1 <- do.call(rbind, unlist(results, recursive=FALSE))
res=rbind(res,res1)
res$pvalue=as.numeric(res$pvalue)
res$Cramers_v=as.numeric(res$Cramers_v)
res$chi_squard=as.numeric(res$chi_squard)
res$OR=as.numeric(res$OR)
res$PR=as.numeric(res$PR)

res_sig=res %>% filter(pvalue<0.1)
write.table(res,"clinical_info_loca_signature_cramer_or.txt",row.names = F,quote = F,
            sep = "\t")
write.table(res_sig,"clinical_info_loca_signature_cramer_sig_or_0.1.txt",row.names = F,quote = F,
            sep = "\t")


clinical_info=unique(res_sig$clin)
sig=sort(unique(res_sig$signature))
or=matrix(NA,nrow=length(clinical_info),ncol = length(sig),
          dimnames = list(clinical_info,sig))

for(i in 1:nrow(res)) {
  clin_idx <- match(res$clin[i], clinical_info)
  sig_idx <- match(res$signature[i], sig)
  or[clin_idx,sig_idx] <- res$OR[i]
}

p=matrix(NA,nrow=length(clinical_info),ncol = length(sig),
         dimnames = list(clinical_info,sig))

for(i in 1:nrow(res)) {
  clin_idx <- match(res$clin[i], clinical_info)
  sig_idx <- match(res$signature[i], sig)
  p[clin_idx,sig_idx] <- res$pvalue[i]
}
p=as.matrix(p)
p[p<0.01]="**"
p[p>0.01 & p<0.05]="*"
p[p>0.05 &p<0.1]="·"
p[p>0.1]=""
library(ComplexHeatmap)

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
                title = "** pvalue < 0.01\n*  0.01 < pvalue < 0.05\n·  0.05 < pvalue < 0.1\n\nOR"))

draw(ht)
pdf("heatmap_clin_signature_OR.pdf",width = 13,height = 5)
draw(ht)
dev.off()

###bar plot####
data2$SBS6=ifelse(data2$SBS6==1,"SBS6","Others")
data2$SBS6=factor(data2$SBS6,levels = c("SBS6","Others"))
data2$MMR=factor(data2$MMR,levels = c("dMMR","pMMR"))
ggbarstats(
  data         = data2,
  x            = MMR,
  y            = SBS6,
  xlab="",
  legend.title = "MMR status", 
  ggtheme          = theme_classic(),
)+theme(axis.text = element_text(color="black",size = 14))
ggsave("MMR_SBS6.pdf",width = 3.2,height = 4)



####clinical phenotype-microbes association####
clin1=readxl::read_xlsx("01data/ours/clinical_infomation_20240413.xlsx")
clin1$stage1=ifelse(grepl("III",clin1$Stage),1,0)
clin1$Positive_Lymph_Nodes1=ifelse(clin1$Positive_Lymph_Nodes>0,1,0)
colnames(clin1)[c(3,4)]=colnames(xx)[c(4,5)]
clin1$`Her-2`=ifelse(clin1$`Her-2`==0,0,1)
data11=read.delim("count/species.filter.tsv",sep = "\t")
data11=as.data.frame(t(data11))
data11=data11/rowSums(data11)
data11=as.data.frame(t(data11))
data11$var=apply(data11,1,var)
data11=data11[order(-data11$var),]
top50=head(data11,50)
top50=top50 %>% select(-var)
top50=as.matrix(top50)
top50[which(top50==0)]="Absent"
top50[which(top50!="Absent")]="Present"
res_ba=as.data.frame(t(top50))
res_ba$sample=rownames(res_ba)
data2=merge(clin1,res_ba,by="sample")


x2=unique(rownames(top50)) ##microbe

data2$Vascular_Cancer_Thrombus=factor(data2$Cancer_Nodules,levels = c(1,0))
data2$`Lymphovascular Invasion(LVI)`=factor(data2$`Lymphovascular Invasion(LVI)`,
                                            levels = c(1,0))
data2$`Perineural Invasion`=factor(data2$`Perineural Invasion`,
                                   levels = c(1,0))
data2$Sex=factor(data2$Sex)
data2$`Lymph_Nodes>12`=factor(data2$`Lymph_Nodes>12`,levels = c(1,0))
data2$MMR=factor(data2$MMR,levels = c("dMMR","pMMR"))
data2$BRAF=factor(data2$BRAF,levels = c(1,0))
data2$Positive_Lymph_Nodes1=factor(data2$Positive_Lymph_Nodes1,levels = c(1,0))
data2$stage1=factor(data2$stage1,levels = c(1,0))
data2$`Her-2`=factor(data2$`Her-2`,levels = c(1,0))
res <- data.frame(clin=character(), Bacteria=character(), Cramers_v=numeric(), 
                  pvalue=numeric(),chi_squard=numeric(),OR=numeric(),PR=numeric(),
                  note=character(),explose=character(),
                  stringsAsFactors=FALSE)

results <- lapply(x1, function(x) {
  lapply(x2, function(y) {
    xx <- data2[, c(x, y)]
    colnames(xx) <- c("clin", "bacteria")
    xx$bacteria=factor(xx$bacteria,levels = c("Present","Absent"))
    corr1 <- crosstable_statistics(xx, x1 = clin, x2 = bacteria, statistics = "cramer")
    Cramers_v1=corr1$estimate
    pvalue1=corr1$p.value
    chi_squard1=as.numeric(corr1$statistic)
    explose1=as.character(levels(xx$clin)[1])
    if(all(table(xx$clin,xx$bacteria)!=0))
    {
      try(re <- epi.2by2(table(xx$clin,xx$bacteria)),TRUE)
      pr=re$massoc.summary[1,2]
      or=re$massoc.summary[2,2]
      note1=""
    }else{
      try(re <- epi.2by2(table(xx$clin,xx$bacteria)+0.5),TRUE)
      pr=re$massoc.summary[1,2]
      or=re$massoc.summary[2,2]
      note1="adjust"
    }
    c(clin=x, Bacteria=y, Cramers_v=Cramers_v1,pvalue=pvalue1,chi_squard=chi_squard1,
      OR=or,PR=pr,explose=levels(xx$clin)[1],note=note1)
  })
})


res1 <- do.call(rbind, unlist(results, recursive=FALSE))
res=rbind(res,res1)
res[which(res$clin=="stage1"),]$clin="Stage"
res[which(res$clin=="Positive_Lymph_Nodes1"),]$clin="Positive_Lymph_Nodes"
res[which(res$clin=="MMR"),]$clin="MMR_Status"
res$OR=as.numeric(res$OR)
res$PR=as.numeric(res$PR)
res$Cramers_v=as.numeric(res$Cramers_v)
res$pvalue=as.numeric(res$pvalue)
res$chi_squard=as.numeric(res$chi_squard)


####Location#####
da=data2[which(data2$Tumor_Location!="Rectum"),]
da$Tumor_Location=factor(da$Tumor_Location,levels = c("Left hemicolon","Right hemicolon"))
x3="Tumor_Location"
results <- lapply(x3, function(x) {
  lapply(x2, function(y) {
    xx <- da[, c(x, y)]
    colnames(xx) <- c("clin", "bacteria")
    xx$bacteria=factor(xx$bacteria,levels = c("Present","Absent"))
    corr1 <- crosstable_statistics(xx, x1 = clin, x2 = bacteria, statistics = "cramer")
    Cramers_v1=corr1$estimate
    pvalue1=corr1$p.value
    chi_squard1=as.numeric(corr1$statistic)
    explose1=as.character(levels(xx$clin)[1])
    if(all(table(xx$clin,xx$bacteria)!=0))
    {
      try(re <- epi.2by2(table(xx$clin,xx$bacteria)),TRUE)
      pr=re$massoc.summary[1,2]
      or=re$massoc.summary[2,2]
      note1=""
    }else{
      try(re <- epi.2by2(table(xx$clin,xx$bacteria)+0.5),TRUE)
      pr=re$massoc.summary[1,2]
      or=re$massoc.summary[2,2]
      note1="adjust"
    }
    c(clin=x, Bacteria=y, Cramers_v=Cramers_v1,pvalue=pvalue1,chi_squard=chi_squard1,
      OR=or,PR=pr,explose=levels(xx$clin)[1],note=note1)
  })
})


res1 <- do.call(rbind, unlist(results, recursive=FALSE))
res=rbind(res,res1)
res$pvalue=as.numeric(res$pvalue)
res$Cramers_v=as.numeric(res$Cramers_v)
res$chi_squard=as.numeric(res$chi_squard)
res$OR=as.numeric(res$OR)
res$PR=as.numeric(res$PR)

res_sig=res %>% filter(pvalue<0.1)
write.table(res,"clinical_info_loca_top50_bacteria_cramer_or.txt",row.names = F,quote = F,
            sep = "\t")
write.table(res_sig,"clinical_info_loca_top50__bacteria_cramer_sig_or_0.1.txt",row.names = F,quote = F,
            sep = "\t")


clinical_info=unique(res_sig$clin)
Bacteria=sort(unique(res_sig$Bacteria))
clinical_info=unique(res_sig$clin)
or=matrix(NA,nrow=length(clinical_info),ncol = length(Bacteria),
          dimnames = list(clinical_info,Bacteria))
# 填充矩阵
for(i in 1:nrow(res)) {
  clin_idx <- match(res$clin[i], clinical_info)
  ba_idx <- match(res$Bacteria[i], Bacteria)
  hm[clin_idx,ba_idx] <- res$OR[i]
}
clinical_info=unique(res_sig$clin)
p=matrix(NA,nrow=length(clinical_info),ncol = length(Bacteria),
         dimnames = list(clinical_info,Bacteria))
# 填充矩阵
for(i in 1:nrow(res)) {
  clin_idx <- match(res$clin[i], clinical_info)
  ba_idx <- match(res$Bacteria[i], Bacteria)
  p[clin_idx,ba_idx] <- res$pvalue[i]
}
p=as.matrix(p)
p[p<0.01]="**"
p[p>0.01 & p<0.05]="*"
p[p>0.05 &p<0.1]="·"
p[p>0.1]=""


ht <- Heatmap(or,
              cluster_rows = F,
              cluster_columns = F,
              name = "OR",
              row_title = NULL,
              col=circlize::colorRamp2(c(0,1,10),c("#619cff", "white", "red")),
              show_row_names = TRUE,
              show_column_names = TRUE,
              width = ncol(hm)*unit(6,"mm"),
              height = nrow(hm)*unit(6,"mm"),
              row_names_gp = gpar(fontsize=10),
              column_names_rot = 45,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(p[i,j], x = x, y = y,gp = gpar(fontsize = 10,col="black"))},
              heatmap_legend_param = list(
                title = "** pvalue < 0.01\n*  0.01 < pvalue < 0.05\n·  0.05 < pvalue < 0.1\n\nOR"))

draw(ht)
pdf("heatmap_clin_top50_bacteria_OR.pdf",width = 13,height = 5)
draw(ht)
dev.off()
