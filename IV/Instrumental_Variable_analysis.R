###IV###Instrumental Variable analysis
setwd("D:/Postgraduation/microbe_crc/")
library(ivreg)
library(sandwich)
library(data.table)
library(parallel)
library(haven)
library(lmtest)
library(readxl)
library(ggthemes)
library(RColorBrewer)

load("WES_mutation_gene_no_artifacts.Rdata")
mat=mat[,c(1:50)]
mat1=mat
mat1$sample=rownames(mat1)
mat1 <- data.frame(sapply(mat1, function(x) as.factor(as.character(x))))
###sginature###
activity=read.csv("02res/ours/WES_signature.csv",row.names = 1)
activity$sample=rownames(activity)


###caculating the gene and signature correlation####
xx=merge(activity,mat1)
results=data.frame()
for (i in 1:50) {
  gene=colnames(mat)[i]
  for (j in 1:3) {
    sig <- colnames(activity)[j]
    data=xx[,c(gene,sig)]
    colnames(data)=c("gene","sig")
    data$gene=as.numeric(data$gene)
    # Perform linear regression
    fit <- lm(sig ~ gene, data = data)
    # Extract coefficient and p-value
    coeff <- coef(summary(fit))[2, 1]
    pvalue <- coef(summary(fit))[2, 4]
    # Store coefficient and p-value along with gene name
    gene_result <- data.frame(coeff = coeff, pvalue = pvalue, gene = gene,sig=sig)
    results=rbind(results,gene_result)
  }
}

###Remove genes related to signature####
sig1_gene=unique(results[which(results$pvalue>0.05 & results$sig=="Sig1"),]$gene)
sig2_gene=unique(results[which(results$pvalue>0.05 & results$sig=="Sig2"),]$gene)
sig3_gene=unique(results[which(results$pvalue>0.05 & results$sig=="Sig3"),]$gene)

clin=readxl::read_xlsx("01data/ours/clinical_infomation_20240413.xlsx")
activity=read.csv("02res/ours/WES_signature.csv",row.names = 1)
activity$sample=rownames(activity)
clin1=merge(clin1,activity)
###new###
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


###signature 1 SBS5######
load("WES_mutation_gene_no_artifacts.Rdata")
mat=mat[,c(1:50)]
mat=mat[,sig1_gene]
mat1 <- data.frame(sapply(mat, function(x) as.factor(as.character(x))))
mat1$sample=rownames(mat)
data1=merge(clin1,mat1,by="sample")
data1=merge(data1,top50hvb,by="sample")
### sig1~ microbe abundance+Sex+age+location+stage | gene mutational status

res=data.frame()
for(i in 1:length(sig1_gene)){
  for (j in 1:50){
    re2=ivreg(formula = data1[,28]~data1[,43+j]+data1[, 6] + data1[, 10] + data1[, 24]+data1[,27]
              | data1[,30+i]+data1[, 6] + data1[, 10] + data1[, 24]+data1[,27] ,x=T,data=data1)
    try(su <- summary(re2, df = Inf, diagnostics = TRUE),TRUE)
    try(re_p <- data.frame(Sig = colnames(data1[28]), 
                           Gene= colnames(data1[(30+i)]),
                           Bacteria=colnames(data1)[(43+j)],
                           coef= su$coefficients[2,1], 
                           wu_Hausman = su$diagnostics[2,4],
                           iv_p.value = su$diagnostics[1,4], 
                           p.value= su$coefficients[2,4],
                           adjusted_r_squard = su$adj.r.squared),TRUE)
    res=rbind(res,re_p)
  }
}

###sig2###
load("WES_mutation_gene_no_artifacts.Rdata")
mat=mat[,c(1:50)]
mat=mat[,sig2_gene]
# 将数据框df中的所有列转换为数值型
mat1 <- data.frame(sapply(mat, function(x) as.factor(as.character(x))))
mat1$sample=rownames(mat)
data1=merge(clin1,mat1,by="sample")
data1=merge(data1,top50hvb,by="sample")
res2=data.frame()
for(i in 1:length(sig2_gene)){
  for (j in 1:50){
    re2=ivreg(formula = data1[,29]~data1[,55+j] +data1[, 6] + data1[, 10] + data1[, 24]+data1[,27]  
              | data1[,30+i]+data1[, 6] + data1[, 10] + data1[, 24] +data1[,27],x=T,data=data1)
    try(su <- summary(re2, df = Inf, diagnostics = TRUE),TRUE)
    try(re_p <- data.frame(Sig = colnames(data1[29]), 
                           Gene= colnames(data1[(30+i)]),
                           Bacteria=colnames(data1)[(55+j)],
                           coef= su$coefficients[2,1], 
                           wu_Hausman = su$diagnostics[2,4],
                           iv_p.value = su$diagnostics[1,4], 
                           p.value= su$coefficients[2,4],
                           adjusted_r_squard = su$adj.r.squared),TRUE)
    res2=rbind(res2,re_p)
  }
}

###sig3 SBS6###
load("WES_mutation_gene_no_artifacts.Rdata")
mat=mat[,c(1:50)]
mat=mat[,sig3_gene]
mat1 <- data.frame(sapply(mat, function(x) as.factor(as.character(x))))
mat1$sample=rownames(mat)
data1=merge(clin1,mat1,by="sample")
data1=merge(data1,top50hvb,by="sample")
res3=data.frame()
for(i in 1:length(sig3_gene)){
  for (j in 1:50){
    re2=ivreg(formula = data1[,30]~ data1[,61+j]+data1[, 6]+data1[,24] + data1[, 10]+data1[, 27]
              | data1[,30+i]+data1[, 6]+data1[,24] + data1[, 10]+data1[, 27] ,x=T,y=T,data=data1)
    try(su <- summary(re2, df = Inf, diagnostics = TRUE),TRUE)
    try(re_p <- data.frame(Sig = colnames(data1[30]), 
                           Gene= colnames(data1[(30+i)]),
                           Bacteria=colnames(data1)[(61+j)],
                           coef= su$coefficients[2,1], 
                           wu_Hausman = su$diagnostics[2,4],
                           iv_p.value = su$diagnostics[1,4], 
                           p.value= su$coefficients[2,4],
                           adjusted_r_squard = su$adj.r.squared),TRUE)
    res3=rbind(res3,re_p)
  }
}


res_all1=rbind(res,res2,res3)
res_all1 <- na.omit(res_all1)
res_sig_all1 <- subset(res_all1,wu_Hausman<0.05 & iv_p.value < 0.05 & p.value<0.1)

write.csv(res_all1,"02res/ours/IV_res_all.csv",row.names = F,quote = F)
write.csv(res_sig_all1,"02res/ours/IV_res_sig.csv",row.names = F,quote = F)