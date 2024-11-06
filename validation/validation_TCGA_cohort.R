library(TCGAbiolinks)
library(maftools)
library(ComplexHeatmap)
library(dplyr)
query <- GDCquery(
  project = "TCGA-COAD", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)

GDCdownload(query)

GDCprepare(query, save = T,save.filename = "TCGA-COAD_SNP.Rdata")
load(file = "TCGA-COAD_SNP.Rdata")
maf.coad <- data
maf <- read.maf(maf.coad)
M_gene=read.delim("clinical_info_loca_gene_cramer_sig_or_0.1.txt")
M_gene=unique(M_gene$gene)
matrix = matrix(nrow = length(unique(maf@data$Tumor_Sample_Barcode)), ncol = length(M_gene))
colnames(matrix) = M_gene
rownames(matrix) = unique(maf@data$Tumor_Sample_Barcode)

matrix[which(is.na(matrix))] = ""

i = 1
for(i in 1:ncol(matrix)){
  gene = colnames(matrix)[i]
  sub = maf@data[which(maf@data$Hugo_Symbol == gene),]
  samples = sub$Tumor_Sample_Barcode
  matrix[,gene][which(rownames(matrix) %in% samples)] = "MUT"
  matrix[,gene][which(!(rownames(matrix) %in% samples))] = "WT"
}

mat=as.data.frame(matrix)
mat$submitter_id=substr(rownames(mat),1,12)
clinical <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical")

data=merge(clinical,mat)

####TCGA microsatellite instability data#####
load("MSI.Rdata")
MSI$submitter_id=rownames(MSI)
data1=merge(data,MSI[,c(2,3)])
data1$Stage=ifelse(data1$ajcc_pathologic_stage %in% 
                     c("Stage I","Stage IA","Stage II",
                       "Stage IIA",  "Stage IIB",  "Stage IIC"),0,1)
data1$Tumor_Location=ifelse(data1$tissue_or_organ_of_origin %in% 
                              c("Ascending colon","Cecum","Hepatic flexure of colon",
                                "Transverse colon"),"Right",ifelse(
                                  data1$tissue_or_organ_of_origin=="Colon, NOS","Colon, NOS",
                                  "Left"
                                ))
data1$Stage=factor(data1$Stage,levels = c(1,0))
data1$MSI1=factor(data1$MSI1,levels = c("MSI-H","MSS"))
data1$gender=factor(data1$gender,levels = c("female","male"))
data1$Tumor_Location=factor(data1$Tumor_Location,levels = c("Left","Right","Colon, NOS"))
data2=data1[which(data1$Tumor_Location!="Colon, NOS"),]
data2$Tumor_Location=factor(data2$Tumor_Location,levels = c("Left","Right"))
x1=c("Stage","gender","MSI1","Tumor_Location")


M_gene
res <- data.frame(clin=character(), gene=character(), Cramers_v=numeric(), 
                  pvalue=numeric(),chi_squard=numeric(),OR=numeric(),
                  PR=numeric(),explose=character(),note=character(),
                  stringsAsFactors=FALSE)

results <- lapply(x1, function(x) {
  lapply(M_gene, function(y) {
    xx <- data2[, c(x, y)]
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
      try(re <- epi.2by2(table(xx$clin,xx$gene)+0.5),TRUE)
      pr=re$massoc.summary[1,2]
      or=re$massoc.summary[2,2]
      note1="adjust"
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
res_fdr[which(res_fdr$clin=="gender"),]$clin="Sex"
res_fdr[which(res_fdr$clin=="MSI1"),]$clin="MSI"
res_sig=res_fdr %>% filter(fdr<0.1)

write.table(res_fdr,"TCGA_clinical_gene_fdr_res.txt",row.names = F,quote = F,sep = "\t")
write.table(res_sig,"TCGA_clinical_gene_fdr_res_0.1_sig.txt",row.names = F,quote = F,
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
clinical_info=unique(res_sig$clin)
p=matrix(NA,nrow=length(clinical_info),ncol = length(gene),
         dimnames = list(clinical_info,gene))
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
              # row_split = left_list$A,
              row_title = NULL,
              # left_annotation = annotation_col1,
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
pdf("heatmap_TCGA_clin_gene_FDR_OR1.pdf",width = 13,height = 5)
draw(ht)
dev.off()
