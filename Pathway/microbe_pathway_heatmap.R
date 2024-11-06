#####热图######
library(dplyr)
library(ComplexHeatmap)
setwd("D:/Postgraduation/microbe_crc")
res=read.delim("IV_res_bacteria_new_pathway_80%_wilcox_fdr_0.05.txt",
               sep = "\t",header = T)
res1=res %>% filter(FC>1.25)
res_all=read.delim("IV_res_bacteria_new_pathway_80%_wilcox_all.txt")
res2=res_all[which(res_all$bacteria %in% unique(res1$bacteria)),]
res2=res2[which(res2$pathway %in% unique(res1$pathway)),]
res2$pathway=str_wrap(res2$pathway,width = 35)
ba=unique(res2$bacteria)
pa=unique(res2$pathway)
hm=matrix(NA,nrow=length(ba),ncol = length(pa),
          dimnames = list(ba,pa))

for(i in 1:nrow(res2)) {
  bacteria_idx <- match(res2$bacteria[i], ba)
  pathway_idx <- match(res2$pathway[i], pa)
  hm[bacteria_idx, pathway_idx] <- res2$FC[i]
}

hm_fdr=matrix(NA,nrow=length(ba),ncol = length(pa),
              dimnames = list(ba,pa))

for(i in 1:nrow(res2)) {
  bacteria_idx <- match(res2$bacteria[i], ba)
  pathway_idx <- match(res2$pathway[i], pa)
  hm_fdr[bacteria_idx, pathway_idx] <- res2$fdr[i]
}
hm_fdr=as.matrix(t(hm_fdr))

hm_fdr[hm_fdr<0.01]="**"
hm_fdr[hm_fdr>0.01 & hm_fdr<0.05]="*"
hm_fdr[hm_fdr>0.05]=""

left_list = res1 %>% select(category,pathway) %>% distinct()
rownames(left_list)=left_list$pathway
left_list = as.data.frame(left_list)
left_list$A = factor(left_list$category)
annotation_col1 = rowAnnotation(
  category = anno_block(gp=gpar(fill=c(2,4)),
                        labels=c("KEGG","METACYC"),
                        labels_gp=gpar(col="white",fontsize=14))
)
hm=as.matrix(t(hm))

ht <- Heatmap(hm,
              cluster_rows = F,
              cluster_columns = F,
              name = "FC",
              row_split = left_list$A,
              row_title = NULL,
              left_annotation = annotation_col1,
              col=circlize::colorRamp2(c(1,1.25,3),c("#619cff", "white", "#f8766d")),
              show_row_names = TRUE,
              show_column_names = TRUE,
              width = ncol(hm)*unit(14,"mm"),
              height = (nrow(hm)+0.8)*unit(14,"mm"),
              row_names_gp = gpar(fontsize=14),
              column_names_gp = gpar(fontsize=14),
              column_names_rot = 45,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(hm_fdr[i,j], x = x, y = y,gp = gpar(fontsize = 14,col="black"))},
              heatmap_legend_param = list(title = "** FDR < 0.01\n*  0.01 < FDR < 0.05\n\nFlodChange"))


pdf("heatmap_fdr0.05_fc_1.25.pdf",width = 10,height = 10)
draw(ht)
dev.off()
