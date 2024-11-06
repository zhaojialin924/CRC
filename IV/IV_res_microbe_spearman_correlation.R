###IV_res_bacteria correlation with other bacteria####
setwd("D:/Postgraduation/microbe_crc")

IV_res=read.csv("02res/ours/IV_res.csv")
data11=read.delim("count/species.filter.tsv",sep = "\t")
data11=as.data.frame(t(data11))
data11=data11/rowSums(data11)
data11=as.data.frame(t(data11))
data11$var=apply(data11,1,var)
data11=data11[order(-data11$var),]
res=data.frame()
for(ba1 in unique(IV_res$Bacteria) ){
  for (ba2 in colnames(data11)){
    xx=cor.test(data11[[ba1]],data11[[ba2]],method="spearman")
    cor1=xx$estimate
    pvalue1=xx$p.value
    res1=data.frame(ba1=ba1,ba2=ba2,cor=cor1,pvalue=pvalue1)
    res=rbind(res,res1)
  }
}
res_fdr=res %>%
  group_by(ba1) %>%
  mutate(fdr=p.adjust(pvalue)) %>%
  ungroup()

write.table(res_fdr,"IV_res_bacteria_spearman_cor_all.txt",row.names = F,
            sep = "\t",quote = F)

