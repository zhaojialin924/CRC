setwd("D:/Postgraduation/microbe_crc/")
library(dplyr)
library(readxl)
library(egg)

####wilcox####
IV_res=read.csv("02res/ours/IV_res_sig_filter.csv")
IV_res=IV_res[,c(1:3)]
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
iv_ba=top50hvb[,which(colnames(top50hvb) %in% c(IV_res$Bacteria,"sample"))]
files=list.files("pircust2/",pattern = ".final.tsv")
res <- data.frame(category=character(), gene=character(), pathway=character(), W=numeric(), pvalue=numeric(),FC=numeric() ,stringsAsFactors=FALSE)

for (file in files) {
  name <- gsub(".final.tsv", "", file)
  abundance <- read.delim(paste0("pircust2/", file), row.names=1, header=TRUE, sep="\t")
  abundance <- log1p(abundance)
  selecter=rowSums(abundance>0)>78*0.8
  abundance <- abundance[selecter,]
  pathway <- as.data.frame(t(abundance))
  pathway <- as.data.frame(pathway)
  pats <- colnames(pathway)
  pathway$sample <- rownames(pathway)
  sigs=unique(IV_res$Bacteria)
  final1 <- merge(iv_ba, pathway, by="sample")
  results <- lapply(sigs, function(sig) {
    lapply(pats, function(pat1) {
      xx <- final1[, c(sig, pat1)]
      colnames(xx) <- c("sig", "pat")
      xx$sig=factor(ifelse(xx$sig==0,0,1))
      wilcox_test <- with(xx, wilcox.test(pat ~ sig, exact = FALSE))
      pvalue1=wilcox_test$p.value
      W=wilcox_test$statistic
      FC=mean(xx[which(xx$sig=="1"),]$pat)/mean(xx[which(xx$sig=="0"),]$pat)
      c(category=name, bacteria=sig, pathway=pat1, W=W, pvalue=pvalue1,FC=FC)
    })
  })
  
  res1 <- do.call(rbind, unlist(results, recursive=FALSE))
  res=rbind(res,res1)
}
res$pvalue=as.numeric(res$pvalue)
res$FC=as.numeric(res$FC)
res$W.W=as.numeric(res$W.W)
res_fdr=res %>%
  group_by(category,bacteria) %>%
  mutate(fdr=p.adjust(pvalue)) %>%
  ungroup()

write.table(res_fdr,"IV_res_bacteria_new_pathway_80%_wilcox_all.txt",row.names = F,
            sep = "\t",quote = F)

res_fdr=read.delim("IV_res_bacteria_new_pathway_80%_wilcox_all.txt")
res_fdr1=res_fdr %>%
  filter(category %in% c("Kegg","metacyc")) %>%
  filter(pvalue<0.05)
write.table(res_fdr1,"IV_res_bacteria_new_pathway_80%_wilcox_pvalue_sig_0904.txt",row.names = F,sep = "\t",quote = F)


sig=res_fdr2[which(res_fdr2$category=="metacyc"),]
for(i in c(1:nrow(sig))){
  xx=final1[c(sig$bacteria[i],sig$pathway[i])]
  colnames(xx)=c("ba","pathway")
  options(repr.plot.width = 3, repr.plot.height =4)
  p=ggplot(xx,aes(as.factor(ifelse(ba==0,0,1)),pathway,fill=as.factor(ifelse(ba==0,0,1))))+
    geom_boxplot()+
    xlab(sig$bacteria[i])+
    ylab(yulab.utils::str_wrap(sig$pathway[i],28))+
    theme_bw()+
    stat_compare_means(#aes(group =  group),
      comparisons = list(c("0","1")),
      size=2,
      label = "p.signif",
      method = "wilcox.test",
    )+
    scale_fill_manual(values = rev(col[-2]))+
    scale_x_discrete(breaks=c("0","1"),
                     labels= c("Absent","Present"))+
    theme(axis.text.x = element_text(color = "black",size=14,angle = 45,hjust = 1),
          axis.text.y = element_text(color="black",size=14),
          axis.title.x = element_text(size = 16, colour = "black"),
          axis.title.y = element_text(size = 16, colour = "black",hjust = 0.5,lineheight = 1.2),
          legend.position = "none",
    )+
    geom_jitter()
  ggsave(paste0("metacyc_",sig$bacteria[i],"_",sig$pathway[i],"_boxplot_plot5_0.05_2.pdf"),width = 3.8,height = 4.5,
         egg::set_panel_size(p,width=unit(1.8, "in"), height=unit(2.5, "in")))
  ggsave(paste0("metacyc_",sig$bacteria[i],"_",sig$pathway[i],"_boxplot_plot5_0.05_2.png"),width = 3.8,height = 4.5,
         egg::set_panel_size(p,width=unit(1.8, "in"), height=unit(2.5, "in")))
}
