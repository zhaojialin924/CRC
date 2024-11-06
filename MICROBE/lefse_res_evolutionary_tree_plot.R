####lefse######
####evolutionary tree####
library(MicrobiotaProcess)
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(tidytree)
library(ggstar)
library(forcats)
###prepare the data ####
data=read.table("micro_sig1_sig3_signature_group_lefse_count.txt",row.names = 1,header = T,sep = "\t")
seqtabfile <- system.file("extdata", "seqtab.nochim.rds", package="MicrobiotaProcess")
seqtab=read.delim("count/species.tsv",row.names = 1)
seqtab=as.data.frame(t(seqtab))
seqtab=seqtab[which(rownames(seqtab) %in% colnames(data)),]
colnames(seqtab)=gsub(" ","_",colnames(seqtab))
taxa <- read.delim("count/annotation.tsv")
taxa=taxa %>%
  select(-taxid) %>%
  distinct() %>%
  as.data.frame()
rownames(taxa)=taxa$species
rownames(taxa)=gsub(" ","_",rownames(taxa))
sampleda <- data[c(1),]
sampleda=as.data.frame(t(sampleda))
write.table(sampleda,"sampleda.txt",quote = F,sep = "\t")
data=data[-1,]
rownames(data)=gsub("_","__",rownames(data))
write.table(data,"evolutionary_tree.txt",quote = F,sep = "\t")

###creat the MPSE object####
file1="evolutionary_tree.txt"
sample.file="sampleda.txt"
mpse3 <- mp_import_metaphlan(profile=file1, mapfilename=sample.file)
mpse3 %<>%
  mp_diff_analysis(
    .abundance = Abundance,
    .group = Sig_group,
    first.test.alpha = 0.05,
    force = T
  )

###plot###
taxa.tree <- mpse3 %>% 
  mp_extract_tree(type="taxatree")
taxa.tree

##replace lefse_result#####
lefres=read.delim("02res/ours/micro_sig1_sig3_signature_group_lefse_count.res",sep = "\t",header = F)
lefres=lefres %>% subset(!is.na(V4))
colnames(lefres)=c("names","Logarith value","Group","LDA_value","P_value")

lefres$names=unlist(lapply(lefres$names,function(x) tail(strsplit(x,"\\.")[[1]],1)))
lefres=lefres %>% filter(P_value<0.05)
lefres$names=gsub("_","__",lefres$names)
lefres=lefres[grepl("^s_",lefres$names),]
lefres$Group=ifelse(lefres$Group=="C1-Sig1","SBS5","SBS6")
node=taxa.tree %>% 
  select(label, node,nodeClass, LDAupper, LDAmean, LDAlower, Sign_Sig_group, pvalue, fdr) %>% 
  dplyr::filter(!is.na(fdr))
colnames(lefres)[1]=colnames(node)[1]
lefres %>% as.tibble()
node2=merge(node,lefres,by="label")
extra=taxa.tree@extraInfo
head(node2)
node3=node2[,c(2,12,12,12,11,13,13)]
colnames(node3)=c("node","LDAupper" ,"LDAmean", 
                  "LDAlower", "Sign_Sig_group", "pvalue","fdr")
head(extra)
extra1=extra[which(!(extra$node %in% node3$node)),]
extra1[which(!is.na(extra1$Sign_Sig_group)),]$Sign_Sig_group=NA
extra2=extra[which((extra$node %in% node3$node)),]
extra3=merge(extra2[,c(1,2)],node3,by="node") %>% as.tibble()
extra4=rbind(extra1,extra3) %>% as.tibble()
extra4$pvalue=as.numeric(extra4$pvalue)
extra4$fdr=as.numeric(extra4$fdr)
extra4$Sign_Sig_group=factor(extra4$Sign_Sig_group,levels = rev(c("SBS6","SBS5")))

##replace###
taxa.tree@extraInfo=extra4
phylo1=taxa.tree@phylo
x=phylo1[["node.label"]]
x[c(4,7:12,14:18)]="Others"
phylo1[["node.label"]]=x
taxa.tree@phylo=phylo1

###plot###
p1 <- ggtree(
  taxa.tree,
  layout="radial",
  size = 0.3
) +
  geom_point(
    data = td_filter(!isTip),
    fill="white",
    size=1,
    shape=21
  )

p2 <- p1 +
  geom_hilight(
    data = td_filter(nodeClass == "Phylum"),
    mapping = aes(node = node, 
                  fill = label)
  )+
  scale_fill_manual(values = c("#f3b169","#3abf99","#c99bff","gray70"),
                    breaks=c("p__Bacillota","p__Bacteroidota",
                             "p__Pseudomonadota","Others"))
p2

# display the tip labels of taxa tree
p3 <- p2 + geom_tiplab(mapping = aes(subset= !is.na(Sign_Sig_group)),size=2, offset=.2)

# display the LDA of significant OTU.
p4 <- p3 +
  ggnewscale::new_scale("fill") +
  geom_fruit(
    geom = geom_col,
    mapping = aes(x = LDAmean,
                  fill = Sign_Sig_group,
                  subset = !is.na(Sign_Sig_group)
    ),
    orientation = "y",
    offset = .7,
    pwidth = 0.6,
    axis.params = list(axis = "x",
                       title = "Log10(LDA)",
                       title.height = 0.01,
                       title.size = 2,
                       text.size = 1.8,
                       vjust = 1),
    grid.params = list(linetype = 2)
  )

# display the significant (FDR) taxonomy after kruskal.test (default)
p5 <- p4 +
  ggnewscale::new_scale("size") +
  geom_point(
    data=td_filter(!is.na(Sign_Sig_group)),
    mapping = aes(size = -log10(pvalue),
                  fill = Sign_Sig_group,
    ),
    shape = 21,
  ) +
  scale_size_continuous(range=c(2, 4)) +
  scale_fill_manual(values=col[-2])+
  labs(fill="Group")

p=p5 + theme(
  legend.key.height = unit(0.3, "cm"),
  legend.key.width = unit(0.3, "cm"),
  legend.spacing.y = unit(0.02, "cm"),
  legend.text = element_text(size = 7,color = "black"),
  legend.title = element_text(size = 9,color = "black"),
)

ggsave(plot = p,"lefse_res_evolutionary_tree.pdf",width = 12,height = 12)
