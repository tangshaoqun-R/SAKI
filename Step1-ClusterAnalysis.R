library(ConsensusClusterPlus)
out="clusterResult"
results=ConsensusClusterPlus(NET_disea_exp,maxK=6,reps=1000,pItem=0.8,pFeature=1,title=out,clusterAlg="km",distance="euclidean",seed=2024,plot="pdf",writeTable=TRUE)
result_cluster<-read.csv("clusterResult/clusterResult.k=2.consensusClass.csv",check.names=F,row.names=1,header=F)
colnames(result_cluster)<-"Group"
list<-result_cluster$Group %>% factor(., levels=c("Cluster1","Cluster2"),ordered=F)
list <- model.matrix(~factor(list)+0)
colnames(list) <- c("Cluster1","Cluster2")
df.fit <- lmFit(disea_exp, list)
df.matrix <- makeContrasts(Cluster1-Cluster2,levels=list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
tempOutput <- topTable(fit,n = Inf, adjust = "fdr")
head(tempOutput)
nrDEG = na.omit(tempOutput)
#write.table(nrDEG,"NET_limmaOut.txt",quote=F,sep="\t")

#差异分析可视化
sort_DEG<-arrange(DEG,logFC)
sort_DEG<-sort_DEG[which(sort_DEG$adj.P.Val<0.05),]
top10<-tail(rownames(sort_DEG),10)
down10<-head(rownames(sort_DEG),10)
choose_gene<-c(top10,down10)
choose_matrix <- disea_exp[choose_gene,]
heat_matrix <- t(scale(t(choose_matrix)))
ann_colors=list(Group=c(Cluster1='#fa6d1d',Cluster2='#0780cf'))
pheatmap(heat_matrix,color = colorRampPalette(c("#D2EEF9","#C01858"))(100),border_color = NA,fontsize = 10,show_rownames = T,annotation_col=result_cluster,annotation_colors=ann_colors,cluster_cols = FALSE,show_colnames=F)

#富集分析
go <- enrichGO(gene = gene_id$ENTREZID,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05,qvalueCutoff=0.05, keyType="ENTREZID",
                ont="all",pAdjustMethod = "BH",readable = TRUE)
write.table(go,file="GO.txt",sep="\t",quote=F,row.names = F)  
go_enrich<-read.table("GO_top10.txt",header=T,check.names=F,sep="\t")
go_enrich$term <- paste(go_enrich$ID, go_enrich$Description, sep = ': ') #将ID与Description合并成新的一列
go_enrich$term <- factor(go_enrich$term, levels = go_enrich$term,ordered = T) #转成因子，防止重新排列
p1 <- ggplot(go_enrich,
       aes(x=term,y=Count, fill=ONTOLOGY)) + #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) +  #柱状图填充颜色
  coord_flip() +  #让柱状图变为纵向
  xlab("GO term") +  #x轴标签
  ylab("Gene Number") +  #y轴标签
  labs(title = "GO Terms Enrich")+  #设置标题
  theme_bw()
p1

#KEGG
kegg <- enrichKEGG(gene = gene_id$ENTREZID,keyType="kegg",organism = "mmu", pvalueCutoff = 0.05,qvalueCutoff=0.05)
write.table(kegg,file="KEGG.txt",sep="\t",quote=F,row.names = F)
kegg_enrich<-read.table("KEGG_top20.txt",header=T,check.names=F,sep="\t")
#ggplot可视化
kegg_enrich$term <- paste(kegg_enrich$ID, kegg_enrich$Description, sep = ': ') #将ID与Description合并成新的一列
kegg_enrich$term <- factor(kegg_enrich$term, levels = kegg_enrich$term,ordered = T) #转成因子，防止重新排列
p <- ggplot(data = kegg_enrich,mapping = aes(x = Count,y = term))+
  geom_point(aes(color= -log10(p.adjust),size = Count)) +
  scale_colour_gradient(high = '#C01858',low = '#AEE2FA') +
  theme_bw()+
  labs(title = 'KEGG_enrich',
       x = 'Count',
       y = 'Description')