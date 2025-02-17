data_sort <- DEG %>%
    arrange(desc(logFC))
gene_list <- data_sort$logFC
names(gene_list) <- rownames(data_sort)

library(msigdbr)
reactome_gmt<-read.gmt("m2.cp.reactome.v2024.1.Mm.symbols.gmt")
gobp_gmt<-read.gmt("m5.go.bp.v2024.1.Mm.symbols.gmt")
react_symbol<-clusterProfiler::GSEA(geneList = gene_list,pvalueCutoff = 0.05,pAdjustMethod = "BH",TERM2GENE=reactome_gmt)
gobp_symbol<-clusterProfiler::GSEA(geneList = gene_list,pvalueCutoff = 0.05,pAdjustMethod = "BH",TERM2GENE=gobp_gmt)

#Reactome
react_symbol_cut <- react_symbol[react_symbol$pvalue<0.05 & abs(react_symbol$NES)>1]
react_symbol_cut<-react_symbol_cut[order(react_symbol_cut$NES,decreasing =T),]

up_react <- head(react_symbol_cut,6)
down_react<-tail(react_symbol_cut,6)

#GOBP
gobp_symbol_cut <- gobp_symbol[gobp_symbol$pvalue<0.05 & abs(gobp_symbol$NES)>1]
gobp_symbol_cut<-gobp_symbol_cut[order(gobp_symbol_cut$NES,decreasing =T),]

up_gobp <- head(gobp_symbol_cut,6)
down_gobp<-tail(gobp_symbol_cut,6)

#山峦图
pdf("GSEA_GOBP_DOWN.pdf",width=8,height=7)
gseaplot2(gobp_symbol,
        down_gobp$ID,#富集的ID编号
        title = "Enriched in GO Biological Process",#标题
        color = "red",#GSEA线条颜色
        base_size = 15,#基础字体大小
        rel_heights = c(1.5, 0.5, 0.5),#副图的相对高度
        subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
        ES_geom = "line",#enrichment score用线还是用点"dot"
        pvalue_table = F)
dev.off()


#GSVA分析
reactome_geneset<-readLines("m2.cp.reactome.v2024.1.Mm.symbols.gmt")
reactome_geneset<-strsplit(reactome_geneset, "\t")
names(reactome_geneset) <- vapply(reactome_geneset, function(y) y[1], character(1))
reactome_geneset <- lapply(reactome_geneset , "[", -c(1:2))

gobp_geneset<-readLines("m5.go.bp.v2024.1.Mm.symbols.gmt")
gobp_geneset<-strsplit(gobp_geneset, "\t")
names(gobp_geneset) <- vapply(gobp_geneset, function(y) y[1], character(1))
gobp_geneset <- lapply(gobp_geneset , "[", -c(1:2))

react_gsva<-gsva(disea_exp,reactome_geneset,kcdf = "Poisson", min.sz = 10)
gobp_gsva<-gsva(disea_exp,gobp_geneset,kcdf = "Poisson", min.sz = 10)

#差异分析GSVA结果
react_DEG <- TCGAanalyze_DEA(
  mat1 = react_gsva[, rownames(result_cluster[which(result_cluster$Group=="Cluster1"),,drop=F])],
  mat2 = react_gsva[, rownames(result_cluster[which(result_cluster$Group=="Cluster2"),,drop=F])],
  metadata = FALSE,
  pipeline = "limma",
  Cond1type = "Cluster1",
  Cond2type = "Cluster2",
  fdr.cut = 0.05,
  logFC.cut = 0,
)

gobp_DEG <- TCGAanalyze_DEA(
  mat1 = gobp_gsva[, rownames(result_cluster[which(result_cluster$Group=="Cluster1"),,drop=F])],
  mat2 = gobp_gsva[, rownames(result_cluster[which(result_cluster$Group=="Cluster2"),,drop=F])],
  metadata = FALSE,
  pipeline = "limma",
  Cond1type = "Cluster1",
  Cond2type = "Cluster2",
  fdr.cut = 0.05,
  logFC.cut = 0,
)

react_DEG<-react_DEG[order(react_DEG$logFC,decreasing=T),]
top10_react<-head(rownames(react_DEG),10)
down10_react<-tail(rownames(react_DEG),10)
choose_pathway<-c(top10_react,down10_react)
choose_matrix <- react_gsva[choose_pathway,]
heat_matrix <- t(scale(t(choose_matrix)))
ann_colors=list(Group=c(Cluster1='#fa6d1d',Cluster2='#0780cf'))
heat_matrix<-heat_matrix[,rownames(result_cluster)]
pheatmap(heat_matrix,color = colorRampPalette(c("#D2EEF9","#C01858"))(100),border_color = NA,fontsize = 10,show_rownames = T,annotation_col=result_cluster,annotation_colors=ann_colors,cluster_cols = FALSE,show_colnames=F,cellwidth=3,cellheight=13)

gobp_DEG<-gobp_DEG[order(gobp_DEG$logFC,decreasing=T),]
top10_gobp<-head(rownames(gobp_DEG),10)
down10_gobp<-tail(rownames(gobp_DEG),10)
choose_pathway<-c(top10_gobp,down10_gobp)
choose_matrix <- gobp_gsva[choose_pathway,]
heat_matrix <- t(scale(t(choose_matrix)))
ann_colors=list(Group=c(Cluster1='#fa6d1d',Cluster2='#0780cf'))
heat_matrix<-heat_matrix[,rownames(result_cluster)]
pheatmap(heat_matrix,color = colorRampPalette(c("#D2EEF9","#C01858"))(100),border_color = NA,fontsize = 10,show_rownames = T,annotation_col=result_cluster,annotation_colors=ann_colors,cluster_cols = FALSE,show_colnames=F,cellwidth=3,cellheight=13)
