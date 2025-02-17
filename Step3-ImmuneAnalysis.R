sig_matrix<-"mice.txt"
mixture_file<-"source_65.txt"
cibersort_file<-"re.txt"
colour=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(12,"Set3"))
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE,cibersort_file)
data<-t(res_cibersort[,-c(26,27,28)])
for(i in c(1:ncol(data))){data[,i]<-data[,i]/sum(data[,i])}

par(las = 1, mar = c(8,5,4,16),  # las用于定义轴标签的样式。mar用于设置图形的边距。
     mgp = c(3,0.1,0), cex.axis = 1.5)
 a1 <- barplot(data, col = colour,  # 设置颜色
               yaxt = "n", xaxt = "n",  # 不显示x、y坐标轴
               ylab = "Relative Percent",  # y轴标签
               cex.lab = 1.8)  # 标签文字大小
 a2 <- axis(side = 2, tick = FALSE, labels = FALSE)  # side用于指定要在绘图的哪一侧绘制轴。tick用于指定是否绘制刻度线和轴线。labels用于指定是否在刻度线处进行（数字）注释
 axis(side = 2, a2, paste0(a2*100, "%"))  # side = 2指定在绘图的左侧绘制轴
 axis(side = 1, a1, labels = FALSE)
 par(srt = 60, xpd = TRUE)  # srt表示以度为单位的字符串旋转。xpd用于设置图形打印裁剪方式
 text(x = a1, y = -0.02, labels = colnames(data), adj = 1, cex = 0.6)  # adj(adjustment)指定标签的x调整，0 表示左/下，1 表示右/上。cex(character expansion)，字体大小。
 # 添加图例
 par(srt = 0)  # srt表示graphics以度为单位的字符串旋转。
 ytick2 <- cumsum(data[, ncol(data)])  # cumsum(Cumulative Sums，累计总和)，将data的最后一列累计总和。
 ytick1 <- c(0, ytick2[-length(ytick2)])  # 将上述总和去掉最后一位，第一位前加0
 legend(x = par("usr")[2]*0.98, y = par("usr")[4],
        legend = rownames(data),
        col = colour, pch = 15, bty = "n", cex = 1.3)


data<-data[rownames(sample_Cluster),]
data<-cbind(sample_Cluster,data)
b<-melt(data,id.vars="Group")
colnames(b)<-c("Cluster","Immune Cell","Fraction")
p=ggboxplot(b, x="Immune Cell", y="Fraction", color = "Cluster", 
             ylab="Immune Cell Fraction",
             xlab="",
             legend.title="Immune Cell",
             palette = c("#fa6d1d","#0780cf"),
            width=0.6, add = "none")
p
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Cluster),
                         method="wilcox.test",
                         symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                         label = "p.signif")

td<-t(data)
xx_cor <- cor (td, method="spearman")
xx_p<-cor.mtest(td, method="spearman",conf.level = 0.95)
addcol <- colorRampPalette(c("#0073c2", "#dcdcdc", "#dc0000"))
pdf("CIBERSORT_COR.pdf",width=8,height=8)
corrplot(xx_cor,method = "color", col = addcol(100), tl.col = "black", tl.cex = 0.8, tl.srt = 45,tl.pos = "lt",p.mat = xx_p$p, diag = T, type = 'upper',sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.8,insig = 'label_sig', pch.col = 'grey20', order = 'AOE',cl.cex=0.6)
corrplot(xx_cor, method = "number", type = "lower",col = colorRampPalette(c("#0073c2", "#959595", "#dc0000"))(50),tl.col = "n", tl.cex = 0.8, tl.pos = "n",order = 'AOE',add = T,number.cex=0.5)
dev.off()

#mMCPcounter
Immune<- mMCPcounter.estimate(source_log, features =c("Gene.Symbol","ENSEMBL.ID","Probes")[1])
gsva_data<-Immune[,rownames(sample_Cluster)]
 library(limma)
 group_list <- factor(sample_Cluster$Group)
 design<-model.matrix(~0+group_list)
 colnames(design)<-levels(group_list)
 rownames(design)<-colnames(gsva_data)
 contrast.matrix<-makeContrasts('Cluster1 - Cluster2',levels=design)
 fit<-lmFit(gsva_data,design=design)
 fit2<-contrasts.fit(fit,contrast.matrix)
 fit2<-eBayes(fit2)
 alldiff=topTable(fit2,coef = 1,n = Inf)
 alldiff$lab = as.factor(ifelse(alldiff$P.Value>=0.05,"",ifelse(alldiff$P.Value>=0.01&alldiff$P.Value<0.05,"*",ifelse(alldiff$P.Value>=0.001&alldiff$P.Value<0.01,"**","***"))))
 alldiff$new <- paste(rownames(alldiff),alldiff$lab)
 alldiff <- alldiff[rownames(gsva_data),]
 rownames(gsva_data) <- alldiff$new
 pre_heatdata<-t(scale(t(gsva_data)))
 pre_heatdata[pre_heatdata>3]<-3#自己调，好看就行
 pre_heatdata[pre_heatdata<-3]<--3
annColors <- list()
annColors[['Group']] <- c('Cluster1'='#fa6d1d','Cluster2'='#0780cf')
pdf("mMCPCounter_heatmap.pdf",width=7,height=5)
 pheatmap(pre_heatdata,
 color = colorRampPalette(c("#0000ff",'#fae768',"#e30039"))(100),
 annotation_col = sample_Cluster,
 annotation_colors = annColors,
 cellheight=9,
 treeheight_row = 50,
 gaps_col = 52,#这里需要对应前面的分组样本数量
 show_rownames = T,
 show_colnames = F,
 cluster_rows = T,
 cluster_cols = F)
 dev.off()
write.table(Immune,"MCPCounter_result.txt",quote=F,sep="\t")
