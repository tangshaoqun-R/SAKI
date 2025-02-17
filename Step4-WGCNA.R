traitData=read.table("ORG_Score.txt",header=T,check.names=F,sep="\t",row.names=1)
sampleTree2 = hclust(dist(datExpr0), method = "average")
plot(sampleTree2)
enableWGCNAThreads()
powers = c(1:20)
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)


par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model
Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,cex=cex1,col="red")


softPower =sft$powerEstimate
adjacency = adjacency(datExpr0, power = softPower)
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04)
minModuleSize = 300
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = minModuleSize)
dynamicColors = labels2colors(dynamicMods)


plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
MEDissThres = 0.7
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors	#保留17个模块

plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
mergedColors = merge$colors
mergedMEs = merge$newMEs
plotDendroAndColors(geneTree,cbind(dynamicColors,mergedColors),c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE,hang = 0.03, addGuide = TRUE, guideHang = 0.05)



datTraits=traitData
datTraits<-datTraits[-31,,drop=F]#剔除一个离群样本88.P
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
 dim(textMatrix) = dim(moduleTraitCor)
 par(mar = c(5, 10, 3, 3))
 labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


#保留所有模块基因
for (mod in 1:nrow(table(moduleColors)))
 {  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr0)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0(modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
 }

