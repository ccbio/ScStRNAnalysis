
###trajectory analysis of B cells
load("eso.merge.B.Rdata")
eso.merge.B$celltype_sub[eso.merge.B$celltype_top == "B cells" & eso.merge.B$RNA_snn_res.1 %in% c(0,8)] = "Activated B"
eso.merge.B <- NormalizeData(eso.merge.B)
eso.merge.B <- FindVariableFeatures(eso.merge.B)
eso.merge.B@meta.data$dataset = sapply(colnames(eso.merge.B),function(x){paste(strsplit(x,"_")[[1]][1],collapse="-")})
library(SeuratWrappers)
eso.merge.B <- RunFastMNN(object.list = SplitObject(eso.merge.B, split.by = 'dataset'))
eso.merge.B <- RunUMAP(eso.merge.B, reduction="mnn", dims = 1:30)
eso.merge.B <- FindNeighbors(eso.merge.B,  reduction="mnn", dims = 1:30)
eso.merge.B <- FindClusters(eso.merge.B, resolution=res)



Idents(eso.merge.B) = "celltype_sub"
sub = subset(eso.merge.B, downsample = 2000)


#  slingshot
library(slingshot)
library(SingleCellExperiment)
library(tidyverse)
library(RColorBrewer)

eso.merge.B$cluster = paste("B",eso.merge.B$RNA_snn_res.1,sep="")
eso.merge.B <- ScaleData(eso.merge.B)
scale.data <- eso.merge.B@assays$RNA@scale.data
scale.gene <- eso.merge.B$mnn.reconstructed@var.features
counts <- eso.merge.B@assays$RNA@counts
counts <- counts[scale.gene,]

sim <- SingleCellExperiment(assays= List(counts = counts))

umap=eso.merge.B@reductions$umap@cell.embeddings
colnames(umap)=c('UMAP-1','UMAP-2')

reducedDims(sim)=SimpleList(UMAP=umap)
# metadata
meta=eso.merge.B@meta.data

colData(sim)$sampleId=meta$orig.ident
colData(sim)$celltype_sub=meta$cluster


sim = slingshot(sim,clusterLabels = 'celltype_sub',
                reducedDim = 'UMAP',
                start.clus="B1",
                end.clus=NULL)
colnames(colData(sim))


colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)

plotcol <- colors[cut(sim$slingPseudotime_1,breaks=100)]

plotcol2<- colors[cut(sim$slingPseudotime_1,breaks=100)]
plotcol[is.na(plotcol)] <- plotcol2[is.na(plotcol)]
plotcol3<- colors[cut(sim$slingPseudotime_1,breaks=100)]
plotcol[is.na(plotcol)] <- plotcol3[is.na(plotcol)]

pdf("slingshot2.pdf",height = 6,width=7)
plot(reducedDims(sim)$UMAP,col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sim), lwd=2, col=brewer.pal(9,"Set1"))
legend("right", legend = paste0("lineage",1:3), col = unique(brewer.pal(3,"Set1")),inset=0.8, pch = 16)
#PlotTools::SpectrumLegend(x=13,y=6, horiz = FALSE,
#              legend = c("early","late"), palette= rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)),inset=0.05)
dev.off()
