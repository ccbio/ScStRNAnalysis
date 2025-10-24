#mature umap and plot for cluster cell type assignment
library(Seurat)
library(dplyr)
library(ggplot2)
library(paletteer)


#####1. Fibroblasts
#1.load data

load("results/Fibroblasts/esophagus.merge.obj.run.Fibroblasts.Rdata")
eso.merge.run.this = eso.merge.run.Fibroblasts ; rm(eso.merge.run.Fibroblasts)
dbl = read.table("results/scDblFinder.res.dbr0.76.txt",sep="\t",header=T)
eso.merge.run.this$scDblFinder.class = dbl$scDblFinder.class[match(colnames(eso.merge.run.this),dbl$cellbarcode)]


#3. assign cell types to cluster
eso.merge.run.this$celltype_sub = as.character(eso.merge.run.this$RNA_snn_res.1)
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(0,1,2,7,10)] = "MatrixCAFs"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(3,4,5,6,11,12,14)] = "InflammatoryCAFs"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(13)] = "VascularCAFs"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(8,18)] = "Pericytes"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(20)]  = "AnrigenPresentCAF"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(9)] = "MixedCAFs(iCAF&apCAF)"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(16)] = "MixedCAFs(iCAF&mCAF&tCAF)"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(17)] = "MixedCAFs(vCAF&tCAF)"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(15,19,21) |
                                  eso.merge.run.this$scDblFinder.class == "doublets"] = "doublets"
Idents(eso.merge.run.this) = "celltype_sub"
eso.merge.run.this = subset(eso.merge.run.this,idents=c("MatrixCAFs","InflammatoryCAFs","VascularCAFs","AnrigenPresentCAF",
                                                        "MixedCAFs(iCAF&apCAF)","MixedCAFs(iCAF&mCAF&tCAF)","MixedCAFs(vCAF&tCAF)","Pericytes"))
#https://r-charts.com/color-palettes/ paletteer::paletteer_d("ggthemes::Classic_20")
# a = as.character(paletteer::paletteer_d("ggthemes::Classic_20"))   ,  paste(a,collapse="\",\"")
colors = c("#a6cee3","#1f78b4","#b2df8a","#33a02c",'#fdbf6f',"#fb9a99","#e31a1c",'#ff7f00','#cab2d6','#6a3d9a',
           '#ffff99','#b15928')[1:length(levels(factor(eso.merge.run.this$celltype_sub)))]
names(colors) = c("MatrixCAFs","InflammatoryCAFs","VascularCAFs","AnrigenPresentCAF",
                  "MixedCAFs(iCAF&apCAF)","MixedCAFs(iCAF&mCAF&tCAF)","MixedCAFs(vCAF&tCAF)","Pericytes")

#4. DimPlot for subcluster cell type annotation
outdir = "results/Fibroblasts"
p3 = DimPlot(eso.merge.run.this, reduction = "umap",group.by = "RNA_snn_res.1", pt.size=0.1,label = T,repel = TRUE, label.size=6)
p3$layers[[1]]$mapping$alpha <- 0.8
p3 <- p3 + scale_alpha_continuous(range = 0.8, guide = "none") 
ggplot2::ggsave(paste0(outdir,"/clusters.dimplot.","res1.0",".pdf"), p3, width=6,height=6)
p3 = DimPlot(eso.merge.run.this, reduction = "umap",group.by = "celltype_sub", cols = colors,raster=F,
             #cells = colnames(eso.merge.run.this)[eso.merge.run.this$celltype_sub != "doublets"],
             pt.size=0.1,label = T,label.size=6, label.box=F, repel = F) + theme(legend.position = "none")
p3$layers[[1]]$mapping$alpha <- 0.7
p3 <- p3 + scale_alpha_continuous(range = 0.7, guide = "none") 
ggplot2::ggsave(paste0(outdir,"/clusters.dimplot.","celltype_sub",".pdf"), p3, width=6,height=6)

###5. Plot for distribution of subcluster SpecimenGroup 
## cell propotion in groups ESCC

eso.merge.run.this$tmp = eso.merge.run.this$SpecimenGroup
eso.merge.run.this$tmp[eso.merge.run.this$tmp %in% c("Adjacent esophagus",
                                                     "ESCC I","ESCC II","ESCC III","ESCC IVA") == FALSE] = "other"
eso.merge.run.this$tmp = factor(eso.merge.run.this$tmp,levels=c("Adjacent esophagus",
                                                                "ESCC I","ESCC II","ESCC III","ESCC IVA"))
Idents(eso.merge.run.this) = "celltype_sub"

p4 <- DimPlot(eso.merge.run.this, reduction = 'umap', split.by = 'tmp', ncol=3, pt.size=0.1, 
              label = F,repel = TRUE,cols = colors)
p4$layers[[1]]$mapping$alpha <- 0.5
p4 <- p4 + scale_alpha_continuous(range = 0.5, guide = "none")
ggplot2::ggsave(paste0(outdir,"/celltypes_sub_umap_ESCCprogression.pdf"), p4, width=12 ,height=12)

plotC <- reshape2::melt(table(eso.merge.run.this$tmp, eso.merge.run.this$celltype_sub))
colnames(plotC) <- c("Sample", "CellType","Number")
plotC$CellType = factor(plotC$CellType,levels=names(colors))
plotC$Prop = apply(plotC,1,function(x){ paste0(round(as.numeric(x[3]) / sum(plotC$Number[plotC$Sample == x[1]]),4)*100,"%")})
pC2 = ggplot(data = plotC, aes(x = Sample, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill")+
  scale_fill_manual(values=colors) + theme_bw()+
  #geom_text(aes(label = Prop), stat="identity",colour = "black",position = position_fill(vjust = 0.5)) +
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
  scale_y_continuous(labels = c(0,0.25,0.5,0.75,1))+ ####用来将y轴移动位置
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(linewidth = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.8))
pC2
ggplot2::ggsave(paste0(outdir,"/cellprop_sub_ESCCprogression.pdf"), pC2,width=6,height=6)



#6. plot dotplot for markers of cell types (expression level and proportion in cluster)

marker_cell = c(
  #"ACTA2","FAP","S100A4","FN1","COL1A1","VIM",  #Fibro
  "POSTN","CDH11","COL6A1","MMP11",  #mCAFs
  "CXCL12","CD34","C3","CFD","CXCL14", #iCAFs
  "ACTA2","MYH11","MCAM","BCAM", # vCAFs
  "ENO1","PGK1","GAPDH","PDPN", #tCAFs
  "HLA-DRA","HLA-DRB1","CD74",  #apCAFs
  "COX4I2","PDGFRB"   # Pericytes
)
Idents(eso.merge.run.this) = "celltype_sub"
p1= DotPlot(object = eso.merge.run.this, 
            features=unique(marker_cell),  col.min = 0, cluster.idents = TRUE,
            assay = "RNA") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) 
ggplot2::ggsave(paste0(outdir,"/markers.exp.dotplot.","celltype_sub",".pdf"), p1, width=10 ,height=4)








###########2. T cells
#1.load data

load("results/TNK/esophagus.merge.obj.run.TNK.Rdata")
eso.merge.run.this = eso.merge.run.TNK ; rm(eso.merge.run.TNK)
dbl = read.table("results/scDblFinder.res.dbr0.76.txt",sep="\t",header=T)
eso.merge.run.this$scDblFinder.class = dbl$scDblFinder.class[match(colnames(eso.merge.run.this),dbl$cellbarcode)]


#3. assign cell types to cluster
eso.merge.run.this$celltype_sub = as.character(eso.merge.run.this$RNA_snn_res.1)

eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(11,22,8,6,25)] = "Treg"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(12)] = "CD4 CTL"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(10)]  = "Th22"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(18)] = "CD4 Tex"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(26)] = "dnT"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(23)] = "gdT"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(7,15)] = "CD4 T Naive"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(14)] = "CD8 T Naive"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(5)] = "Th17"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(1,9,13)] = "CD8 Tex"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(20)] = "NKT"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(0,4)] = "CD8 Tem"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(2)] = "CD8 Tcm"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(16,7,3)] = "CD8 T cytotoxic"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(19,21,24,27,28) |
                                  eso.merge.run.this$scDblFinder.class == "doublets"] = "doublets"
Idents(eso.merge.run.this) = "celltype_sub"
eso.merge.run.this = subset(eso.merge.run.this,idents=c("Treg","CD4 CTL","Th22", "CD4 Tex",
                                                        "dnT","gdT","CD4 T Naive","CD8 T Naive","Th17","CD8 Tex","NKT",
                                                        "CD8 Tem","CD8 Tcm","CD8 T cytotoxic"))
colors = c("#4E79A7FF","#A0CBE8FF","#F28E2BFF","#FFBE7DFF","#59A14FFF","#8CD17DFF","#B6992DFF","#F1CE63FF","#499894FF",
           "#86BCB6FF","#E15759FF","#FF9D9AFF","#79706EFF","#BAB0ACFF","#D37295FF","#FABFD2FF","#B07AA1FF","#D4A6C8FF",
           "#9D7660FF")[1:length(levels(factor(eso.merge.run.this$celltype_sub)))]
names(colors) = c("Treg","CD4 CTL","Th22", "CD4 Tex",
                  "dnT","gdT","CD4 T Naive","CD8 T Naive","Th17","CD8 Tex","NKT",
                  "CD8 Tem","CD8 Tcm","CD8 T cytotoxic")


#4. DimPlot for subcluster cell type annotation
outdir = "results/TNK"
p3 = DimPlot(eso.merge.run.this, reduction = "umap",group.by = "RNA_snn_res.1", pt.size=0.1,label = T,repel = TRUE)
p3$layers[[1]]$mapping$alpha <- 0.5
p3 <- p3 + scale_alpha_continuous(range = 0.5, guide = "none") 
ggplot2::ggsave(paste0(outdir,"/clusters.dimplot.","res1.0",".pdf"), p3, width=6,height=6)
p3 = DimPlot(eso.merge.run.this, reduction = "umap",group.by = "celltype_sub", cols = colors,
             #cells = colnames(eso.merge.run.this)[eso.merge.run.this$celltype_sub != "doublets"],
             pt.size=0.1,label = T,label.size=3, label.box=T, repel = F) + theme(legend.position = "none")
p3$layers[[1]]$mapping$alpha <- 0.5
p3 <- p3 + scale_alpha_continuous(range = 0.5, guide = "none") 
ggplot2::ggsave(paste0(outdir,"/clusters.dimplot.","celltype_sub",".pdf"), p3, width=8,height=8)

###5. Plot for distribution of subcluster SpecimenGroup 
## cell propotion in groups ESCC

eso.merge.run.this$tmp = eso.merge.run.this$SpecimenGroup
eso.merge.run.this$tmp[eso.merge.run.this$tmp %in% c("Healthy esophagus","Adjacent esophagus",
                                                     "ESCC I","ESCC II","ESCC III","ESCC IVA") == FALSE] = "other"
eso.merge.run.this$tmp = factor(eso.merge.run.this$tmp,levels=c("Healthy esophagus","Adjacent esophagus",
                                                                "ESCC I","ESCC II","ESCC III","ESCC IVA"))
Idents(eso.merge.run.this) = "celltype_sub"

p4 <- DimPlot(eso.merge.run.this, reduction = 'umap', split.by = 'tmp', ncol=3, pt.size=0.1, 
              label = F,repel = TRUE,cols = colors)
#p4$layers[[1]]$mapping$alpha <- 0.5
#p4 <- p4 + scale_alpha_continuous(range = 0.5, guide = "none")
ggplot2::ggsave(paste0(outdir,"/celltypes_sub_umap_ESCCprogression.pdf"), p4, width=12 ,height=12)

plotC <- reshape2::melt(table(eso.merge.run.this$tmp, eso.merge.run.this$celltype_sub))
colnames(plotC) <- c("Sample", "CellType","Number")
plotC$CellType = factor(plotC$CellType,levels=names(colors))
plotC$Prop = apply(plotC,1,function(x){ paste0(round(as.numeric(x[3]) / sum(plotC$Number[plotC$Sample == x[1]]),4)*100,"%")})
pC2 = ggplot(data = plotC, aes(x = Sample, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill")+
  scale_fill_manual(values=colors) + theme_bw()+
  #geom_text(aes(label = Prop), stat="identity",colour = "black",position = position_fill(vjust = 0.5)) +
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
  scale_y_continuous(labels = c(0,0.25,0.5,0.75,1))+ ####用来将y轴移动位置
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(linewidth = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.8))
pC2
ggplot2::ggsave(paste0(outdir,"/cellprop_sub_ESCCprogression.pdf"), pC2,width=6,height=6)



#6. plot dotplot for markers of cell types (expression level and proportion in cluster)

marker_cell = c(
  "CD3E",
  "CD4","CD8A",
  "FOXP3","IL2RA","IL7R",  #Treg all IL7R-low
  "MKI67","TOP2A", #Treg Proliferating  22
  "GATA3",              #Treg GATA3   6,8,25
  "AHR",             #Treg GATA3 AHR  11
  "GZMA","GZMB","PRF1",   #CD4 CTL    12
  "AHR","CCR6","CXCL13",   #Th22  10
  "AHR","CCR6","CXCL13","TIGIT","TOX","PDCD1","LAG3","HAVCR2","CTLA4",  #Th22 exausted   18
  "CD4","CD8A","GZMA","GZMB","GZMK","GZMH","PRF1",   #dnT         23
  "TRDC","TRGC1","TRGC2",       #gdT    26
  "CCR7","SELL","IL2RA",      #CD4 T naive   15
  "CCR7","SELL","IL2RA","CD69",     #CD4 T naive CD69+  7 
  "RORC","IL17A","KLRB1",       #Th17 5
  "TIGIT","TOX","PDCD1","LAG3","HAVCR2","CTLA4","IL2","IFNG","TNF", "MKI67","TOP2A",       #CD8 Tex proliferating  13
  "TIGIT","TOX","PDCD1","LAG3","HAVCR2","CTLA4","IL2","IFNG","TNF", "NELL2","IL7R","NR4A1",   #CD8 Tex   NELL2+  IL7R-  1,9
  "CCL5","NKG7","CST7","GZMH","ZEB2",      #CD8 Tem2   ZEB2high vs CD8 Tem1  0
  "CCL5","NKG7","CST7","GZMH","TCF7",     #CD8 Tem1   TCF7high vs CD8 Tem2   4
  "ANXA1","IL7R","TRAC","GZMABHK","CD69",    #CD8 Tcm         2
  "GZMABHK","PRF1","CD69","IL2RA" ,     #CD8 T cytotoxic IL2RA+    16,17
  "GZMABHK","PRF1","CD69","CD69",      #CD8 T cytotoxic CD69+    3
  "NCAM1","NCR1","GNLY","FGFBP2","FCGR3A"  #NKT  20
)
Idents(eso.merge.run.this) = "celltype_sub"
p1= DotPlot(object = eso.merge.run.this, 
            features=unique(marker_cell), cluster.idents = TRUE,
            assay = "RNA") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) 
ggplot2::ggsave(paste0(outdir,"/markers.exp.dotplot.","celltype_sub",".pdf"), p1, width=15 ,height=8)



















###########3. B cells
#1.load data
load("results/B/esophagus.merge.obj.run.B.Rdata")
eso.merge.run.this = eso.merge.run.B #; rm(eso.merge.run.B)
dbl = read.table("results/scDblFinder.res.dbr0.76.txt",sep="\t",header=T)
eso.merge.run.this$scDblFinder.class = dbl$scDblFinder.class[match(colnames(eso.merge.run.this),dbl$cellbarcode)]


#3. assign cell types to cluster
eso.merge.run.this$celltype_sub = as.character(eso.merge.run.this$RNA_snn_res.1)
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(13,9)] = "GC B dark"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(12,7)] = "GC B light"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(1,3,6,10)] = "Naive B"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(2,4)]  = "Memory B"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(0,8)]  = "Activated B"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(5,11,14,15) |
                                  eso.merge.run.this$scDblFinder.class == "doublets"] = "doublets"
Idents(eso.merge.run.this) = "celltype_sub"
eso.merge.run.this = subset(eso.merge.run.this,idents=c("GC B dark","GC B light","Naive B","Memory B","Activated B"))
colors = c("#4E79A7FF","#A0CBE8FF","#F28E2BFF","#FFBE7DFF","#59A14FFF","#8CD17DFF","#B6992DFF","#F1CE63FF","#499894FF",
           "#86BCB6FF","#E15759FF","#FF9D9AFF","#79706EFF","#BAB0ACFF","#D37295FF","#FABFD2FF","#B07AA1FF","#D4A6C8FF",
           "#9D7660FF")[1:length(levels(factor(eso.merge.run.this$celltype_sub)))]
names(colors) = c("GC B dark","GC B light","Naive B","Activated B","Memory B")

Idents(eso.merge.run.B) = "seurat_clusters"
allmarkers <- FindAllMarkers(eso.merge.run.B, logfc.threshold = 0.2, min.pct = 0.2, only.pos = FALSE)
write.table(allmarkers,"/data_group/cunyupeng/songjing/project/01.xiaohua.tumors/geodata/esophagus/results/B/B_allmarkers.csv",sep=",",row.names=F,quote=F)
avg.b <- log1p(AverageExpression(eso.merge.run.B, verbose = FALSE)$RNA)
write.table(avg.b,"/data_group/cunyupeng/songjing/project/01.xiaohua.tumors/geodata/esophagus/results/B/B_allgene_pseudobulk.csv",sep=",",quote=F)

features = c("MKI67","TOP2A","MS4A1","CD38","TNFRSF13C","IGHD","EZH2","BACH2","BCL6","CD40","REL","CXCR4","CXCR5","CD83","CD86","MYC","IRF4")
p1= DotPlot(object = eso.merge.run.B, 
            features=features, cluster.idents = TRUE,  
            assay = "RNA") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) 
ggplot2::ggsave("cellmarkers_GC.B_res1.pdf", p1, width=7,height=5)

#4. DimPlot for subcluster cell type annotation
outdir = "results/B"
p3 = DimPlot(eso.merge.run.this, reduction = "umap",group.by = "RNA_snn_res.0.3", pt.size=0.1,label = T,repel = TRUE,label.size=6)
p3$layers[[1]]$mapping$alpha <- 0.8
p3 <- p3 + scale_alpha_continuous(range = 0.8, guide = "none") 
ggplot2::ggsave(paste0(outdir,"/clusters.dimplot.","res1.0",".pdf"), p3, width=6,height=6)
p3 = DimPlot(eso.merge.run.this, reduction = "umap",group.by = "celltype_sub", cols = colors,raster=F,
             #cells = colnames(eso.merge.run.this)[eso.merge.run.this$celltype_sub != "doublets"],
             pt.size=0.1,label = T,label.size=6, label.box=F, repel = F) + theme(legend.position = "none")
p3$layers[[1]]$mapping$alpha <- 0.7
p3 <- p3 + scale_alpha_continuous(range = 0.7, guide = "none") 
ggplot2::ggsave(paste0(outdir,"/clusters.dimplot.","celltype_sub",".pdf"), p3, width=6,height=6)

###5. Plot for distribution of subcluster SpecimenGroup 
## cell propotion in groups ESCC

eso.merge.run.this$tmp = eso.merge.run.this$SpecimenGroup
eso.merge.run.this$tmp[eso.merge.run.this$tmp %in% c("Healthy esophagus","Adjacent esophagus",
                                                     "ESCC I","ESCC II","ESCC III","ESCC IVA") == FALSE] = "other"
eso.merge.run.this$tmp = factor(eso.merge.run.this$tmp,levels=c("Healthy esophagus","Adjacent esophagus",
                                                                "ESCC I","ESCC II","ESCC III","ESCC IVA"))
Idents(eso.merge.run.this) = "celltype_sub"

p4 <- DimPlot(eso.merge.run.this, reduction = 'umap', split.by = 'tmp', ncol=3, pt.size=0.1, 
              label = F,repel = TRUE,cols = colors)
#p4$layers[[1]]$mapping$alpha <- 0.5
#p4 <- p4 + scale_alpha_continuous(range = 0.5, guide = "none")
ggplot2::ggsave(paste0(outdir,"/celltypes_sub_umap_ESCCprogression.pdf"), p4, width=12 ,height=12)

plotC <- reshape2::melt(table(eso.merge.run.this$tmp, eso.merge.run.this$celltype_sub))
colnames(plotC) <- c("Sample", "CellType","Number")
plotC$CellType = factor(plotC$CellType,levels=names(colors))
plotC$Prop = apply(plotC,1,function(x){ paste0(round(as.numeric(x[3]) / sum(plotC$Number[plotC$Sample == x[1]]),4)*100,"%")})
pC2 = ggplot(data = plotC, aes(x = Sample, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill")+
  scale_fill_manual(values=colors) + theme_bw()+
  #geom_text(aes(label = Prop), stat="identity",colour = "black",position = position_fill(vjust = 0.5)) +
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
  scale_y_continuous(labels = c(0,0.25,0.5,0.75,1))+ ####用来将y轴移动位置
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(linewidth = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.8))
pC2
ggplot2::ggsave(paste0(outdir,"/cellprop_sub_ESCCprogression.pdf"), pC2,width=6,height=6)



#6. plot dotplot for markers of cell types (expression level and proportion in cluster)

marker_cell = c(
  "MS4A1","CD19","CD79A",
  "IL4R","TCL1A","CXCR5","CXCR4","CD83","CD86","STMN1","LMO2",    #GCB light
  "TOP2A","MKI67",   #GCB dark
  "IGHM","IGHD","CD79A",  #naive B
  "CD80","CD86",     #B activated
  "BANK1","CD82"    #Memory B
)
Idents(eso.merge.run.this) = "celltype_sub"
p1= DotPlot(object = eso.merge.run.this, 
            features=unique(marker_cell), cluster.idents = TRUE,
            assay = "RNA") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) 
ggplot2::ggsave(paste0(outdir,"/markers.exp.dotplot.","celltype_sub",".pdf"), p1, width=7,height=3)


Idents(eso.merge.BPlasma) <- 'celltype_sub'
DefaultAssay(eso.merge.run.B) <- "RNA"
eso.merge.run.B = NormalizeData(eso.merge.run.B, verbose = FALSE)
avg.eso.merge.run.B <- AggregateExpression(object = eso.merge.run.B, group.by = c('RNA_snn_res.1'))$RNA
avg.eso.merge.run.B <-  log1p(AverageExpression(eso.merge.run.B, verbose = FALSE)$RNA)
write.table(avg.eso.merge.run.B,"B_allgene_pseudobulk.csv",sep=",",row.names=T,quote=F)

p1= DotPlot(object = eso.merge.BPlasma, 
            features=unique(marker_cell), cluster.idents = TRUE,
            assay = "RNA") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) 
ggplot2::ggsave("plasma.exp.res1.dotplot.pdf", p1, width=30,height=7)
p4 <- DimPlot(eso.merge.run.this, reduction = 'umap', pt.size=0.1, 
              label = F,repel = TRUE)
ggplot2::ggsave("B.umap.celltype_sub.pdf", p4, width=6 ,height=4)

allmarkers <- FindAllMarkers(eso.merge.run.this, logfc.threshold = 0.3, min.pct = 0.25, only.pos = TRUE)
write.table(allmarkers,"B_celltpye_sub_allmarkers.csv",sep=",",row.names=F,quote=F)











###########5. Plasma cells
#1.load data

load("results/Plasma/esophagus.merge.obj.run.Plasma.Rdata")
eso.merge.run.this = eso.merge.run.Plasma #; rm(eso.merge.run.Plasma)
dbl = read.table("results/scDblFinder.res.dbr0.76.txt",sep="\t",header=T)
eso.merge.run.this$scDblFinder.class = dbl$scDblFinder.class[match(colnames(eso.merge.run.this),dbl$cellbarcode)]


#3. assign cell types to cluster
eso.merge.run.this$celltype_sub = as.character(eso.merge.run.this$RNA_snn_res.0.1)
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(2,11,14)] = "IgG Plasma"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(0,9,10,13)] = "IgA Plasma"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(3,6,8,7)] = "Plasmablast"


eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(5,12,1,4,15) |
                                  eso.merge.run.this$scDblFinder.class == "doublets"] = "doublets"
Idents(eso.merge.run.this) = "celltype_sub"
eso.merge.run.this = subset(eso.merge.run.this,idents=c("IgG Plasma","IgA Plasma","Plasmablast"))
colors = c("#F28E2BFF","#59A14FFF","#E15759FF","#B6992DFF","#F1CE63FF","#499894FF",
           "#86BCB6FF","#FF9D9AFF","#79706EFF","#BAB0ACFF","#D37295FF","#FABFD2FF","#B07AA1FF","#D4A6C8FF",
           "#9D7660FF")[1:length(levels(factor(eso.merge.run.this$celltype_sub)))]
names(colors) = c("IgG Plasma","IgA Plasma","Plasmablast")

#4. DimPlot for subcluster cell type annotation
outdir = "results/Plasma"
p3 = DimPlot(eso.merge.run.this, reduction = "umap",group.by = "RNA_snn_res.1", pt.size=0.1,label = T,repel = TRUE)
p3$layers[[1]]$mapping$alpha <- 0.5
p3 <- p3 + scale_alpha_continuous(range = 0.5, guide = "none") 
ggplot2::ggsave(paste0(outdir,"/clusters.dimplot.","res1.0",".pdf"), p3, width=6,height=6)
p3 = DimPlot(eso.merge.run.this, reduction = "umap",group.by = "celltype_sub", cols = colors,
             #cells = colnames(eso.merge.run.this)[eso.merge.run.this$celltype_sub != "doublets"],
             pt.size=0.1,label = T,label.size=3, label.box=T, repel = F) + theme(legend.position = "none")
p3$layers[[1]]$mapping$alpha <- 0.5
p3 <- p3 + scale_alpha_continuous(range = 0.5, guide = "none") 
ggplot2::ggsave(paste0(outdir,"/clusters.dimplot.","celltype_sub",".pdf"), p3, width=8,height=8)

###5. Plot for distribution of subcluster SpecimenGroup 
## cell propotion in groups ESCC

eso.merge.run.this$tmp = eso.merge.run.this$SpecimenGroup
eso.merge.run.this$tmp[eso.merge.run.this$tmp %in% c("Adjacent esophagus",
                                                     "ESCC I","ESCC II","ESCC III","ESCC IVA") == FALSE] = "other"
eso.merge.run.this$tmp = factor(eso.merge.run.this$tmp,levels=c("Adjacent esophagus",
                                                                "ESCC I","ESCC II","ESCC III","ESCC IVA"))
Idents(eso.merge.run.this) = "celltype_sub"

p4 <- DimPlot(eso.merge.run.this, reduction = 'umap', split.by = 'tmp', ncol=3, pt.size=0.1, 
              label = F,repel = TRUE,cols = colors)
#p4$layers[[1]]$mapping$alpha <- 0.5
#p4 <- p4 + scale_alpha_continuous(range = 0.5, guide = "none")
ggplot2::ggsave(paste0(outdir,"/celltypes_sub_umap_ESCCprogression.pdf"), p4, width=12 ,height=12)

plotC <- reshape2::melt(table(eso.merge.run.this$tmp, eso.merge.run.this$celltype_sub))
colnames(plotC) <- c("Sample", "CellType","Number")
plotC$CellType = factor(plotC$CellType,levels=names(colors))
plotC$Prop = apply(plotC,1,function(x){ paste0(round(as.numeric(x[3]) / sum(plotC$Number[plotC$Sample == x[1]]),4)*100,"%")})
pC2 = ggplot(data = plotC, aes(x = Sample, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill")+
  scale_fill_manual(values=colors) + theme_bw()+
  #geom_text(aes(label = Prop), stat="identity",colour = "black",position = position_fill(vjust = 0.5)) +
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
  scale_y_continuous(labels = c(0,0.25,0.5,0.75,1))+ ####用来将y轴移动位置
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(linewidth = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.8))
pC2
ggplot2::ggsave(paste0(outdir,"/cellprop_sub_ESCCprogression.pdf"), pC2,width=6,height=6)



#6. plot dotplot for markers of cell types (expression level and proportion in cluster)

marker_cell = c(
  "MZB1","PRDM1",
  "IGHG1","IGHG2",   #IgG Plasma
  "IGHA2",    #IgA Plasma
  "TNFRSF17","POU2AF1","DERL3"   #Plasmablast
  
)
Idents(eso.merge.run.this) = "celltype_sub"
p1= DotPlot(object = eso.merge.run.this, 
            features=unique(marker_cell), cluster.idents = TRUE,
            assay = "RNA") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) 
ggplot2::ggsave(paste0(outdir,"/markers.exp.dotplot.","celltype_sub",".pdf"), p1, width=5 ,height=3)










###########6. Endothelial cells
#1.load data

load("results/Endothelial/esophagus.merge.obj.run.Endothelial.Rdata")
eso.merge.run.this = eso.merge.run.Endothelial #; rm(eso.merge.run.Endothelial)
dbl = read.table("results/scDblFinder.res.dbr0.76.txt",sep="\t",header=T)
eso.merge.run.this$scDblFinder.class = dbl$scDblFinder.class[match(colnames(eso.merge.run.this),dbl$cellbarcode)]


#3. assign cell types to cluster
eso.merge.run.this$celltype_sub = as.character(eso.merge.run.this$RNA_snn_res.1)
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(9)] = "Lymphatic Endo"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(1,8,2,6,14)] = "Hypoxia&Inflammatory Endo"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(3,4)] = "Microvascular Endo"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(0,7)] = "Pluripotent Endo"

eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(10,11,13,16,15) |
                                  eso.merge.run.this$scDblFinder.class == "doublets"] = "doublets"
Idents(eso.merge.run.this) = "celltype_sub"
eso.merge.run.this = subset(eso.merge.run.this,idents=c("Lymphatic Endo","Hypoxia&Inflammatory Endo",
                                                        "Microvascular Endo","Pluripotent Endo"))
colors = c("#F28E2BFF","#D37295FF","#59A14FFF","#499894FF","#8CD17DFF","#B6992DFF","#F1CE63FF",
           "#86BCB6FF","#E15759FF","#FF9D9AFF","#79706EFF","#BAB0ACFF","#FABFD2FF","#B07AA1FF","#D4A6C8FF",
           "#9D7660FF")[1:length(levels(factor(eso.merge.run.this$celltype_sub)))]
names(colors) = c("Lymphatic Endo","Hypoxia&Inflammatory Endo",
                  "Microvascular Endo","Pluripotent Endo")

#4. DimPlot for subcluster cell type annotation
outdir = "results/Endothelial"
p3 = DimPlot(eso.merge.run.this, reduction = "umap",group.by = "RNA_snn_res.1", pt.size=0.1,label = T,repel = TRUE,label.size = 6)
p3$layers[[1]]$mapping$alpha <- 0.8
p3 <- p3 + scale_alpha_continuous(range = 0.8, guide = "none") 
ggplot2::ggsave(paste0(outdir,"/clusters.dimplot.","res1.0",".pdf"), p3, width=6,height=6)
p3 = DimPlot(eso.merge.run.this, reduction = "umap",group.by = "celltype_sub", cols = colors,raster=F,
             #cells = colnames(eso.merge.run.this)[eso.merge.run.this$celltype_sub != "doublets"],
             pt.size=0.1,label = T,label.size=6, label.box=F, repel = F) + theme(legend.position = "none")
p3$layers[[1]]$mapping$alpha <- 0.7
p3 <- p3 + scale_alpha_continuous(range = 0.7, guide = "none") 
ggplot2::ggsave(paste0(outdir,"/clusters.dimplot.","celltype_sub",".pdf"), p3, width=6,height=6)

###5. Plot for distribution of subcluster SpecimenGroup 
## cell propotion in groups ESCC

eso.merge.run.this$tmp = eso.merge.run.this$SpecimenGroup
eso.merge.run.this$tmp[eso.merge.run.this$tmp %in% c("Adjacent esophagus",
                                                     "ESCC I","ESCC II","ESCC III","ESCC IVA") == FALSE] = "other"
eso.merge.run.this$tmp = factor(eso.merge.run.this$tmp,levels=c("Adjacent esophagus",
                                                                "ESCC I","ESCC II","ESCC III","ESCC IVA"))
Idents(eso.merge.run.this) = "celltype_sub"

p4 <- DimPlot(eso.merge.run.this, reduction = 'umap', split.by = 'tmp', ncol=3, pt.size=0.1, 
              label = F,repel = TRUE,cols = colors)
#p4$layers[[1]]$mapping$alpha <- 0.5
#p4 <- p4 + scale_alpha_continuous(range = 0.5, guide = "none")
ggplot2::ggsave(paste0(outdir,"/celltypes_sub_umap_ESCCprogression.pdf"), p4, width=12 ,height=12)

plotC <- reshape2::melt(table(eso.merge.run.this$tmp, eso.merge.run.this$celltype_sub))
colnames(plotC) <- c("Sample", "CellType","Number")
plotC$CellType = factor(plotC$CellType,levels=names(colors))
plotC$Prop = apply(plotC,1,function(x){ paste0(round(as.numeric(x[3]) / sum(plotC$Number[plotC$Sample == x[1]]),4)*100,"%")})
pC2 = ggplot(data = plotC, aes(x = Sample, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill")+
  scale_fill_manual(values=colors) + theme_bw()+
  #geom_text(aes(label = Prop), stat="identity",colour = "black",position = position_fill(vjust = 0.5)) +
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
  scale_y_continuous(labels = c(0,0.25,0.5,0.75,1))+ ####用来将y轴移动位置
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(linewidth = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.8))
pC2
ggplot2::ggsave(paste0(outdir,"/cellprop_sub_ESCCprogression.pdf"), pC2,width=6,height=6)


#6. plot dotplot for markers of cell types (expression level and proportion in cluster)

marker_cell = c(
  "PECAM1","VMF","ACVRL1",
  "PROX1","LYVE1","FLT4","PDPN",  #lymphatic endo
  "HIF1A","SELE","ICAM1","VCAM1","IL6",  # Inflammatory and hypoxia endo
  "CLDN5","JAM1","VWF","CDH5","CD34",    #microvascular endo VWF low CDH5- CD34-
  "TEK","KDR",   #endo progenitor
  "CXCL12","HIF1A","ACTB","VIM","PLVAP"   #多能内皮细胞
)
Idents(eso.merge.run.this) = "celltype_sub"
p1= DotPlot(object = eso.merge.run.this, 
            features=unique(marker_cell), cluster.idents = TRUE,
            assay = "RNA") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) 
ggplot2::ggsave(paste0(outdir,"/markers.exp.dotplot.","celltype_sub",".pdf"), p1, width=15 ,height=8)




























###########7. Myeloid cells
#1.load data

load("results/Myeloid/esophagus.merge.obj.run.Myeloid.Rdata")
eso.merge.run.this = eso.merge.run.Myeloid #; rm(eso.merge.run.Myeloid)
dbl = read.table("results/scDblFinder.res.dbr0.76.txt",sep="\t",header=T)
eso.merge.run.this$scDblFinder.class = dbl$scDblFinder.class[match(colnames(eso.merge.run.this),dbl$cellbarcode)]


#3. assign cell types to cluster
eso.merge.run.this$celltype_sub = as.character(eso.merge.run.this$RNA_snn_res.1)
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(3,5,20)] = "cDC2"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(16)] = "cDC1"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(17)] = "tDC"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(9)] = "pDC"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(15)] = "Megakaryocyte"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(0,1,6,21,22)] = "CD16 Mono"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(11)] = "M2 Macrophage"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(8,18)] = "TRNs"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(14)] = "M1 Macrophage"
eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(2,4,7)] = "CD14 Mono"

eso.merge.run.this$celltype_sub[eso.merge.run.this$celltype_sub %in% c(19,10,12,13) |
                                  eso.merge.run.this$scDblFinder.class == "doublets"] = "doublets"
Idents(eso.merge.run.this) = "celltype_sub"
eso.merge.run.this = subset(eso.merge.run.this,idents=c("cDC2","cDC1","tDC","pDC","Megakaryocyte","CD16 Mono","CD14 Mono",
                                                        "M1 Macrophage","M2 Macrophage","TRNs"))
colors = c("#F28E2BFF","#D37295FF","#59A14FFF","#499894FF","#8CD17DFF","#B6992DFF","#F1CE63FF",
           "#86BCB6FF","#E15759FF","#FF9D9AFF","#FABFD2FF","#B07AA1FF","#79706EFF","#BAB0ACFF","#D4A6C8FF",
           "#9D7660FF")[1:length(levels(factor(eso.merge.run.this$celltype_sub)))]
names(colors) = c("cDC2","cDC1","tDC","pDC","Megakaryocyte","CD16 Mono","CD14 Mono",
                  "M1 Macrophage","M2 Macrophage","TRNs")

#4. DimPlot for subcluster cell type annotation
outdir = "results/Myeloid"
p3 = DimPlot(eso.merge.run.this, reduction = "umap",group.by = "RNA_snn_res.1", pt.size=0.1,label = T,repel = TRUE,raster=T,label.size = 8)
#p3$layers[[1]]$mapping$alpha <- 0.8
#p3 <- p3 + scale_alpha_continuous(range = 0.8, guide = "none") 
ggplot2::ggsave(paste0(outdir,"/clusters.dimplot.","res1.0",".pdf"), p3, width=7,height=6)
p3 = DimPlot(eso.merge.run.this, reduction = "umap",group.by = "celltype_sub", cols = colors,raster=F,
             #cells = colnames(eso.merge.run.this)[eso.merge.run.this$celltype_sub != "doublets"],
             pt.size=0.1,label = T,label.size=6, label.box=F, repel = F) + theme(legend.position = "none")
p3$layers[[1]]$mapping$alpha <- 0.5
p3 <- p3 + scale_alpha_continuous(range = 0.5, guide = "none") 
ggplot2::ggsave(paste0(outdir,"/clusters.dimplot.","celltype_sub",".pdf"), p3, width=6,height=6)

###5. Plot for distribution of subcluster SpecimenGroup 
## cell propotion in groups ESCC

eso.merge.run.this$tmp = eso.merge.run.this$SpecimenGroup
eso.merge.run.this$tmp[eso.merge.run.this$tmp %in% c("Healthy esophagus","Adjacent esophagus",
                                                     "ESCC I","ESCC II","ESCC III","ESCC IVA") == FALSE] = "other"
eso.merge.run.this$tmp = factor(eso.merge.run.this$tmp,levels=c("Healthy esophagus","Adjacent esophagus",
                                                                "ESCC I","ESCC II","ESCC III","ESCC IVA"))
Idents(eso.merge.run.this) = "celltype_sub"

p4 <- DimPlot(eso.merge.run.this, reduction = 'umap', split.by = 'tmp', ncol=3, pt.size=0.1, 
              label = F,repel = TRUE,cols = colors)
#p4$layers[[1]]$mapping$alpha <- 0.5
#p4 <- p4 + scale_alpha_continuous(range = 0.5, guide = "none")
ggplot2::ggsave(paste0(outdir,"/celltypes_sub_umap_ESCCprogression.pdf"), p4, width=12 ,height=12)

plotC <- reshape2::melt(table(eso.merge.run.this$tmp, eso.merge.run.this$celltype_sub))
colnames(plotC) <- c("Sample", "CellType","Number")
plotC$CellType = factor(plotC$CellType,levels=names(colors))
plotC$Prop = apply(plotC,1,function(x){ paste0(round(as.numeric(x[3]) / sum(plotC$Number[plotC$Sample == x[1]]),4)*100,"%")})
pC2 = ggplot(data = plotC, aes(x = Sample, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill")+
  scale_fill_manual(values=colors) + theme_bw()+
  #geom_text(aes(label = Prop), stat="identity",colour = "black",position = position_fill(vjust = 0.5)) +
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
  scale_y_continuous(labels = c(0,0.25,0.5,0.75,1))+ ####用来将y轴移动位置
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(linewidth = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.8))
pC2
ggplot2::ggsave(paste0(outdir,"/cellprop_sub_ESCCprogression.pdf"), pC2,width=6,height=6)



#6. plot dotplot for markers of cell types (expression level and proportion in cluster)

marker_cell = c(
  "ITGAM","CD68","CD14","FUT4",
  "CD1C",               #cDC2 
  "XCR1", "CLEC9A",  #cDC1
  "CCR7","CD274",     #tDC  CCR7 CD274   6,8,25
  "IL3RA",             #pDC  11
  "CCL3L1",   #Megakayocyte    12
  "CD14",    #CD14 Mono
  "CD14","FCGR3A",   #CD16 Mono
  "MRC1","CD163",  #M2 macrophage
  "CXCR2","G0S2","FCGR3B","CSF3R",   #TRNs
  "CXCL9","CXCL10"     #M1 Macrophage  CXCL9 CXCL10
  
)
Idents(eso.merge.run.this) = "celltype_sub"
p1= DotPlot(object = eso.merge.run.this, 
            features=unique(marker_cell), cluster.idents = TRUE, col.min = 0,
            assay = "RNA") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) 
ggplot2::ggsave(paste0(outdir,"/markers.exp.dotplot.","celltype_sub",".pdf"), p1, width=15 ,height=8)

