
setwd("path/to/workdirectory")

#GSE160269
library(Seurat)
CD45neg <- read.table("GSE160269/GSE160269_CD45neg_UMIs.txt.gz")
CD45neg <- CreateSeuratObject(counts=CD45neg,project="CD45neg")
CD45neg <- PercentageFeatureSet(CD45neg, pattern = "^MT-", col.name = "percent.mt")
CD45neg@meta.data$orig.ident = sapply(colnames(CD45neg),
                                      function(x){paste(strsplit(x,"\\.")[[1]][c(1,2)],collapse="-")})
CD45pos <- read.table("GSE160269/GSE160269_CD45pos_UMIs.txt.gz")
CD45pos2 <- CreateSeuratObject(counts=CD45pos,project="CD45pos")
CD45pos <- PercentageFeatureSet(CD45pos, pattern = "^MT-", col.name = "percent.mt")
CD45pos@meta.data$orig.ident = sapply(colnames(CD45pos),function(x){paste(strsplit(x,"\\.")[[1]][c(1,2)],collapse="-")})
GSE160269 = merge(x=CD45neg,y=CD45pos)



#GSE199654
counts = read.table("GSE199654/GSE199654_scTDN_UMI_matrix_epithelial_cells.txt.gz")
counts <- CreateSeuratObject(counts=counts)
counts <- PercentageFeatureSet(counts, pattern = "^MT-", col.name = "percent.mt")
GSE199654 = counts

#GSE221561
myRead10X = function (data.dir, sample = sample, gene.column = 2, cell.column = 1, unique.features = TRUE, 
                      strip.suffix = FALSE) 
{
  full.data <- list()
  has_dt <- requireNamespace("data.table", quietly = TRUE) && 
    requireNamespace("R.utils", quietly = TRUE)
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, paste(sample,"barcodes.tsv",sep="_"))
    gene.loc <- file.path(run, paste(sample,"genes.tsv",sep="_"))
    features.loc <- file.path(run, paste(sample,"features.tsv.gz",sep="_"))
    matrix.loc <- file.path(run, paste(sample,"matrix.mtx",sep="_"))
    pre_ver_3 <- file.exists(gene.loc)
    if (!pre_ver_3) {
      addgz <- function(s) {
        return(paste0(s, ".gz"))
      }
      barcode.loc <- addgz(s = barcode.loc)
      matrix.loc <- addgz(s = matrix.loc)
    }
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
    }
    if (!pre_ver_3 && !file.exists(features.loc)) {
      stop("Gene name or features file missing. Expecting ", 
           basename(path = features.loc))
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing. Expecting ", 
           basename(path = matrix.loc))
    }
    data <- readMM(file = matrix.loc)
    if (has_dt) {
      cell.barcodes <- as.data.frame(data.table::fread(barcode.loc, 
                                                       header = FALSE))
    }
    else {
      cell.barcodes <- read.table(file = barcode.loc, header = FALSE, 
                                  sep = "\t", row.names = NULL)
    }
    if (ncol(x = cell.barcodes) > 1) {
      cell.names <- cell.barcodes[, cell.column]
    }
    else {
      cell.names <- readLines(con = barcode.loc)
    }
    if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
      cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, 
                                                          FUN = ExtractField, field = 1, delim = "-")))
    }
    if (is.null(x = names(x = data.dir))) {
      if (length(x = data.dir) < 2) {
        colnames(x = data) <- cell.names
      }
      else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    }
    else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], 
                                   "_", cell.names)
    }
    if (has_dt) {
      feature.names <- as.data.frame(data.table::fread(ifelse(test = pre_ver_3, 
                                                              yes = gene.loc, no = features.loc), header = FALSE))
    }
    else {
      feature.names <- read.delim(file = ifelse(test = pre_ver_3, 
                                                yes = gene.loc, no = features.loc), header = FALSE, 
                                  stringsAsFactors = FALSE)
    }
    if (any(is.na(x = feature.names[, gene.column]))) {
      warning("Some features names are NA. Replacing NA names with ID from the opposite column requested", 
              call. = FALSE, immediate. = TRUE)
      na.features <- which(x = is.na(x = feature.names[, 
                                                       gene.column]))
      replacement.column <- ifelse(test = gene.column == 
                                     2, yes = 1, no = 2)
      feature.names[na.features, gene.column] <- feature.names[na.features, 
                                                               replacement.column]
    }
    if (unique.features) {
      fcols = ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0("gene.column was set to ", gene.column, 
                    " but feature.tsv.gz (or genes.tsv) only has ", 
                    fcols, " columns.", " Try setting the gene.column argument to a value <= to ", 
                    fcols, "."))
      }
      rownames(x = data) <- make.unique(names = feature.names[, 
                                                              gene.column])
    }
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) == 
          0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) {
        lvls <- c(expr_name, lvls[-which(x = lvls == 
                                           expr_name)])
      }
      data <- lapply(X = lvls, FUN = function(l) {
        return(data[data_types == l, , drop = FALSE])
      })
      names(x = data) <- lvls
    }
    else {
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, 
                                               FUN = `[[`, j))
    list_of_data[[j]] <- as.sparse(x = list_of_data[[j]])
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  }
  else {
    return(list_of_data)
  }
}


samplecombine = function(thisdir){
  library(Seurat)
  library(Matrix)
  mtxs = dir(thisdir)[grepl("mtx.gz",dir(thisdir))]
  samplelist = as.vector(sapply(mtxs,function(x){paste(strsplit(x,"_")[[1]][c(1,2)],collapse = "_")}))
  alldata = NULL
  for (sample in samplelist){
    counts1 <- myRead10X(data.dir = thisdir, sample = sample)
    data0 <- CreateSeuratObject(counts = counts1,
                                project = sample)
    data0 <- PercentageFeatureSet(data0, pattern = "^MT-", col.name = "percent.mt")
    alldata = c(alldata, list(data0))
  }
  combine = merge(x = alldata[[1]], y = unlist(alldata)[2:length(alldata)],
                  add.cell.ids = samplelist)
  combine
}

GSE221561 = samplecombine("GSE221561")

#GSE203115
GSE203115 = samplecombine("GSE203115")

#GSE188900
samplecombine2 = function(thisdir){
  # with seperate dir for each sample
  library(Seurat)
  library(Matrix)
  
  samplelist = NULL
  for (x in dir(thisdir)){
    if(file_test("-d",paste(thisdir,x,sep="/"))){
      samplelist = c(samplelist,x)
    }
  }
  
  alldata = NULL
  for (sample in samplelist){
    counts1 <- Read10X(data.dir = paste(thisdir,sample,sep="/"))
    data0 <- CreateSeuratObject(counts = counts1,
                                project = paste(thisdir,sample,sep="_"))
    data0 <- PercentageFeatureSet(data0, pattern = "^MT-", col.name = "percent.mt")
    alldata = c(alldata, list(data0))
  }
  combine = merge(x = alldata[[1]], y = unlist(alldata)[2:length(alldata)],
                  add.cell.ids = samplelist)
  combine
}
GSE188900 = samplecombine2("GSE188900")

#GSE145370
GSE145370 = samplecombine2("GSE145370")


#GSE196756  
samplecombine = function(dir,samplelist){
  datacell = NULL
  ndiscard = NULL
  nretain = NULL
  for (sample in samplelist){
    data = filterCells(dir=dir,sample=sample)
    nretain = c(nretain,data[3][[1]][1])
    ndiscard = c(ndiscard,data[3][[1]][2])
    maxcount = data[4]
    maxgene = data[5]
    medianmt = data[6]
  }
  
}

GSE196756 = samplecombine("GSE196756")



##merge cells from cohorts

load("GSE145370.Rdata")
load("GSE188900.Rdata")
load("GSE160269.Rdata")
load("GSE196756.Rdata")
load("GSE199654.Rdata")
load("GSE203115.Rdata")
load("GSE221561.Rdata")

eso.merge = merge(x=GSE196756,
                  y=c(GSE199654,GSE221561,GSE203115,GSE188900,GSE160269,GSE145370,Karol2021,GSE222078,GSE251990),
                  add.cell.ids = c("GSE196756","GSE199654","GSE221561","GSE203115","GSE188900","GSE160269",
                                   "GSE145370","Karol2021","GSE222078","GSE251990"))
eso.merge <- PercentageFeatureSet(eso.merge, pattern = "^MT-", col.name = "percent.mt")
eso.merge@meta.data$dataset = sapply(colnames(eso.merge),function(x){paste(strsplit(x,"_")[[1]][1],collapse="-")})
eso.merge <- PercentageFeatureSet(eso.merge, pattern = "^RP[SL]",col.name = "percent.rb")

eso.merge = CreateSeuratObject(counts = GetAssayData(eso.merge),
                               meta.data = eso.merge@meta.data,
                               min.cells = 3,
                               min.features = 200)


eso.merge = subset(x = eso.merge, subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt <20 )

eso.merge.run = eso.merge
table(eso.merge$nFeature_RNA > 200 & eso.merge$nCount_RNA > 500 & eso.merge$percent.mt <20)
##Seurat clustering without batch removal
#normalization
eso.merge.run <- NormalizeData(eso.merge.run)
#high variable genes
eso.merge.run <- FindVariableFeatures(eso.merge.run)
#scale
eso.merge.run <- ScaleData(eso.merge.run)

# These are now standard steps in the Seurat workflow for visualization and clustering
eso.merge.run <- RunPCA(eso.merge.run)  
eso.merge.run <- RunUMAP(eso.merge.run, dims = 1:30)

eso.merge.run <- FindNeighbors(eso.merge.run, dims = 1:30)
eso.merge.run <- FindClusters(eso.merge.run, resolution=0.1)




###cluster umap, batch factor
library(ggplot2)
p1 <- DimPlot(eso.merge.run, reduction = "umap",group.by = "seurat_clusters", pt.size=0.2,label = T,repel = TRUE)
p1$layers[[1]]$mapping$alpha <- 0.6
p1 <- p1 + scale_alpha_continuous(range = 0.6, guide = "none")
ggplot2::ggsave("results/clusters_res0.1_alpha0.6.pdf", p1, width=12 ,height=12)

eso.merge.run@meta.data$dataset = sapply(colnames(eso.merge.run),function(x){paste(strsplit(x,"_")[[1]][1],collapse="-")})
p3 <- DimPlot(eso.merge.run, reduction = "umap",group.by = "dataset", pt.size=0.2, label = F,repel = TRUE)
ggplot2::ggsave("results/clusters_res0.1_alpha0.6_dataset.pdf", p3, width=12 ,height=12)

##fastMNN batch removal
eso.merge.run <- NormalizeData(eso.merge.run)
eso.merge.run <- FindVariableFeatures(eso.merge.run)
eso.merge.run@meta.data$dataset = sapply(colnames(eso.merge.run),
                                         function(x){paste(strsplit(x,"_")[[1]][1],collapse="-")})
library(SeuratWrappers)
eso.merge.run <- RunFastMNN(object.list = SplitObject(eso.merge.run, split.by = 'dataset'))
eso.merge.run <- RunUMAP(eso.merge.run, reduction="mnn", dims = 1:30)
eso.merge.run <- FindNeighbors(eso.merge.run,  reduction="mnn", dims = 1:30)
for (res in seq(0.1,1,by=0.1)){
  eso.merge.run <- FindClusters(eso.merge.run, resolution=res)
}
library(clustree);library(patchwork)
p1 <- clustree(eso.merge.run, prefix = 'RNA_snn_res.') + coord_flip()
ggplot2::ggsave("results/clustersfastMNN_clustree.pdf",p1,width=12,height=8)

p1 <- DimPlot(eso.merge.run, reduction = "umap",group.by = "seurat_clusters", pt.size=0.1,label = T,repel = TRUE)
p1$layers[[1]]$mapping$alpha <- 0.3
p1 <- p1 + scale_alpha_continuous(range = 0.3, guide = "none")
ggplot2::ggsave("results/clustersfastMNN_res0.1_alpha0.6.pdf", p1, width=12 ,height=12)
p3 <- DimPlot(eso.merge.run, reduction = 'umap', group.by = 'dataset', pt.size=0.2, label = F,repel = TRUE, raster=T)
#p3$layers[[1]]$mapping$alpha <- 0.3
#p3 <- p3 + scale_alpha_continuous(range = 0.3, guide = "none")
ggplot2::ggsave("clustersfastMNN_res0.1_alpha0.6_dataset.pdf", p3, width=12 ,height=12)
p4 <- DimPlot(eso.merge.run, reduction = 'umap', split.by = 'dataset', group.by =  ncol=3, pt.size=0.2, label = F,repel = TRUE, raster= T)
p4$layers[[1]]$mapping$alpha <- 0.3
p4 <- p4 + scale_alpha_continuous(range = 0.3, guide = "none")
ggplot2::ggsave("clustersfastMNN_res0.1_alpha0.6_splitdataset.pdf", p4, width=12 ,height=16)

eso.cells = eso.merge.run@meta.data[!duplicated(eso.merge.run@meta.data$orig.ident),c("dataset","orig.ident")]
write.table(eso.cells,"results/cells_unique.csv",sep=",",quote=F,row.names=F)


#marker genes if each cluster for manually determine the cell types
Idents(eso.merge.run) <- 'seurat_clusters'
allmarkers <- FindAllMarkers(eso.merge.run, logfc.threshold = 0.5, min.pct = 0.25, only.pos = TRUE)
write.table(allmarkers,"/data_group/cunyupeng/songjing/project/01.xiaohua.tumors/geodata/esophagus/results/eso_30pc_allmarkers.csv",sep=",",row.names=F,quote=F)



####cell types determined

eso.merge.run$celltype_top = as.character(eso.merge.run$seurat_clusters)
eso.merge.run$celltype_top[eso.merge.run$celltype_top %in% c(0,8,16,22,26)] = "Epithelial cells"
eso.merge.run$celltype_top[eso.merge.run$celltype_top %in% c(1,19)] = "Myeloid"
eso.merge.run$celltype_top[eso.merge.run$celltype_top %in% c(2)] = "Basal cells"
eso.merge.run$celltype_top[eso.merge.run$celltype_top %in% c(12)] = "Endothelial cells"
eso.merge.run$celltype_top[eso.merge.run$celltype_top %in% c(7,10,20)] = "Fibroblasts"
eso.merge.run$celltype_top[eso.merge.run$celltype_top %in% c(17)] = "Mast cells"
eso.merge.run$celltype_top[eso.merge.run$celltype_top %in% c(3)] = "B cells"
eso.merge.run$celltype_top[eso.merge.run$celltype_top %in% c(15,18,24)] = "Plasma cells"
eso.merge.run$celltype_top[eso.merge.run$celltype_top %in% c(4,5,6,9,11,13,14,21,23,25)] = "TNK cells"
eso.merge.run$celltype_top[eso.merge.run$dataset %in% c("GSE199654")] = "Epithelial cells"

colors = c("Epithelial cells"='#CE0071',"Myeloid"='#FF3900',"Basal cells"='#00B454',
           "Endothelial cells"='#9BED00',"Fibroblasts"='#6699ff',"Mast cells"='#ff9900',
           "B cells"='#66cccc','Plasma cells'='#FFDA00','TNK cells'='#ff9999')
p1 <- DimPlot(eso.merge.run, reduction = "umap",group.by = "celltype_top", pt.size=0.1,label = T,repel = TRUE,cols = colors) + NoLegend()
p1$layers[[1]]$mapping$alpha <- 0.3
p1 <- p1 + scale_alpha_continuous(range = 0.3, guide = "none")
ggplot2::ggsave("results/celltype_toplevel.pdf", p1, width=12 ,height=12)




### Potenital doublets removal (first round)
### For second round (in cell subset annotation, e.g., CD4 T, CD8 T), you need manually determine the doublet clusters and remove them
## BoubletFinder run by dataset
library(Seurat)
library(dplyr)
###scDblFinder
load("esophagus.merge.obj.Rdata")
sce = as.SingleCellExperiment(eso.merge)
rm(esp.merge)
sce@metadata = sce@colData@listData
library(scDblFinder)
sce <- scDblFinder(sce, dbr=0.076)
write.table(sce@colData@listData,"results/scDblFinder.res.dbr0.76.txt",sep="\t",quote=F)

#######in subsequent sub celltype annotation
dbl = read.table("results/scDblFinder.res.dbr0.76.txt",sep="\t",header=T)
eso.merge.run.Myeloid$scDblFinder.class = dbl$scDblFinder.class[match(colnames(eso.merge.run.Myeloid),dbl$cellbarcode)]
Idents(eso.merge.run.Myeloid) = "scDblFinder.class"
eso.merge.run.Myeloid0 = subset(eso.merge.run.Myeloid, idents="singlet")
  




##cell propotion in groups ESCC
eso.merge.run$tmp = eso.merge.run$SpecimenGroup
eso.merge.run$tmp[eso.merge.run$tmp %in% c("Adjacent esophagus",
                                           "ESCC I","ESCC II","ESCC III","ESCC IVA") == FALSE] = "other"
eso.merge.run$tmp = factor(eso.merge.run$tmp,levels=c("Adjacent esophagus",
                                                      "ESCC I","ESCC II","ESCC III","ESCC IVA"))
eso.merge.run$tmp = as.character(eso.merge.run$tmp)

Idents(eso.merge.run) = "celltype_top"
p4 <- DimPlot(eso.merge.run, reduction = 'umap',  ncol=3, pt.size=0.1, label = F,repel = TRUE)#split.by = 'tmp',
p4$layers[[1]]$mapping$alpha <- 0.3
p4 <- p4 + scale_alpha_continuous(range = 0.3, guide = "none")
ggplot2::ggsave("tmp.pdf", p4, width=9 ,height=8)

plotC <- reshape2::melt(table(eso.merge.run$TreatmentResponse, eso.merge.run$celltype_top))
colnames(plotC) <- c("Sample", "CellType","Number")
plotC$CellType = factor(plotC$CellType,levels=c("Epithelial cells","Basal cells","Fibroblasts","Endothelial cells",
                                                'TNK cells',"Myeloid","B cells",'Plasma cells',
                                                "Mast cells"))
plotC$Prop = apply(plotC,1,function(x){ paste0(round(as.numeric(x[3]) / sum(plotC$Number[plotC$Sample == x[1]]),4)*100,"%")})
pC2 = ggplot(data = plotC, aes(x = Sample, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill")+
  scale_fill_manual(values=colors) + theme_bw()+
  #geom_text(aes(label = Prop), stat="identity",colour = "black",position = position_fill(vjust = 0.5)) +
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
  scale_y_continuous(labels = c(0,0.25,0.5,0.75,1))+ 
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.8))
pC2
ggplot2::ggsave("tmp.pdf",pC2,width=6,height=6)



##marker expression
marker_cell = c(
  "EPCAM","KRT14","TP63", #EPI
  "KRT4","KRT13","GJB2", #Basal
  "DCN","ACTA2","FBLN1", #Fibro
  "PECAM1","ENG","PLVAP", #Endo
  "CD3D","CD3E","CD8A","NKG7", #TNK
  "CD68","FCGR3A","CD14", #Myeloid
  "CD19","MS4A1","BANK1", #B
  "MZB1","JCHAIN","DERL3", #Plasma
  "FCER1A","TPSB2","CPA3" #mast
  
)
p1= DotPlot(object = eso.merge, 
            features=unique(marker_cell),  col.min = 0, cluster.idents = FALSE,
            assay = "RNA") + theme(axis.text.x = element_text(angle = 90, 
                                                              vjust = 0.5, hjust=0.5))
ggplot2::ggsave(paste0("tmp.pdf"), p1, width=10 ,height=6)






