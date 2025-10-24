# packages required for Visium HD
install.packages("hdf5r")
install.packages("arrow")

library(hdf5r)
library(arrow)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
setwd("Visium HD sequencing/analysis/")


Load10X_Spatial = function (data.dir, filename = "filtered_feature_bc_matrix.h5", 
                            assay = "Spatial", slice = "slice1", bin.size = NULL, filter.matrix = TRUE, 
                            to.upper = FALSE, image = NULL, ...) 
{
  if (length(x = data.dir) > 1) {
    data.dir <- data.dir[1]
    warning(paste0("`data.dir` expects a single value but recieved multiple - ", 
                   "continuing using the first: '", data.dir, "'."), 
            immediate. = TRUE, )
  }
  if (!file.exists(data.dir)) {
    stop(paste0("No such file or directory: ", "'", data.dir, 
                "'"))
  }
  if (is.null(bin.size) & file.exists(paste0(data.dir, "/binned_outputs"))) {
    bin.size <- c(16, 8)
  }
  if (!is.null(bin.size)) {
    bin.size.pretty <- paste0(sprintf("%03d", bin.size), 
                              "um")
    data.dirs <- paste0(data.dir, "/binned_outputs/", "square_", 
                        bin.size.pretty)
    assay.names <- paste0(assay, ".", bin.size.pretty)
    slice.names <- paste0(slice, ".", bin.size.pretty)
  }
  else {
    data.dirs <- data.dir
    assay.names <- assay
    slice.names <- slice
  }
  counts.paths <- lapply(data.dirs, file.path, filename)
  counts.list <- lapply(counts.paths, Read10X_h5, ...)
  if (to.upper) {
    rownames(counts) <- lapply(rownames(counts), toupper)
  }
  if (is.null(image)) {
    image.list <- mapply(Read10X_Image, file.path(data.dirs, 
                                                  "spatial"), assay = assay.names, slice = slice.names, 
                         MoreArgs = list(filter.matrix = filter.matrix))
  }
  else {
    image.list <- c(image)
  }
  if (length(image.list) != length(counts.list)) {
    stop(paste0("The number of images does not match the number of counts matrices. ", 
                "Ensure each spatial dataset has a corresponding image."))
  }
  object.list <- mapply(CreateSeuratObject, counts.list, assay = assay.names)
  object.list <- mapply(function(.object, .image, .assay, .slice) {
    .image <- .image[Cells(.object)]
    .object[[.slice]] <- .image
    return(.object)
  }, object.list, image.list, assay.names, slice.names)
  object <- merge(object.list[[1]], y = object.list[-1])
  return(object)
}

read.spatialHD.10X = function(thisdir){
  library(Seurat)
  library(Matrix)
  dirs = NULL
  for (x in dir(thisdir)){
    if(file_test("-d",paste(thisdir,x,sep="/"))){
      dirs = c(dirs,x)
    }
  }
  dirs = dirs[grep("run",dirs)]
  alldata = NULL
  for (sample in 1:length(dirs)){
    show(dirs[sample])
    spdata1 <- Load10X_Spatial(data.dir = paste(thisdir,dirs[sample],"outs",sep="/"), 
                               slice = dirs[sample], bin.size = c(8))
    
    alldata = c(alldata,list(spdata1))
  }
  names(alldata) = dirs
  alldata
}

object = read.spatialHD.10X(getwd())

eso.sp = merge(x=object[[1]], y=c(object[[2]],object[[3]],object[[4]]),
               add.cell.ids = c("D202301729","X202300185","X202316042","X202319854"))
eso.sp[["Spatial.008um"]] <- as(object = eso.sp[["Spatial.008um"]], Class = "Assay")

filter = subset(eso.sp, cells = colnames(eso.sp)[eso.sp@meta.data$nFeature_Spatial.008um >=10])
eso.sp = filter

getCoord = function(eso.sp, sam){
  for (sam in substr(names(eso.sp@images),1,10)){
    polygons = read.csv(paste(getwd(),paste0(sam,"run"),"outs",paste(sam,"_geom.csv",sep=""),sep="/"),header=T,stringsAsFactors = F)
    tmp = lapply(polygons$centroid, function(x){ stringr::str_extract(strsplit(x," ")[[1]], "\\d+")[-1] })
    centroid_coord = data.frame(x = as.numeric(data.frame(tmp)[seq(1,length(tmp[[1]]),2),]),
                                y = as.numeric(data.frame(tmp)[seq(2,length(tmp[[1]]),2),]) )
    geometry_coord = data.frame( lapply(polygons$geometry, function(x){ stringr::str_extract(strsplit(x," ")[[1]], "\\d+")[-1] }) )
    
  }
  coords_polygons = list(geometry_coord=geometry_coord, centroid_coord=centroid_coord)
  
  allpolys = list()
  for (i in 1:length(coords_polygons$geometry_coord)){
    show(i)
    xy = coords_polygons$geometry_coord[[i]]
    x = xy[seq(1,length(xy),2)]
    y = xy[seq(2,length(xy),2)]
    center = c(coords_polygons$centroid_coord$x[i], coords_polygons$centroid_coord$y[i])
    x_enlarge = as.numeric(x) + 25
    less = which(x < center[1])
    x_enlarge[less] = x_enlarge[less] - 50
    y_enlarge = as.numeric(y) + 25
    less = which(y < center[2])
    y_enlarge[less] = y_enlarge[less] - 50
    poly.pnts = as.matrix(data.frame(x=x_enlarge, y=y_enlarge))
    allpolys = c(allpolys, list(poly.pnts))
  }
  allpolys
}

getBCsOfTissueRegion = function(eso.sp, sam){
  allpolys = getCoord(eso.sp, sam)
  library(DescTools)
  pnts = eso.sp@images[grep(sam,names(eso.sp@images))][[1]]@boundaries$centroids@coords
  cells = eso.sp@images[grep(sam,names(eso.sp@images))][[1]]@boundaries$centroids@cells
  
  allhits = NULL
  for (x in allpolys){
    ptinpoly = sp::point.in.polygon(pnts[,1],pnts[,2],x[,1],x[,2])
    ind = which(ptinpoly >= 1)
    if(length(ptinpoly[ind]) > 0){
      allhits = rbind(allhits, cells[ind])
    }
  }
  unique(allhits)
}


hits.X202319854 = getBCsOfTissueRegion(eso.sp, sam)

eso.sp <- PercentageFeatureSet(eso.sp, pattern = "^MT-", col.name = "percent.mt")
eso.sp <- PercentageFeatureSet(eso.sp, pattern = "^RP[SL]",col.name = "percent.rb")







####RCTD
library(spacexr)
refdir = "results"
load(paste(refdir,"eso.merge.run.Rdata",sep="/"))
celltype_for_ref = factor(eso.merge.run@meta.data$celltype_sub)
names(celltype_for_ref) = colnames(eso.merge.run)
counts_for_ref = GetAssayData(eso.merge.run,assay = "RNA",slot = "counts")

reference <- Reference(counts_for_ref, celltype_for_ref)
save(reference,file="RCTD_reference.Rdata")
####run RCTD infering cell weights

save(eso.sp,file="RCTD_eso.sp.Rdata")
for (sam in names(eso.sp@images)){
  show(sam)
  coords_for_spacexr = data.frame(row=eso.sp@images[grep(sam,names(eso.sp@images))][[1]]@boundaries$centroids@coords[,1],
                                  col=eso.sp@images[grep(sam,names(eso.sp@images))][[1]]@boundaries$centroids@coords[,2],
                                  row.names = eso.sp@images[grep(sam,names(eso.sp@images))][[1]]@boundaries$centroids@cells)
  counts_for_spacexr = GetAssayData(eso.sp,assay = "Spatial.008um",slot = "counts")
  spaceRNA = SpatialRNA(coords_for_spacexr, counts_for_spacexr)
  
  myRCTD <- create.RCTD(spaceRNA, reference)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
  write.table(myRCTD@results$weights,paste0("res_spacexr_all_",sam,".tsv"),sep="\t",quote=F)
  
}

load("RCTD_eso.sp.Rdata")
eso.sp$orig.ident = sapply(rownames(eso.sp@meta.data),function(x){strsplit(x,"_")[[1]][1]}) 
cellinfer = data.frame()
for (i in dir("spacexr")){
  if(grepl("res_spacexr_all",i) == FALSE){ next }
  show(i)
  f = read.table(paste0("spacexr",i),sep="\t",header=T,check.names=F, stringsAsFactors = FALSE)
  cellinfer = rbind(cellinfer,round(f,2))
}
cellinfer = cellinfer[match(colnames(eso.sp),rownames(cellinfer)),]

cellinfer$celltype = apply(cellinfer,1,function(x){
  x = na.omit(x)
  names = colnames(cellinfer)
  ind = which.max(x)
  val = max(x)
  if(length(ind) == 0){
    NA
  }else{
    if(is.numeric(val) & val >= 0.1){
      names[ind]
    }else{
      NA
    }
  }
  
})
table(eso.sp$orig.ident,eso.sp$celltype_comb)
apply(table(eso.sp$orig.ident,eso.sp$celltype_comb),1,function(x){x/sum(x)})
apply(table(eso.sp$orig.ident,eso.sp$celltype_comb)[,c(2,6,8)],1,function(x){x/sum(x)})
apply(table(eso.sp$orig.ident,eso.sp$celltype_comb)[,c(1,3:5,7,9:13)],1,function(x){x/sum(x)})
top = rbind(apply(table(eso.sp$orig.ident,eso.sp$celltype_comb),1,function(x){x/sum(x)}),
            apply(table(eso.sp$orig.ident,eso.sp$celltype_comb)[,c(2,6,8)],1,function(x){x/sum(x)}))


# T,B,CAF,Endo,Myeloid,Mast,Plasma,Epi,Malignant,CHolesterol.cell
cellinfer$celltype.comb = cellinfer$celltype
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("CD4.Tex","CD4.CTL","CD4.T.Naive","Treg","Th22","Th17")] = "T helper"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("NKT","CD8.T.cytotoxic","CD8.Tcm.CD69.","CD8.T.Naive","gdT","CD8.Tex","dnT","CD8.Tem","gdT")] = "T cytotoxic"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("Activated.B","Naive.B","Memory.B","GC.B.dark","GC.B.light")] = "B"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("Pluripotent.Endo","Microvascular.Endo","Hypoxia.Inflammatory.Endo","Lymphatic.Endo")] = "Endo"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("Pericytes","VascularCAFs","MatrixCAFs","AnrigenPresentCAF","MixedCAFs.vCAF.tCAF.","MixedCAFs.iCAF.mCAF.tCAF.","MixedCAFs.iCAF.apCAF.","InflammatoryCAFs")] = "Fibroblasts"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("M2.Macrophage","CD16.Mono","M1.Macrophage","CD14.Mono","Megakaryocyte")] = "Mono"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("cDC1","cDC2","pDC","tDC")] = "DC"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("Basal.Stem..BS.","Basal.Keratinocytes..BK.","Differentiated.Keratinocytes..DK.")] = "NormalEpi"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("ESCCmalignant")] = "Malignant"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("CHolesterol.cell")] = "CHolesterol.cell"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("IgA.Plasma","Plasmablast","IgG.Plasma")] = "Plasma"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("TRNs")] = "TRNs"

eso.sp@meta.data$celltype_comb = cellinfer$celltype.comb
#rev(paletteer::paletteer_d("ggthemes::Classic_20")[seq(1,20,by=2)])
col_assign <- c("B"="#17BECFFF","Plasma"="#9EDAE5FF",
                "T helper"="#FF7F0EFF","T cytotoxic"="#FFBB78FF",
                "Endo"="#BCBD22FF","Fibroblasts"="#E377C2FF",
                "Malignant"="#8C564BFF","NormalEpi"="#2CA02CFF",
                "Mono"="#9467BDFF","TRNs" = "#AEC7E8FF","DC" = "#1F77B4FF",
                "CHolesterol.cell"="#D62728FF"#,NA="#7F7F7FFF",
)

library(ggplot2)
library(Seurat)
eso.sp.noNA = subset(eso.sp, cells = colnames(eso.sp)[is.na(eso.sp$celltype_comb) == FALSE])
p2 = SpatialDimPlot(eso.sp.noNA, group.by = "celltype_comb", cols = col_assign, 
                    image.alpha = 0.8, pt.size.factor = 3,
                    shape = 22, ncol =2) 

ggplot2::ggsave(paste0("SpatialDimplot_celltype_comb",".pdf"),p2, height=20,width=20)

cellinfer$celltype.comb = cellinfer$celltype
eso.sp@meta.data$celltype_comb = cellinfer$celltype.comb


eso.sp@meta.data$celltype_comb2 = ifelse(eso.sp@meta.data$celltype_comb %in% c("CHolesterol.cell","GC.B.dark","NA") == TRUE,eso.sp@meta.data$celltype_comb,"other" )
col_assign <- c("GC.B.dividing"="#176e5FfF","other"="grey80",
                "CHolesterol.cell"="#D62728FF"
)
library(ggplot2)
library(Seurat)
eso.sp.noNA = subset(eso.sp, cells = colnames(eso.sp)[is.na(eso.sp$celltype_comb) == FALSE])
p2 = SpatialDimPlot(eso.sp.noNA, group.by = "celltype_comb2", cols = col_assign, 
                    image.alpha = 0.8, pt.size.factor = 3,
                    shape = 22, ncol =2)

ggplot2::ggsave(paste0("SpatialDimplot_celltype_sub2",".pdf"),p2, height=20,width=20)


####TLS score
chemokine12 = c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19","CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13")
eso.sp.noNA = AddModuleScore(eso.sp.noNA, features = list(chemokine12),  name = "chemokine12")

library(ggplot2)
library(Seurat)
p2 = SpatialFeaturePlot(eso.sp.noNA, features = "chemokine12",  
                        image.alpha = 0.8, pt.size.factor = 3,
                        shape = 22, ncol =2)

ggplot2::ggsave(paste0("SpatialDimplot_TLSscore",".pdf"),p2, height=20,width=20)


####cxcr4-mif CD74-MIF
cellinfer$celltype.comb = cellinfer$celltype
eso.sp@meta.data$celltype_comb = cellinfer$celltype.comb
data = GetAssayData(eso.sp,slot="data")
lrexp = data.frame(t(data.frame(data[rownames(data) %in% c("MIF","CXCR4","CD74"),])))
lrexp$MIF_CXCR4 = ifelse(lrexp$CXCR4>0 & lrexp$MIF >0, 1, 0)
lrexp$MIF_CD74 = ifelse(lrexp$CD74>0 & lrexp$MIF >0, 1, 0)
lrexp$cellint = ifelse(eso.sp@meta.data$celltype_comb %in% c("CHolesterol.cell","GC.B.dark","NA") == TRUE,eso.sp@meta.data$celltype_comb,"other" )
table(lrexp$cellint, lrexp$MIF_CXCR4)
table(lrexp$cellint, lrexp$MIF_CD74)
lrexp$MIF_CXCR4 = ifelse(lrexp$MIF_CXCR4 == 1 & lrexp$cellint != "other", 1, 0)
lrexp$MIF_CD74 = ifelse(lrexp$MIF_CD74 == 1 & lrexp$cellint != "other", 1, 0)
lrexp$orig = eso.sp$orig.ident
eso.sp@meta.data$celltype_comb3 = lrexp$MIF_CXCR4

col_assign <- c("1"="#D62728FF","0"="grey80")
library(ggplot2)
library(Seurat)
eso.sp.noNA = subset(eso.sp, cells = colnames(eso.sp)[is.na(eso.sp$celltype_comb3) == FALSE])
p2 = SpatialDimPlot(eso.sp.noNA, group.by = "celltype_comb3", cols = col_assign, 
                    image.alpha = 0.8, pt.size.factor = 3,
                    shape = 22, ncol =2)

ggplot2::ggsave(paste0("SpatialDimplot_LR_CXCR4",".pdf"),p2, height=20,width=20)





### niche of chol. cell 
getDirectNeigborBCs = function(obj,celltype,micron=8,res=8){
  '%!in%' <- function(x,y)!('%in%'(x,y))
  micron = sqrt(2*micron^2) / res * 3 +0.2  #box(x,y) increase by 3
  celltype.list = na.omit(unique(eso.sp@meta.data$celltype_comb))
  celltype.list = celltype.list[celltype.list %!in% c("Malignant","CHolesterol.cell","NormalEpi","ESCCmalignant","Basal.Stem..BS.","Basal.Keratinocytes..BK.","Differentiated.Keratinocytes..DK.")]
  BCs.celltype = na.omit(rownames(obj@meta.data)[obj@meta.data$celltype_comb == celltype])
  res = NULL
  for (sam in names(obj@images)){
    show(sam)
    scalef = obj@images[grep(sam,names(obj@images))][[1]]@scale.factors$spot
    coords = data.frame(row=floor(obj@images[grep(sam,names(obj@images))][[1]]@boundaries$centroids@coords[,1] * scalef),
                        col=floor(obj@images[grep(sam,names(obj@images))][[1]]@boundaries$centroids@coords[,2] * scalef),
                        row.names = obj@images[grep(sam,names(obj@images))][[1]]@boundaries$centroids@cells)
    celltype.in = coords[rownames(coords) %in% BCs.celltype,]
    celltype.other = coords[rownames(coords) %!in% BCs.celltype ,]
    
    allneighbors = NULL; celltype.in$neighbors = NA
    for(j in 1:length(celltype.in$row)){
      Distance = sqrt((celltype.other$row - celltype.in$row[j])^2 + (celltype.other$col - celltype.in$col[j])^2)
      neighbor.BCs = rownames(celltype.other)[Distance < micron]
      #celltype.in$neighbors[j] = paste(neighbor.BCs,collapse = ",")
      allneighbors = c(allneighbors, neighbor.BCs)
    }
    
    neighbor.stat = table(obj@meta.data$celltype_comb[ match(allneighbors,rownames(obj@meta.data)) ])
    neighbor.df = data.frame(
      label = celltype.list,
      count = rep(0,length(celltype.list)),
      stringsAsFactors = F
    )
    neighbor.df$count = as.vector(neighbor.stat[match(celltype.list,names(neighbor.stat))])
    neighbor.df$count[is.na(neighbor.df$count)] = 0
    res = cbind(res, neighbor.df$count)
  }
  colnames(res) = names(obj@images)
  rownames(res) = celltype.list
  res.prop = apply(res,2,function(x){x/sum(x)})
  list(count=data.frame(res),prop=data.frame(res.prop))
}
eso.sp@meta.data$celltype_comb = cellinfer$celltype.comb
res = getDirectNeigborBCs(obj = eso.sp, celltype = "CHolesterol.cell",micron=8,res=8)

eso.sp@meta.data$celltype_comb = cellinfer$celltype
res.sep = getDirectNeigborBCs(obj = eso.sp, celltype = "CHolesterol.cell",micron=8,res=8)

res.sep$prop$Resistant2vs1 = res.sep$prop[,colnames(res.sep$prop) == "D202301729run.008um"] - res.sep$prop[,colnames(res.sep$prop) == "X202300185run.008um"]
res.sep$prop$Response2vs1 = res.sep$prop[,colnames(res.sep$prop) == "X202319854run.008um"] - res.sep$prop[,colnames(res.sep$prop) == "X202316042run.008um"]
res.sep$prop = res.sep$prop[rev(order(res.sep$prop$Resistant2vs1,res.sep$prop$Response2vs1)),]






#########
###########spatial : show MIF+ cells and CD4/CD8 TEX


load("RCTD_eso.sp.Rdata")
eso.sp$orig.ident = sapply(rownames(eso.sp@meta.data),function(x){strsplit(x,"_")[[1]][1]}) 
cellinfer = data.frame()
for (i in dir("spacexr")){
  if(grepl("res_spacexr_all",i) == FALSE){ next }
  show(i)
  f = read.table(paste0("spacexr/",i),sep="\t",header=T,check.names=F, stringsAsFactors = FALSE)
  cellinfer = rbind(cellinfer,round(f,2))
}
cellinfer = cellinfer[match(colnames(eso.sp),rownames(cellinfer)),]

cellinfer$celltype = apply(cellinfer,1,function(x){
  x = na.omit(x)
  names = colnames(cellinfer)
  ind = which.max(x)
  val = max(x)
  if(length(ind) == 0){
    NA
  }else{
    if(is.numeric(val) & val >= 0.1){
      names[ind]
    }else{
      NA
    }
  }
  
})

cellinfer$celltype.comb = cellinfer$celltype

eso.sp@meta.data$celltype_comb = "other"
eso.sp@meta.data$celltype_comb[cellinfer$celltype.comb %in% c("CD4.Tex")] = "CD4.Tex"
eso.sp@meta.data$celltype_comb[cellinfer$celltype.comb %in% c("CD8.Tex")] = "CD8.Tex"
eso.sp@meta.data$celltype_comb[cellinfer$celltype.comb %in% c("ESCCmalignant","CH.cell")] = "Malignant"
eso.sp@meta.data$celltype_comb[cellinfer$celltype.comb %in% c("GC.B","GC.B.dividing")] = "GCB"

data = GetAssayData(eso.sp,slot="data")
lrexp = data.frame(t(data.frame(data[rownames(data) %in% c("VEGFA","MIF"),])))
eso.sp@meta.data$celltype_comb[lrexp$MIF>0 ] = "MIFpos"

library(ggplot2)
library(Seurat)
col_assign <- c("CD4.Tex"="#fff740","CD8.Tex"="#d1f","GCB" = "#176e5FfF",
                "other"="grey80","Malignant" = "black",
                "MIFpos"="red"
)
eso.sp.noNA = subset(eso.sp, cells = colnames(eso.sp)[is.na(eso.sp$celltype_comb) == FALSE])
p2 = SpatialDimPlot(eso.sp.noNA, group.by = "celltype_comb", cols = col_assign, 
                    image.alpha = 0.8, pt.size.factor = 3,
                    shape = 22, ncol =2)

ggplot2::ggsave(paste0("SpatialDimplot_celltype_Tex_MIF",".pdf"),p2, height=20,width=20)



##macrophages
eso.sp@meta.data$celltype_comb = "other"
eso.sp@meta.data$celltype_comb[cellinfer$celltype.comb %in% c("M2.Macrophage")] = "M2.Macrophage"
eso.sp@meta.data$celltype_comb[cellinfer$celltype.comb %in% c("M1.Macrophage")] = "M1.Macrophage"
eso.sp@meta.data$celltype_comb[cellinfer$celltype.comb %in% c("ESCCmalignant","CH.cell")] = "Malignant"
eso.sp@meta.data$celltype_comb[cellinfer$celltype.comb %in% c("GC.B","GC.B.dividing")] = "GCB"



library(ggplot2)
library(Seurat)
col_assign <- c("M1.Macrophage"="#fff740","M2.Macrophage"="#d1f","GCB" = "#176e5FfF",
                "other"="grey80","Malignant" = "black"
)
eso.sp.noNA = subset(eso.sp, cells = colnames(eso.sp)[is.na(eso.sp$celltype_comb) == FALSE])
p2 = SpatialDimPlot(eso.sp.noNA, group.by = "celltype_comb", cols = col_assign, 
                    image.alpha = 0.8, pt.size.factor = 3,
                    shape = 22, ncol =2)

ggplot2::ggsave(paste0("SpatialDimplot_celltype_MIFpos_macrophages",".pdf"),p2, height=20,width=20)
