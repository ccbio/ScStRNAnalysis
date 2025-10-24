
###single-cell RNA-seq
###HMGCR MIF CholesterolScore
Idents(eso.merge.run) = "SpecimenGroup"
eso.merge.run = subset(eso.merge.run, idents = c("ESCC I","ESCC II","ESCC III","ESCC IVA"))
load("results/epi/ESCCepi.obj.Rdata")
eso.merge.run@meta.data$celltype_top2 = eso.merge.run@meta.data$celltype_top
eso.merge.run@meta.data$celltype_top2[colnames(eso.merge.run) %in% colnames(ESCCepi)] = "Malignant"
eso.merge.run@meta.data$celltype_top2[colnames(eso.merge.run) %in% colnames(ESCCepi)[ESCCepi$CH.cell == "TRUE"]] = "CH.cell"
genelist = c("MSMO1","DHCR7","MVD","INSIG1","IDI1","CYP51A1","HMGCR","HMGCS1","HSD17B7","TM7SF2","FDFT1","FDPS","CES1","MVK",
             "SC5D","LSS","DHCR24","G6PD","SREBF2","ACLY","SQLE","ACAT2")
eso.merge.run <- AddModuleScore(
  object = eso.merge.run,
  features = list(cholesterol_pathway = genelist),  
  name = "CholesterolScore",                        # generate CholesterolScore1
  ctrl = 100,                                       
  seed = 42                                         
)

ggplot2::ggsave("CholesterolScore_umap.pdf", p5, width=5, height=5)

load("ESCCepi.obj.Rdata")
p <- DotPlot(eso.merge.run, features=c("MIF","CholesterolScore1","HMGCR"), group.by = 'celltype_top2') + RotatedAxis() +
  scale_color_gradientn(
    colors = c("grey", "#48DD00", "#FFC600", "red"),
    values = scales::rescale(c(0, 0.33, 0.67, 1))  
  )

ggplot2::ggsave("CholesterolScore_Dotplot.pdf",height = 4,width=6)


Idents(eso.merge.run) = "celltype_top2"
p5 = FeaturePlot(eso.merge.run, features = c("MIF","CholesterolScore1","HMGCR"), pt.size = 0.1,
                 reduction = "umap", raster=T, cols = c("lightgrey", "darkgreen","red"), ncol=3) 
p5$layers[[1]]$mapping$alpha <- 0.3
p5 <- p5 + scale_alpha_continuous(range = 0.3, guide = "none")
ggplot2::ggsave("MIF_HMGCR_gene_exp.pdf", p5, width=30, height=10)




#####spatial RNA-seq
#######CHscore, HMGCR, MIF expression proportion
load("RCTD_eso.sp.Rdata")
genelist = c("MSMO1","DHCR7","MVD","INSIG1","IDI1","CYP51A1","HMGCR","HMGCS1","HSD17B7","TM7SF2","FDFT1","FDPS","CES1","MVK",
             "SC5D","LSS","DHCR24","G6PD","SREBF2","ACLY","SQLE","ACAT2")
eso.sp <- AddModuleScore(
  object = eso.sp,
  features = list(cholesterol_pathway = genelist), 
  name = "CholesterolScore",                        # generate CholesterolScore1
  ctrl = 100,                                       
  seed = 10                                         
)

library(ggplot2)
library(Seurat)

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

### GC.B.dividing = GC.B dark;  GC.B = GC.B light; CH.cell = Cholesterol biosynthesis tumor cell
cellinfer$celltype.comb = cellinfer$celltype
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("CD4.Tex","CD4.CTL","CD4.T.Naive","Treg","Th22","Th17")] = "T helper"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("NKT","CD8.T.cytotoxic","CD8.Tcm.CD69.","CD8.T.Naive","gdT","CD8.Tex","dnT","CD8.Tem","gdT")] = "T cytotoxic"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("Activated.B","Naive.B","Memory.B","GC.B.dividing","GC.B")] = "B"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("Pluripotent.Endo","Microvascular.Endo","Hypoxia.Inflammatory.Endo","Lymphatic.Endo")] = "Endo"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("Pericytes","VascularCAFs","MatrixCAFs","AnrigenPresentCAF","MixedCAFs.vCAF.tCAF.","MixedCAFs.iCAF.mCAF.tCAF.","MixedCAFs.iCAF.apCAF.","InflammatoryCAFs")] = "Fibroblasts"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("M2.Macrophage","CD16.Mono","M1.Macrophage","CD14.Mono","Megakaryocyte")] = "Mono"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("cDC1","cDC2","pDC","tDC")] = "DC"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("Basal.Stem..BS.","Basal.Keratinocytes..BK.","Differentiated.Keratinocytes..DK.")] = "NormalEpi"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("ESCCmalignant")] = "Malignant"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("CH.cell")] = "CHolesterol.cell"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("IgA.Plasma","Plasmablast","IgG.Plasma")] = "Plasma"
cellinfer$celltype.comb[cellinfer$celltype.comb %in% c("TRNs")] = "TRNs"

eso.sp@meta.data$celltype_comb = cellinfer$celltype.comb

eso.sp$MIF_expression <- GetAssayData(eso.sp)["MIF", ]
eso.sp$HMGCR_expression = GetAssayData(eso.sp)["HMGCR", ]


counts_matrix <- GetAssayData(eso.sp, assay = "Spatial.008um", slot = "counts")
metadata <- eso.sp@meta.data


sc_obj <- CreateSeuratObject(
  counts = counts_matrix,
  meta.data = metadata,
  project = "Converted",
  assay = "RNA"
)
sc_obj = subset(sc_obj, cells = colnames(sc_obj)[is.na(sc_obj$celltype_comb) == FALSE])
p1 <- DotPlot(sc_obj, features=c("MIF_expression"), group.by = 'celltype_comb') + RotatedAxis() +
  scale_color_gradientn(
    colors = c("grey", "#48DD00", "#FFC600", "red"),
    values = scales::rescale(c(0, 0.33, 0.67, 1))
  )
p2<- DotPlot(sc_obj, features=c("CholesterolScore1"), group.by = 'celltype_comb') + RotatedAxis() +
  scale_color_gradientn(
    colors = c("grey", "#48DD00", "#FFC600", "red"),
    values = scales::rescale(c(0, 0.33, 0.67, 1))
  )
p3<- DotPlot(sc_obj, features=c("HMGCR_expression"), group.by = 'celltype_comb') + RotatedAxis() +
  scale_color_gradientn(
    colors = c("grey", "#48DD00", "#FFC600", "red"),
    values = scales::rescale(c(0, 0.33, 0.67, 1))  
  )
ggplot2::ggsave("CholesterolScore_spatial_Dotplot.pdf",p1+p2+p3, height = 4,width=12)
