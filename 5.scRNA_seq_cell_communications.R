

###cellchat
load("../esophagus.merge.obj.run.Rdata")
load("../results/epi/ESCCepi.Rdata")

eso.merge.run$celltype_sub = NA
eso.merge.run$celltype_sub[match(colnames(ESCCepi),colnames(eso.merge.run))] = "ESCCmalignant"
eso.merge.run$celltype_sub[match(colnames(ESCCepi)[ESCCepi$CHolesterol.cell=="TRUE"],colnames(eso.merge.run))] = "CHolesterol.cell"


Idents(eso.merge.run) = "celltype_sub"
eso.merge.run = subset(eso.merge.run, idents = na.omit(unique(eso.merge.run$celltype_sub)) )
Idents(eso.merge.run) = "SpecimenGroup"
eso.merge.run = subset(eso.merge.run, idents = c("ESCC I","ESCC II","ESCC III","ESCC IVA") )


cellchat <- createCellChat(object = eso.merge.run, group.by = "celltype_sub", assay = "RNA")

CellChatDB <- CellChatDB.human

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)
# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 2) # do parallel
options(future.globals.maxSize=4000000000)
cellchat <- identifyOverExpressedGenes(cellchat,do.fast = TRUE)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)



Bcell_net <- subset(cellchat, target=='GC.B.dark')
Bcell_net$cell_inter <- paste0(Bcell_net$source,"<-",Bcell_net$target)

p2 = ggplot(Bcell_net,aes(x=cell_inter,y=interaction_name)) +
  geom_point(aes(size=prob,color=prob)) +
  geom_point(shape=21,aes(size=prob))+
  facet_wrap(~group)+
  scale_color_gradientn('Communication\nProbability', 
                        colors=colorRampPalette(rev(brewer.pal(9, "PRGn")))(100)) +
  theme_bw() +
  theme(axis.text=element_text(size=10, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, size=10),
        axis.text.y = element_text(size=8, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
ggplot2::ggave("GCB.cellchat.pdf",p2,heigth=8,width=5)



####to run with manuall added LRs
library(CellChat)
library(plyr)
library(dplyr)
library(stringr)

# load DB
CellChatDB <- CellChatDB.human

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
db <- subsetDB(CellChatDB)
mif_pathways <- subset(db$interaction, ligand == "MIF")

# 
print(unique(mif_pathways$receptor))

mif_cd74 <- mif_pathways[str_detect(mif_pathways$receptor, "CD74") & str_detect(mif_pathways$receptor, "CXCR4"), ]

# 
if(nrow(mif_cd74) == 0){
  stop("No receptor matching 'CD74+CXCR4' found in MIF pathways!")
}

mif_cd74$receptor <- "CD74"
mif_cd74$pathway_name <- "MIF_CD74"

mif_cxcr4 <- mif_pathways[str_detect(mif_pathways$receptor, "CD74") & str_detect(mif_pathways$receptor, "CXCR4"), ]
mif_cxcr4$receptor <- "CXCR4"
mif_cxcr4$pathway_name <- "MIF_CXCR4"

#
db_custom <- db
db_custom$interaction <- rbind(db$interaction, mif_cd74, mif_cxcr4)

subset(db_custom$interaction, ligand == "MIF")[, c("ligand", "receptor", "pathway_name")]
cellchat@DB <- db_custom

# 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
