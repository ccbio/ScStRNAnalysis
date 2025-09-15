library(Seurat)

######inferCNV
malignant2run = c('ESCC','ESCC I', 'ESCC II', 'ESCC III','ESCC IVA')
eso.merge.run$pheno.inferCNV = "other"
eso.merge.run$pheno.inferCNV[eso.merge.run$SpecimenGroup %in% malignant2run & eso.merge.run$celltype_top == "Epithelial cells"] = "MalignantEpi"
eso.merge.run$pheno.inferCNV[eso.merge.run$SpecimenGroup %in% c("Healthy esophagus") & eso.merge.run$celltype_top == "Epithelial cells"] = "NormalEpi"
eso.merge.run$pheno.inferCNV[eso.merge.run$celltype_top == "B cells"] = "B cells"
eso.merge.run$pheno.inferCNV[eso.merge.run$celltype_top == "Myeloid"] = "Myeloid"
eso.merge.run$pheno.inferCNV[eso.merge.run$celltype_top == "TNK cells"] = "TNK cells"
eso.merge.run$pheno.inferCNV[eso.merge.run$celltype_top == "Mast cells"] = "Mast cells"
anno = data.frame(V1=colnames(eso.merge.run),V2=eso.merge.run$pheno.inferCNV)

allmalignant = which(eso.merge.run$pheno.inferCNV == "MalignantEpi")
selectnormal = sample(which(eso.merge.run$pheno.inferCNV == "NormalEpi"),200)
selectB = sample(which(eso.merge.run$celltype_top == "B cells"),200)
selectMyeloid = sample(which(eso.merge.run$celltype_top == "Myeloid"),200)
selectTNK = sample(which(eso.merge.run$celltype_top == "TNK cells"),200)
selectMast = sample(which(eso.merge.run$celltype_top == "Mast cells"),200)
counts <- GetAssayData(eso.merge.run, slot = 'counts')
for(i in 1:ceiling(length(allmalignant) / 5000)){
  assign.batch = sapply(1:ceiling(length(allmalignant) / 5000),function(x){rep(paste0("batch",x),5000)})[1:length(allmalignant)]
  selectmalignant = which(assign.batch == paste0("batch",i))
  runanno = anno[c(allmalignant[selectmalignant],selectnormal,selectB,selectMyeloid,selectTNK,selectMast),]
  write.table(runanno,"inferCNVfile/tmpanno.txt",row.names=F,col.names=F,sep="\t",quote=F)
  gene.pos.file = "inferCNVfile/esophagus_gene_pos.txt"
  
  library(infercnv)
  
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts,
                                      annotations_file="inferCNVfile/tmpanno.txt",
                                      delim="\t",
                                      gene_order_file=gene.pos.file,
                                      ref_group_names=c("NormalEpi","B cells","Myeloid","TNK cells","Mast cells"))
  rm(eso.merge.run)
  rm(counts)
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir=paste0("inferCNVfile/","batch",i),
                               cluster_by_groups=TRUE,
                               denoise=TRUE,
                               HMM=TRUE)
  
  
  saveobj = paste0(paste0("inferCNVfile/","infercnv_batch",i),".Rdata")
  save(infercnv_obj,file=saveobj)
}

infercnvres = function(datadir=datadir){
  cnvScore <- function(data){
    library(dplyr)
    library(scales)
    data <- data %>% as.matrix() %>%
      t() %>% 
      scale() %>% 
      rescale(to=c(-1, 1)) %>% 
      t()
    
    cnv_score <- as.data.frame(colSums(data * data))
    return(cnv_score)
  }
  allscores = NULL; allcells = NULL
  dirnames = dir(datadir)[grep("^batch",dir(datadir),perl=T)]
  for (dirname in dirnames){
    show(dirname)
    expr <- read.table(paste0(datadir,dirname,"/infercnv.observations.txt"), header=T,check.names = F)
    cnv_score <- cnvScore(expr)
    score0 =  cnv_score[,1]
    cell0  =  rownames(cnv_score)
    allscores = c(allscores,score0)
    allcells = c(allcells,cell0)
  }
  data.frame(cell=allcells,score=allscores)
}
cnv_score = infercnvres(datadir="results/epi/")
write.table(cnv_score,"results/epi/cnv_scores_epithelial_cells.infercnv.txt",sep="\t",row.names=F,quote=F)



############RUN NMF
##functions
###定义相应函数
GetData = function(
    srt,
    genes = NULL,
    slot = 'scale.data',
    assay = NULL
){
  if (is.null(assay)){
    assay = DefaultAssay(srt)
  }
  if (is.null(genes)){
    if ('RNA' %in% names(srt@assays)){
      genes = rownames(GetAssayData(srt, assay = 'RNA', slot = 'counts'))
    } else if  ('Spatial' %in% names(srt@assays)){
      genes = rownames(GetAssayData(srt, assay = 'Spatial', slot = 'counts'))
    } else {
      genes = rownames(GetAssayData(srt, assay = assay, slot = 'counts'))
    }
  }
  data = GetAssayData(srt, assay = assay, slot = slot)
  missing = setdiff(genes, rownames(data))
  add = matrix(0, nrow = length(missing), ncol = ncol(data))
  rownames(add) = missing
  data = rbind(data, add)
  data = data[genes, ]
  return(data)
}
NMFToModules = function(
    res,
    gmin = 5
){
  
  scores = basis(res)
  coefs = coefficients(res)
  
  # Remove if fewer than gmin genes
  ranks_x = t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
  ranks_y = apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)
  for (i in 1:ncol(scores)){
    ranks_y[ranks_x[,i] > 1,i] = Inf
  }
  modules = apply(ranks_y, 2, function(m){
    a = sort(m[is.finite(m)])
    a = a[a == 1:length(a)]
    names(a)
  })
  l = sapply(modules, length)
  keep = (l >= gmin)
  scores = scores[, keep]
  coefs = coefs[keep, ]
  
  # Find modules
  ranks_x = t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
  ranks_y = apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)
  for (i in 1:ncol(scores)){
    ranks_y[ranks_x[,i] > 1,i] = Inf
  }
  modules = apply(ranks_y, 2, function(m){
    a = sort(m[is.finite(m)])
    a = a[a == 1:length(a)]
    names(a)
  })
  
  names(modules) = sapply(modules, '[', 1)
  names(modules) = paste('m', names(modules), sep = '_')
  names(modules) = gsub('-','_',names(modules))
  
  return(modules)
}
BuildEnrichmentMatrix = function(
    genes,
    type = 'GO',
    db = NULL
){
  if (is.null(db)){
    db = FindMSigDB(type)
  }
  terms = names(db)
  enrichment.matrix = sapply(db, function(term){
    genes %in% term
  })
  rownames(enrichment.matrix) = genes
  #enrichment.matrix = enrichment.matrix[, colSums(enrichment.matrix) > 0]
  enrichment.matrix[, colSums(enrichment.matrix) == 0] = NA
  return(enrichment.matrix)
}
FindMSigDB = function(
    type
){
  library(msigdbr)
  if (type == 'GO'){
    gene_sets <- msigdbr(species = "Homo sapiens",category = 'C5')%>% dplyr::select(gs_name,human_entrez_gene, gene_symbol)
    db = split(gene_sets$gene_symbol,f = gene_sets$gs_name)
  } else if (type == 'HALLMARK'){
    gene_sets <- msigdbr(species = "Homo sapiens",category = 'H')%>% dplyr::select(gs_name,human_entrez_gene, gene_symbol)
    db = split(gene_sets$gene_symbol,f = gene_sets$gs_name)
  } else if (type == 'MOTIF'){
    gene_sets <- msigdbr(species = "Homo sapiens",category = 'C3')%>% dplyr::select(gs_name,human_entrez_gene, gene_symbol)
    db = split(gene_sets$gene_symbol,f = gene_sets$gs_name)
    db = db[-grep('UNKNOWN|MIR',names(db))]
  } else if (type == 'PATHWAYS'){
    gene_sets <- msigdbr(species = "Homo sapiens",category = 'C2')%>% dplyr::select(gs_name,human_entrez_gene, gene_symbol)
    db = split(gene_sets$gene_symbol,f = gene_sets$gs_name)
    db = db[grep('BIOCARTA|REACTOME|KEGG',names(db))]
  } else if (type == 'BIOCARTA'){
    gene_sets <- msigdbr(species = "Homo sapiens",category = 'C2')%>% dplyr::select(gs_name,human_entrez_gene, gene_symbol)
    db = split(gene_sets$gene_symbol,f = gene_sets$gs_name)
    db = db[grep('BIOCARTA',names(db))]
  } else if (type == 'KEGG'){
    gene_sets <- msigdbr(species = "Homo sapiens",category = 'C2')%>% dplyr::select(gs_name,human_entrez_gene, gene_symbol)
    db = split(gene_sets$gene_symbol,f = gene_sets$gs_name)
    db = db[grep('KEGG',names(db))]
  } else if (type == 'REACTOME'){
    gene_sets <- msigdbr(species = "Homo sapiens",category = 'C2')%>% dplyr::select(gs_name,human_entrez_gene, gene_symbol)
    db = split(gene_sets$gene_symbol,f = gene_sets$gs_name)
    db = db[grep('REACTOME',names(db))]
  } else {
    gene_sets <- msigdbr(species = "Homo sapiens",category = 'C5')%>% dplyr::select(gs_name,human_entrez_gene, gene_symbol)
    db = split(gene_sets$gene_symbol,f = gene_sets$gs_name)
    db = db[grep(type, names(db), ignore.case = TRUE, value = TRUE)]
  }
  return(db)
}
loadRData = function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}
Heatmap = function(
    matrix,
    colors = NULL,
    breaks = NULL,
    col = NULL,
    clustering_method_rows = 'ward.D2', 
    clustering_method_columns = 'ward.D2',
    clustering_distance_rows = 'pearson',
    clustering_distance_columns = 'pearson',
    show_row_names = FALSE,
    show_column_names = FALSE,
    is.symmetric = FALSE,
    ...
){
  if (is.null(breaks)){
    mi = min(matrix, na.rm = TRUE)
    ma = max(matrix, na.rm = TRUE)
    if (mi < 0 & ma > 0){
      breaks = c(mi, 0, ma)
    } else {
      breaks = c(mi, (mi + ma)/2, ma)
    }
  }
  if (is.null(colors)){
    colors = c('blue','white','red')
  }
  if (is.null(col)){
    col = colorRamp2(breaks = breaks, colors = colors)
  }
  if (is.symmetric == TRUE){
    h = ComplexHeatmap::Heatmap(matrix, 
                                col = col,
                                clustering_method_rows = clustering_method_rows, 
                                clustering_distance_rows = clustering_distance_rows,
                                cluster_columns = FALSE,
                                show_row_names = show_row_names,
                                show_column_names = show_column_names,
                                ...)
    o = unlist(row_order(h))
    return(ComplexHeatmap::Heatmap(matrix, 
                                   col = col,
                                   clustering_method_rows = clustering_method_rows, 
                                   clustering_distance_rows = clustering_distance_rows,
                                   cluster_columns = FALSE,
                                   column_order = o,
                                   show_row_names = show_row_names,
                                   show_column_names = show_column_names, 
                                   ...))
  } else {
    return(ComplexHeatmap::Heatmap(matrix, 
                                   col = col,
                                   clustering_method_rows = clustering_method_rows, 
                                   clustering_method_columns = clustering_method_columns,
                                   clustering_distance_rows = clustering_distance_rows,
                                   clustering_distance_columns = clustering_distance_columns,
                                   show_row_names = show_row_names,
                                   show_column_names = show_column_names,
                                   ...))
  }
  
}
GeneToEnrichment <- function (srt, type = "GO", db = NULL, method = "rand", genes = NULL, 
                              assay = NULL, do.rescale = FALSE, min.cells = 0, min.genes = 0, 
                              min.var = 0, min.var.rescaled = 0, auc_percentile = 0.05, 
                              db_rand = NULL, nrand = 4, nbin = 25, ...) 
{
  if (is.null(assay)) {
    assay = DefaultAssay(srt)
  }
  if (is.null(db)) {
    db = FindMSigDB(type)
  }
  counts = as.matrix(GetData(srt, assay = assay, slot = "counts"))
  genes = rownames(counts)
  genes.expr = rownames(counts)[rowSums(counts) > min.cells]
  if (method == "metagene") {
    data = as.matrix(GetAssayData(srt, assay = assay, slot = "scale.data"))
    db = lapply(db, intersect, genes.expr)
    enrichment.profile = t(sapply(names(db), function(m) {
      colMeans(data[db[[m]], ], na.rm = TRUE)
    }))
    enrichment.profile = enrichment.profile[sapply(names(db), 
                                                   function(x) {
                                                     v = var(enrichment.profile[x, ])
                                                     l = length(db[[x]])
                                                     return(l > min.genes && v > min.var && v * l^2 > 
                                                              min.var.rescaled)
                                                   }), ]
    if (do.rescale) {
      mn = apply(enrichment.profile, 1, mean)
      v = apply(enrichment.profile, 1, var)
      enrichment.profile = (enrichment.profile - mn)/sqrt(v)
    }
    srt = AddMetaData(srt, t(enrichment.profile), col.name = rownames(enrichment.profile))
  }
  if (method == "auc") {
    data = as.matrix(GetData(srt, assay = assay, slot = "data"))
    cells_rankings = AUCell_buildRankings(data)
    cells_AUC = AUCell_calcAUC(db, cells_rankings, aucMaxRank = nrow(cells_rankings) * 
                                 auc_percentile)
    enrichment.profile = getAUC(cells_AUC)
    if (do.rescale) {
      mn = apply(enrichment.profile, 1, mean)
      v = apply(enrichment.profile, 1, var)
      enrichment.profile = (enrichment.profile - mn)/sqrt(v)
    }
    srt = AddMetaData(srt, t(enrichment.profile), col.name = rownames(enrichment.profile))
  }
  if (method == "score") {
    temp = AddModuleScore(srt, features = db, assay = assay, 
                          name = names(db), nbin = nbin, ...)
    enrichment.profile = t(temp@meta.data[, names(db)])
    if (do.rescale) {
      mn = apply(enrichment.profile, 1, mean)
      v = apply(enrichment.profile, 1, var)
      enrichment.profile = (enrichment.profile - mn)/sqrt(v)
    }
    srt = AddMetaData(srt, t(enrichment.profile), col.name = rownames(enrichment.profile))
  }
  if (method == "rand") {
    data = as.matrix(GetData(srt, assay = assay, slot = "scale.data"))
    db = lapply(db, intersect, genes)
    if (is.null(db_rand)) {
      db_rand = MakeRand(srt, db, nrand = nrand, nbin = nbin)
    }
    else {
      nrand = log10(length(db_rand[[1]]))
    }
    enrichment.profile = t(sapply(names(db), function(m) {
      ra = sapply(db_rand[[m]], function(i) {
        colMeans(data[i, ], na.rm = TRUE)
      })
      re = colMeans(data[db[[m]], ], na.rm = TRUE)
      p = rowMeans(ra >= re)
      p = -log10(p)
      return(p)
    }))
    enrichment.profile[is.infinite(enrichment.profile)] = nrand
    enrichment.profile = enrichment.profile/nrand
    srt = AddMetaData(srt, t(enrichment.profile), col.name = rownames(enrichment.profile))
  }
  return(srt)
}
MakeRand <-function (srt, db, assay = NULL, nrand = 3, nbin = 25) 
{
  if (is.null(assay)) {
    assay = DefaultAssay(srt)
  }
  data = GetData(srt, slot = "data")
  db = lapply(db, intersect, rownames(data))
  data.avg = sort(rowMeans(x = data))
  data.cut = cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
                        n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) = names(x = data.avg)
  binned = split(names(data.cut), data.cut)
  db_rand = lapply(names(db), function(m) {
    lapply(1:10^nrand, function(i) {
      used = vector()
      unused = binned
      for (g in db[[m]]) {
        pool = data.cut[g]
        new = sample(unused[[pool]], 1)
        used = c(used, new)
        unused[[pool]] = setdiff(unused[[pool]], new)
      }
      return(used)
    })
  })
  names(db_rand) = names(db)
  return(db_rand)
}




##run NMF by sample
#ESCCepi: cnv score > 0.04
DefaultAssay(ESCCepi) <- 'RNA'
malignant_list = SplitObject(object = ESCCepi, split.by = "orig.ident")
malignant_list_names = names(malignant_list)
saveRDS(malignant_list,file='results/epi/ESCCepi/runNMF/malignant_list.RDS')
saveRDS(malignant_list_names,file='results/epi/ESCCepi/runNMF/malignant_list_names.RDS')

newnames = NULL
for(i in malignant_list_names){
  if(length(ESCCepi$orig.ident[ESCCepi$orig.ident == i]) >= 100 & length(ESCCepi$orig.ident[ESCCepi$orig.ident == i]) <= 1000 ){
    newnames = c(newnames, i)
    saveRDS(malignant_list[[i]],file = paste0('results/epi/ESCCepi/runNMF/malignant_RDS/',i,'_malignant_Seurat.RDS'))
  }
}
malignant_list_names = newnames


finished <- c()
for(i in 1:length(dir("./sample_list/"))){
  finished <- c(finished, strsplit(dir("results/epi/ESCCepi/runNMF/run/")[i], split = "res.list")[[1]][1])
}

library(dplyr)

malignant_list_names <- readRDS("results/epi/ESCCepi/runNMF/malignant_list_names.RDS")
for(i in malignant_list_names) {
  if(i %in% finished){
    print(paste0(i, " is already over!"))
    next;
  }
  
  obj <- readRDS(file = paste0('results/epi/ESCCepi/runNMF/malignant_RDS/',i,'_malignant_Seurat.RDS'))
  obj = NormalizeData(obj)
  obj = FindVariableFeatures(obj)
  obj = ScaleData(obj)
  data <- GetAssayData(obj, slot = 'scale.data') %>% as.matrix()
  
  data <- data[VariableFeatures(obj),]
  data[data < 0] <- 0
  data <- data[apply(data, 1, var) > 0, ]
  res.list <- parallel::mclapply(5:13, function(r){
    nmf(data, rank = r, nrun = 1, seed = 'ica', method = 'nsNMF')
  }, mc.cores = 4)
  names(res.list) = 5:13
  saveRDS(res.list,file = paste0("results/epi/ESCCepi/runNMF/run/",i,"res.list.rds"))
  show(paste0(i, " is already over!"))
}



###read NMF results，extract modules
malignant_list_names <- readRDS("results/epi/ESCCepi/runNMF/malignant_list_names.RDS")
respath <- "results/epi/ESCCepi/runNMF/run/"
res.list<- list()
newnames = NULL
for(i in malignant_list_names){
  file = paste0(respath,i ,"res.list.rds")
  if(file.exists(file)){
    newnames = c(newnames, i)
    a <- readRDS(file)
    res.list[[i]] <- a
  }
}
names(res.list) <- newnames
malignant_list_names = newnames

#### obtain the best NMF clustering rank，extract modules
modules_L <- list()
res_optima_L <- list()
for (i in 1:length(malignant_list_names)){ 
  modules.list = lapply(res.list[[i]], NMFToModules, gmin = 5)
  print(sapply(modules.list,length))
  comp = as.numeric(names(modules.list)) - sapply(modules.list, length)
  mi = min(comp)
  r = names(which(comp == mi))
  r = r[length(r)]
  print(r)
  res = res.list[[i]][[r]]
  modules = NMFToModules(res, gmin = 5)   
  modules_L[[i]] <- modules
  res_optima_L[[i]] <- res  
}
saveRDS(object = modules_L,file = "results/epi/ESCCepi/runNMF/run/modules_L.rds")
saveRDS(object = res_optima_L,file = "results/epi/ESCCepi/runNMF/run/res_optima_L.rds")



##obtain all modules with overlap >5% in >=3 other modules
library(NMF)
modules_L = readRDS("results/epi/ESCCepi/runNMF/run/modules_L.rds")
res_optima_L = readRDS("results/epi/ESCCepi/runNMF/run/res_optima_L.rds")
names(res_optima_L) <- malignant_list_names

modules.list <- lapply(res_optima_L, NMFToModules)

all = unlist(modules.list, recursive = FALSE, use.names = FALSE)
names(all) = unlist(sapply(modules.list, names))
ta = table(unlist(all))
genes.use = names(ta)[ta > 1]

all = unlist(modules.list, recursive = FALSE, use.names = TRUE)
all = lapply(all, intersect, genes.use)
sim = sapply(all, function(x){
  sapply(all, function(y){
    length(intersect(x,y))/length(union(x,y))
  })
})
keep = rownames(sim)[apply(sim, 1, function(x){
  sum(x > 0.05) >= 3 
})]

all = all[keep]
modules.list = lapply(names(modules.list), function(x){
  li = modules.list[[x]]
  li[names(li)[paste(x,names(li),sep='.') %in% keep]]
})
names(modules.list) =  names(res_optima_L)
ta = table(unlist(all))
genes.use = names(ta)[ta > 1] 
print(length(all))

options(repr.plot.height = 7.5, repr.plot.width = 8)
png("modules.png")
pheatmap(sim, show_rownames = F, show_colnames = F)
dev.off()

#final modules

sub = matrix(0, nrow = length(genes.use), ncol = length(genes.use))
rownames(sub) = genes.use
colnames(sub) = genes.use
for (s in names(modules.list)){
  for (mod in modules.list[[s]]){
    mod = intersect(mod, genes.use)
    for (x in mod){
      for (y in mod){
        sub[x,y] = sub[x,y] + 1
      }
    }
  }
}
diag(sub) = 0
# Remove low connections
adj = sub
adj[] = (adj >= 3)
#adj[adj <= 1] = 0
for (i in 1:5){
  keep = names(which(rowSums(adj) >= 5))
  adj = adj[keep,keep]
  print(dim(adj))
}
library(igraph)
g = graph_from_adjacency_matrix(adj, diag = FALSE, mode = 'undirected', weighted = TRUE)
modules = communities(cluster_infomap(g, nb.trials = 100))
names(modules) = paste0('m_', sapply(modules, '[', 1))
save(modules,file = "results/epi/ESCCepi/runNMF/modules.Rdata")







#GO annotation of modules
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
load("results/epi/modules.Rdata")
load("results/epi/ESCCepi/runNMF/modules.Rdata")
# gene symbol to ENTREZID
for(i in 1:length(modules)){
  modules[[i]] <- bitr(modules[[i]], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  modules[[i]] <- modules[[i]][,2]
}
#GO enrichment
ego_GO <- list()
for(i in 1:length(modules)){
  ego_GO[[i]] <- enrichGO(gene = modules[[i]],
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          ont = "ALL", #CC，MF，BP，or ALL for all three
                          pAdjustMethod = "BH",
                          minGSSize = 1,
                          pvalueCutoff = 0.01, #通常为0.01或0.05
                          qvalueCutoff = 0.25,
                          readable = TRUE)
}
for(i in 1:length(modules)){
  name = names(modules)[[i]] 

  write.table(ego_GO[[i]],paste0("results/epi/NMFmodules/module_",i,"_",name,".csv"),sep=",",quote=F,row.names=F)
}
write.table(ego_GO,paste0("results/epi/NMFmodules/modules",".csv"),sep=",",quote=F,row.names=F)





#assign function for modules based on GO annotation summary
module_GO <- c("cEMT","Keratinization","Fatty acid metabolism",
               "Amino acid metabolism","Stress","Cholesterol metabolism","Hypoxia","pEMT","Invasion and metastasis",
               "Tissue homeostasis","Unknown","Interferon response","Invasion and metastasis","Unknown","Invasion and metastasis",
               "Cell cycle","Unknown","Metal response","Cell cycle","Telomere maintenance","Unknown","Unknown")

for(i in 1:length(modules)){
  names(modules)[[i]] <- module_GO[i]
}
sim = sapply(all, function(x){
  sapply(modules, function(y){
    pval = phyper(length(intersect(x, y)), length(x), 2*10^4 - length(x), length(y), lower.tail = FALSE)
    return(-log10(pval))
  })
})
sim[is.infinite(sim)] = max(sim[is.finite(sim)])
sim = sim[,apply(sim, 2, function(x){any(x > 3)})]
df = data.frame('sample' = sapply(colnames(sim), function(x){
  y = sapply(strsplit(x, '.', fixed = TRUE), '[', 1)
}))
df$top = apply(sim, 2, which.max)
df$top = factor(names(modules)[df$top], levels = unique(names(modules)))


#heatmap

ovlp = sapply(all, function(x){
  sapply(all, function(y){
    pval = phyper(length(intersect(x, y)), length(x), 2*10^4 - length(x), length(y), lower.tail = FALSE)
    return(-log10(pval))
  })
})

ovlp[is.infinite(ovlp)] = max(ovlp[is.finite(ovlp) == TRUE])
ovlp = ovlp[colnames(sim),colnames(sim)]


save.image("results/epi/ESCCepi/runNMF/runNMF.Rdata")

set.seed(123)
library(ComplexHeatmap); library(RColorBrewer); library(circlize)
#https://r-charts.com/color-palettes/ paletteer::paletteer_d("ggthemes::Classic_20")
# a = as.character(paletteer::paletteer_d("ggthemes::Classic_20"))   ,  paste(a,collapse="\",\"")
col = c("#1F77B4FF","#AEC7E8FF","#FF7F0EFF","#FFBB78FF","#2CA02CFF","#98DF8AFF","#D62728FF","#FF9896FF","#9467BDFF",
        "#C5B0D5FF","#8C564BFF","#C49C94FF","#E377C2FF","#F7B6D2FF","#7F7F7FFF")
names(col) = unique(names(modules))
side_ann = ComplexHeatmap::HeatmapAnnotation(module = df[, c("top")],  col = list(module=col) ,
                                             which = 'row')

h = Heatmap(name = 'Module overlp p-value', ovlp,column_order = order(df$top),  row_order = order(df$top),
            breaks = seq(0, 6,length = 7),colors =rev(RColorBrewer::brewer.pal(7, "RdBu"))
) + side_ann














###module genes selection
module_signature <- list(
  "cEMT" = c("COL1A1","FN1","COL4A2","COL6A1","MMP2","ACTA2","ACTG2","TP63","SFRP2","LUM"),
  "Keratinization" = c("KLK5","CNFN","LY6D","KRT1","KRT5","KRT6B","KRT17","SPRR1A","TGM1","TMEM79"),
  "Fatty acid metabolism" = c("ADH7","AKR1B1","CBR1","CES1","CYP4F3","EPHX1","GSTA1","GSTM1","GSTM4","PTGR1"),
  "Amino acid metabolism" = c("ASNS","ASS1","BCAT1","CTH","GFPT1","PHGDH","PSAT1","PSPH","PYCR1","SHMT2"),
  "Stress" = c("DNAJA1","HSPA1A","HSPA1B","JUN","FOSB","FOS","PTGS2","RHOB","TNF","TNFAIP3"),
  "Cholesterol metabolism" = c("MSMO1","DHCR7","MVD","INSIG1","IDI1","CYP51A1","HMGCR","HMGCS1","HSD17B7","TM7SF2","FDFT1","FDPS","CES1","MVK",
                                 "SC5D","LSS","DHCR24","G6PD","SREBF2","ACLY","SQLE","ACAT2"),
  "Hypoxia" = c("HILPDA","VEGFA","CA9","PGF","HK2","ENO1","NDRG1","PGK1","PDK1"),
  "pEMT" = c("LAMC2","VIM","ITGB1","RGCC","MMP9","TGFB1I1","THBS1","KRT14","DKK1","WNT10A"),
  "Invasion and metastasis" = c("CEACAM1","MACC1","MUC16","LGALS3","IGFBP5","SERPINE1","TIMP3","ANGPTL4","NUAK2","SOX2"),
  "Tissue homeostasis" = c("ALDH1A1","AZGP1","CLDN3","LTF","PIGR","PROM1","PRR4","SERPINA3","SLC12A2","TFF3"),
  "Interferon response" = c("IFIT1","HLA-A","HLA-DRA","STAT1","IRF1","CD74","GBP2","PARP9","LGALS9","IFITM3"),
  "Cell cycle" = c("MKI67","TOP2A","PCNA","AURKB","BUB1B","KIF2C","CCNB1","CDC25C","TTK","UBE2C"),
  "Metal response" = c("MT1E","MT1F","MT1G","MT1H","MT1M"),
  "Telomere maintenance" = c("DKC1","CCT5","NME2")
)

load("results/epi/ESCCepi.obj.Rdata")
ESCCepi = AddModuleScore(ESCCepi, features = module_signature, name = "Module_")
colnames(ESCCepi@meta.data)[(length(colnames(ESCCepi@meta.data))-length(module_signature)+1):length(colnames(ESCCepi@meta.data))] = paste0("Module_",gsub(" ","_",names(module_signature)))

library(ggplot2)
p = NULL
mo = colnames(ESCCepi@meta.data)[c(32:39, 41:43,45:46,48)]
col = c("#1F77B4FF","#AEC7E8FF","#FF7F0EFF","#FFBB78FF","#2CA02CFF","#98DF8AFF","#D62728FF","#FF9896FF","#9467BDFF",
        "#C5B0D5FF","#8C564BFF","#C49C94FF","#E377C2FF","#F7B6D2FF","#7F7F7FFF")

px = FeaturePlot(ESCCepi, features = mo, 
                 reduction = "umap", cols = c("lightgrey", "darkgreen"),ncol=4)
px$layers[[1]]$mapping$alpha <- 0.4
px <- px + scale_alpha_continuous(range = 0.4, guide = "none")
ggplot2::ggsave("results/epi/ESCCepi/ESCCepi_ModuleScore.pdf", px, width=16 ,height=16)




ESCCepi$tmp = ESCCepi$SpecimenGroup
ESCCepi$tmp = factor(ESCCepi$tmp,levels=c("ESCC I","ESCC II","ESCC III","ESCC IVA"))
Idents(ESCCepi) = "tmp"
for (i in 1:length(mo)){
  px = RidgePlot(ESCCepi, features = mo[i])
  ggplot2::ggsave(paste0("results/epi/ESCCepi/ESCCepi_ModuleScore_SpecimenGroup",mo[i],".pdf"), px, width=10 ,height=6)
}

sampleinfo = read.csv("../../../results/cells_unique.csv",sep=",",header=T,stringsAsFactors = F)
select = match(ESCCepi$orig.ident, sampleinfo$orig.ident)
ESCCepi@meta.data$Treatment = sampleinfo$Treatment[select]
ESCCepi$tmp = ESCCepi$Treatment
ESCCepi$tmp = factor(ESCCepi$tmp,levels=c("ESCC I","ESCC II","ESCC III","ESCC IVA"))
Idents(ESCCepi) = "tmp"
for (i in 1:length(mo)){
  mscore = ESCCepi@meta.data[,31+i]
  mscore = ifelse(mscore>0,"Yes","No")
  plotC <- reshape2::melt(table(ESCCepi$tmp,mscore ))
  colnames(plotC) <- c("Sample", "CellType","Number")
  plotC$CellType = factor(plotC$CellType,levels=c("No","Yes"))
  plotC$Prop = apply(plotC,1,function(x){ paste0(round(as.numeric(x[3]) / sum(plotC$Number[plotC$Sample == x[1]]),4)*100,"%")})
  pC2 = ggplot(data = plotC, aes(x = Sample, y = Number, fill = CellType)) +
    geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill")+theme_bw()+
    #scale_fill_manual(values=colors) + 
    #geom_text(aes(label = Prop), stat="identity",colour = "black",position = position_fill(vjust = 0.5)) +
    theme(panel.grid =element_blank()) +
    labs(x="",y="Cell proportion")+
    scale_y_continuous(labels = c(0,0.25,0.5,0.75,1))+ ####用来将y轴移动位置
    theme(axis.text = element_text(size=12, colour = "black"))+
    theme(axis.title.y = element_text(size=12, colour = "black"))+
    theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
    theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.8))
  pC2
  ggplot2::ggsave(paste0("results/epi/ESCCepi/ESCCepi_ModuleScore",mo[i],".pdf"), pC2, width=6 ,height=6)
}

