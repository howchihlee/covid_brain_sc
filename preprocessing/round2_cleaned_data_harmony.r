library(Seurat)
llibrary(data.table)
library(dplyr) 
library('org.Hs.eg.db')
library(limma)

seurat_recipet = function(expr){
expr_obj = CreateSeuratObject(expr)
expr_obj = NormalizeData(expr_obj, normalization.method = "LogNormalize", scale.factor = 10000)
expr_obj = ScaleData(expr_obj)
expr_obj = FindVariableFeatures(expr_obj, selection.method = "vst", nfeatures = 2000)
expr_obj = RunPCA(expr_obj, features = VariableFeatures(object = expr_obj), npcs = 100)
expr_obj
}

write_csv = function(vec, fn){write.table(vec, file = fn, append = FALSE, quote = FALSE, sep = ",",
                eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                col.names = TRUE, qmethod = c("escape", "double"),
                fileEncoding = "")}

folder_ = paste0('../processed_data/clean_data/')
exprs=read.csv(paste0(folder_, 'gene_expression.csv'), row.names=1, header=T , check.names = F)
exprs=t(exprs)

meta=read.csv(paste0(folder_, 'metadata.csv'), row.names=1, header=T, check.names = F)
abrain = seurat_recipet(exprs)

for (npc in c(20, 50, 100))
{
    
    df_harmony = read.csv(paste0(folder_, 'harmony', npc ,'.csv' ), row.names = 1)

    abrain[["harmony"]] <- CreateDimReducObject(embeddings = as.matrix(df_harmony), 
                                 key = "harmony_", assay = DefaultAssay(abrain))
    
    abrain1=FindNeighbors(object = abrain, dims = 1:npc, k.param = 30, reduction = "harmony") 
    abrain1=FindClusters(object = abrain1, resolution = 1., graph.name = 'RNA_snn')
    abrain1=RunUMAP(object = abrain1, min.dist = 1.5, graph = 'RNA_snn',  n.epochs = 200) 
    
write_csv(abrain1@meta.data$seurat_clusters, paste0(folder_, "seurat_clusters_harmony", npc, '.csv'))
write_csv(Embeddings(object = abrain1, reduction = "umap"), paste0(folder_, 'seurat_umap_harmony', npc, '.csv'))
    
}