library(Seurat)
library(Seurat)
llibrary(data.table)
library(dplyr) 
library('org.Hs.eg.db')
library(limma)

write_csv = function(vec, fn){write.table(vec, file = fn, append = FALSE, quote = FALSE, sep = ",",
                eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                col.names = TRUE, qmethod = c("escape", "double"),
                fileEncoding = "")}

seurat_recipet = function(expr){
expr_obj = CreateSeuratObject(expr)
expr_obj = NormalizeData(expr_obj, normalization.method = "LogNormalize", scale.factor = 10000)
expr_obj = ScaleData(expr_obj)
expr_obj = FindVariableFeatures(expr_obj, selection.method = "vst", nfeatures = 2000)
expr_obj = RunPCA(expr_obj, features = VariableFeatures(object = expr_obj), npcs = 100)
expr_obj
}

folder_to_save = './seurat/'
exprs=read.csv('expr_v10k.csv', row.names=1, header=T , check.names = F)
exprs=t(exprs)

meta=read.csv('meta_v10k', row.names=1, header=T, check.names = F)
abrain = seurat_recipet(exprs)
write_csv(Embeddings(object = abrain, reduction = "pca"), paste0(folder_to_save, 'seurat_pca.csv'))
    
for (npc in c(20, 50, 100))
{
    abrain1=FindNeighbors(object = abrain, dims = 1:npc, k.param = 100) 
    abrain1=FindClusters(object = abrain1, resolution = 1., graph.name = 'RNA_snn')
    abrain1=RunUMAP(object = abrain1, min.dist = 1.5, graph = 'RNA_snn') 

    write_csv(abrain1@meta.data$seurat_clusters, paste0(folder_to_save, "seurat_clusters_pc", npc, '.csv') )
    write_csv(Embeddings(object = abrain1, reduction = "umap"), paste0(folder_to_save, 'seurat_umap_pc', npc, '.csv') )
}

npc = 50

fn_cluster = fn_hpca = paste0('./seurat/harmony_cluster', npc, '.csv')
fn_topgene = 'seurat/topgene_r1.csv'
fn_marker = 'seurat/markers_r1.csv'

fn_hpca = paste0('./seurat/harmony', npc, '.csv')

fn_umap = paste0('./seurat/umap_harmony', npc, '.csv')

df_harmony = read.csv(fn_hpca, row.names = 1)

abrain[["harmony"]] <- CreateDimReducObject(embeddings = as.matrix(df_harmony), 
                                 key = "harmony_", assay = DefaultAssay(abrain))
    
abrain1=FindNeighbors(object = abrain, dims = 1:npc, k.param = 15, reduction = "harmony") 
abrain1=FindClusters(object = abrain1, resolution = .1, graph.name = 'RNA_snn')
write_csv(abrain1@meta.data$seurat_clusters, fn_cluster)

abrain1=RunUMAP(object = abrain1, min.dist = .3, graph = 'RNA_snn',  n.epochs = 200) 
write_csv(Embeddings(object = abrain1, reduction = "umap"), fn_umap)


abrain.markers=FindAllMarkers(object = abrain1, min.pct= 0.3, return.thresh = 1.)
topgene=abrain.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write_csv(topgene,fn_topgene)
write_csv(as.data.frame(abrain.markers), fn_marker)

abrain.markers=FindAllMarkers(object = abrain1, min.pct= 0.3, return.thresh = 1.)
topgene=abrain.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write_csv(topgene,fn_topgene)
write_csv(as.data.frame(abrain.markers), fn_marker)

