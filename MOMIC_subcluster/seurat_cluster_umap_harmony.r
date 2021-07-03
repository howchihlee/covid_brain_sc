library(Seurat)
llibrary(data.table)
library(dplyr) 
library('org.Hs.eg.db')
library(limma)

seurat_recipet = function(expr){
expr_obj = CreateSeuratObject(expr)
#expr_obj = NormalizeData(expr_obj, normalization.method = "LogNormalize", scale.factor = 10000)
#expr_obj = ScaleData(expr_obj)
#expr_obj = FindVariableFeatures(expr_obj, selection.method = "vst", nfeatures = 2000)
#expr_obj = RunPCA(expr_obj, features = VariableFeatures(object = expr_obj), npcs = 100)
expr_obj
}

fn_expr = 'MO_round2.csv'
fn_meta = 'MO_meta_round2.csv'
fn_hpca = 'hpca_round2.csv'
fn_umap = 'umap_round2.csv'
fn_cluster = 'seurat_cluster_round2.csv'
fn_topgene = 'seurat_topgene_round2.csv'
fn_marker = 'seurat_markers_round2.csv'

write_csv = function(vec, fn){write.table(vec, file = fn, append = FALSE, quote = FALSE, sep = ",",
                eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                col.names = TRUE, qmethod = c("escape", "double"),
                fileEncoding = "")}

folder_ = paste0('./')
exprs=read.csv(paste0(folder_, fn_expr), row.names=1, header=T , check.names = F)
exprs=t(exprs)

meta=read.csv(paste0(folder_, fn_meta), row.names=1, header=T, check.names = F)
abrain = seurat_recipet(exprs)

npc = 50

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
