library(Seurat)
library(cowplot)
library(data.table)
library(future)
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr) 

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

exprs=read.table('GSE150903_sub.csv', row.names=1, header=T , check.names = F)
exprs=t(exprs)

meta=read.csv('GSE150903/meta.tsv', sep = '\t')

abrain = seurat_recipet(exprs)
Idents(abrain) = meta$Cluster
abrain.markers=FindAllMarkers(object = abrain)

saveRDS(abrain.markers, file = "seurat.marker.list.rds")
topgene=abrain.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write_csv(topgene, 'seurat_topgene.csv')
write_csv(as.data.frame(abrain.markers), 'seurat_markers.csv' )

