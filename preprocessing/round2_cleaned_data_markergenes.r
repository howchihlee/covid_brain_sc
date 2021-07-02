library(Seurat)
library(ggplot2)
library(cowplot)
library(data.table)
library(DoubletFinder)
library(pheatmap)
library(future)
library(magrittr) # needs to be run every time you start R and want to use %>%
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

folder_to_read = paste0('../processed_data/clean_data/')
exprs=read.csv(paste0(folder_to_read, 'gene_expression.csv'), row.names=1, header=T , check.names = F)
exprs=t(exprs)

meta=read.csv(paste0(folder_to_read, 'metadata.csv'), row.names=1, header=T, check.names = F)
abrain = seurat_recipet(exprs)

Idents(abrain) = meta$celltype

abrain.markers=FindAllMarkers(object = abrain, min.pct= 0.3)

saveRDS(abrain.markers, file = paste0(folder_to_read, "seurat.marker.list.rds"))
topgene=abrain.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write_csv(topgene, paste0(folder_to_read,'seurat_topgene.csv'))
write_csv(as.data.frame(abrain.markers), paste0(folder_to_read, 'seurat_markers.csv') )

bk_genes = mapIds(org.Hs.eg.db, rownames(abrain@assays$RNA), 'ENTREZID', 'SYMBOL')

for (ct in unique(meta$celltype))
{
    fn_go = paste0(folder_to_read, ct, '_go.csv')
    markers = as.data.frame(abrain.markers[abrain.markers[ ,6] == ct, ])
    
    ind = (markers$p_val_adj < 0.05) & (markers$avg_logFC > 0.2)
    #ind = (markers$p_val_adj < 0.05)
    sgenes = markers[ind, ]
    ind = rev(order(sgenes$avg_logFC))[1:50]
    sgenes = sgenes[ind, ]$gene

    sgene_id = mapIds(org.Hs.eg.db, sgenes, 'ENTREZID', 'SYMBOL')

    go = goana(sgene_id, universe = bk_genes)

    write.table(go, file = fn_go, append = FALSE, sep = ",",
                eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                col.names = TRUE, qmethod = c("escape", "double"),
                fileEncoding = "")
    
    }