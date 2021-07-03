library(scran)

write_csv = function(vec, fn){write.table(vec, file = fn, append = FALSE, quote = FALSE, sep = ",",
                eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                col.names = TRUE, qmethod = c("escape", "double"),
                fileEncoding = "")}

folder_to_read = paste0('../processed_data/clean_data/')
meta=read.csv(paste0(folder_to_read, 'metadata.csv'), row.names=1, header=T, check.names = F)

for (ct in levels(meta$celltype)){  
    folder_to_read = paste0('../processed_data/clean_data/')
    exprs=read.csv(paste0('gene_expression_smoothed_k30_d50_', ct,'.csv'), 
                   row.names=1, header=T , check.names = F)
    saveRDS(exprs, file=paste0('gene_expression_smoothed_k30_d50_', ct,'.rds'))
    exprs=t(exprs)


    sce <- SingleCellExperiment(list(counts=exprs))
    sce <- computeSumFactors(sce, BPPARAM=MulticoreParam(workers = 10))
    sf <- sizeFactors(sce)
    normmat <- sweep(exprs, 2, sf,'/')
    normmat <- log2(normmat + 1)


    saveRDS(t(normmat), file=paste0('knnsmoothed_normed_log2_', ct,'.rds'))
    write_csv (as.data.frame(t(normmat)), paste0('knnsmoothed_normed_log2_', ct,'.csv'))

    }