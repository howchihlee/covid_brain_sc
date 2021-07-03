library(Seurat)
library(qusage)

seurat_recipet = function(expr){
expr_obj = CreateSeuratObject(expr)
expr_obj = NormalizeData(expr_obj, normalization.method = "LogNormalize", scale.factor = 10000)
expr_obj = ScaleData(expr_obj)
expr_obj = FindVariableFeatures(expr_obj, selection.method = "vst", nfeatures = 2000)
expr_obj = RunPCA(expr_obj, features = VariableFeatures(object = expr_obj), npcs = 100)
expr_obj
}

write_csv = function(df, fn){
write.table(df, file = fn, append = FALSE, sep = ",",
                eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                col.names = TRUE, qmethod = c("escape", "double"),
                fileEncoding = "")
}


folder_to_read = paste0('../processed_data/clean_data/')
exprs=read.csv(paste0(folder_to_read, 'gene_expression.csv'), row.names=1, header=T , check.names = F)
exprs=t(exprs)

meta=read.csv(paste0(folder_to_read, 'metadata.csv'), row.names=1, header=T, check.names = F)

expr_obj = seurat_recipet(exprs)

gss = read.gmt('./h.all.v7.1.symbols.gmt')

res = c()
for (i in  seq_along(gss))
{
gs = list(gss[[i]])
expr_obj = AddModuleScore(expr_obj, gs, name = 'mod')
res = rbind(res, expr_obj$mod1)
}

rownames(res) = names(gss)


write_csv(res, 'hallmark/score.csv')


gss = read.gmt('./c2.cp.kegg.v7.1.symbols.gmt')

res = c()
for (i in  seq_along(gss))
{
gs = list(gss[[i]])
expr_obj = AddModuleScore(expr_obj, gs, name = 'mod')
res = rbind(res, expr_obj$mod1)
}

rownames(res) = names(gss)
write_csv(res, 'kegg/score.csv')