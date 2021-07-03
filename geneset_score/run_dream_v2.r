#library(ggplot2)
library(data.table)
library(BiocParallel)
library(limma)
library(edgeR)
library(variancePartition)
library(dplyr)
library(Matrix)

pdf(file = "figures.pdf") # The height of the plot in inches

geneExpr = fread('./kegg/score_mean_tissue.csv')
metadata = geneExpr[, 1:7]
geneExpr = as.matrix(geneExpr[, 8:ncol(geneExpr)])

metadata$is_covid = factor(metadata$covid)
metadata$subject_id = factor(metadata$subject_id)
metadata$tissue = factor(metadata$tissue)
metadata$pid = factor(metadata$pid)
metadata$celltype = factor(metadata$celltype)
metadata$Biosample = factor(with(metadata, paste(subject_id, pid, tissue)))
metadata$aggr_unit = factor(with(metadata, paste(subject_id, pid, tissue_id,  celltype))) # tissue_id
metadata$ID = factor(with(metadata, paste(subject_id, pid, tissue, celltype)))

## dream
test_per_celltype = function(CT){
  idx = which(info$celltype == CT)
  form = ~ 0 + is_covid:tissue + (1|tissue) + (1|subject_id) 
  L1  = getContrast(mat[,idx], form, info[idx,], 
    c("is_covid1:tissuePFC",'is_covid0:tissuePFC'))
  L2  = getContrast(mat[,idx], form, info[idx,], 
    c("is_covid1:tissuemedulla",'is_covid0:tissuemedulla'))
  L3  = getContrast(mat[,idx], form, info[idx,], 
    c("is_covid1:tissuechoroidplexus",'is_covid0:tissuechoroidplexus'))
  L = data.frame(PFC = L1, medulla = L2, choroidpleus = L3)
  fit = dream(mat[,idx], form, info[idx,], L = L, BPPARAM=SnowParam(24, progressbar=TRUE))
  fit    
    }

info = metadata
mat = t(geneExpr)

fitList = lapply( levels(info$celltype), test_per_celltype)
names(fitList) = levels(info$celltype)

resTab = lapply( names(fitList), function(CT){
  fit = fitList[[CT]]
  tabList = lapply( colnames(coef(fit))[1:3], function(tissue){
    tab = topTable(fit, coef=tissue, number=Inf)
   
    data.frame( Gene      = rownames(tab),
                Tissue    = tissue, 
                CellType  = CT, 
                logFC     = tab$logFC,
                P.Value   = tab$P.Value,
                z.std     = tab$z.std)
  })
  do.call(rbind, tabList)
})
resTab = data.table(do.call(rbind, resTab))
write.csv(resTab, 'kegg/dream_test.csv')

dev.off()
