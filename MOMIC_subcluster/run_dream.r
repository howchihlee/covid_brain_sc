#library(ggplot2)
library(data.table)
library(BiocParallel)
library(limma)
library(edgeR)
library(variancePartition)
library(dplyr)
library(Matrix)
library(compositions)

pdf(file = "figures.pdf") # The height of the plot in inches

geneExpr = fread('./freq_celltype.csv')
metadata = geneExpr[, 10:15]
geneExpr = as.matrix(geneExpr[, 1:9])

metadata$is_covid = factor(metadata$covid)
metadata$subject_id = factor(metadata$subject)
metadata$tissue = factor(metadata$tissue)
metadata$pid = factor(metadata$pid)
metadata$Biosample = factor(with(metadata, paste(subject_id, pid, tissue)))
metadata$aggr_unit = factor(with(metadata, paste(subject_id, pid, tissue_id))) # tissue_id
metadata$ID = factor(with(metadata, paste(subject_id, pid, tissue)))

## dream
info = metadata
mat = clr(t(geneExpr))

form = ~ 0 + is_covid:tissue + (1|tissue) + (1|subject_id) 

L1  = getContrast(mat, form, info, c("is_covid1:tissuePFC",'is_covid0:tissuePFC'))
L2  = getContrast(mat, form, info, c("is_covid1:tissuemedulla",'is_covid0:tissuemedulla'))
L3  = getContrast(mat, form, info, 
  c("is_covid1:tissuechoroidplexus",'is_covid0:tissuechoroidplexus'))
L = data.frame(PFC = L1, medulla = L2, choroidplexus = L3)
fit = dream(mat, form, info, L = L, BPPARAM=SnowParam(24, progressbar=TRUE))

tabList = lapply( colnames(coef(fit))[1:3], function(tissue){
    tab = topTable(fit, coef=tissue, number=Inf)
   
    data.frame( Gene      = rownames(tab),
                Tissue    = tissue,  
                logFC     = tab$logFC,
                P.Value   = tab$P.Value,
                z.std     = tab$z.std)
  })

resTab = do.call(rbind, tabList)
resTab$fdr = p.adjust(resTab$P.Value, 'fdr')
write.csv(resTab, 'dream_test.csv')

dev.off()
