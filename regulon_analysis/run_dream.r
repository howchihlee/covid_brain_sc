#library(ggplot2)
library(data.table)
library(BiocParallel)
library(limma)
library(edgeR)
library(variancePartition)
library(dplyr)
library(Matrix)

pdf(file = "figures.pdf") # The height of the plot in inches

geneExpr = fread('./output/pyscenic_aucell_output.csv') %>% as.data.frame
rownames(geneExpr) = geneExpr$V1
geneExpr = geneExpr[,-1]
geneExpr = as.matrix(geneExpr)

metadata = fread('[PATH]/metadata.csv')
metadata$is_covid = factor(metadata$is_covid)
metadata$subject_id = factor(metadata$subject_id)
metadata$tissue = factor(metadata$tissue)
metadata$pid = factor(metadata$pid)
metadata$celltype = factor(metadata$celltype)
metadata$Biosample = factor(with(metadata, paste(subject_id, pid, tissue)))
metadata$aggr_unit = factor(with(metadata, paste(subject_id, pid, tissue_id,  celltype))) # tissue_id
metadata$ID = factor(with(metadata, paste(subject_id, pid, tissue, celltype)))

aggregate_counts = function( geneExpr, metadata){
  # aggrecate gene counts
  geneCounts = lapply( levels(metadata$aggr_unit), function(id){
      idx = which(metadata$aggr_unit == id)
      colMeans(geneExpr[idx,,drop=FALSE])
    })
  geneCounts = do.call(cbind, geneCounts)
  colnames(geneCounts) = levels(metadata$aggr_unit)
  # simplify metadata 
  md2 = lapply( levels(metadata$aggr_unit), function(id){
      unique(metadata[metadata$aggr_unit == id,c("subject_id", "Biosample", "clinical", "is_covid", "pid", "tissue", "celltype")])
  })
  md2 = do.call(rbind, md2)
  list( geneCounts  = geneCounts,
        metadata    = md2) 
}
res = aggregate_counts( geneExpr, metadata)
form = ~ (1|tissue) + (1|celltype) + (1|subject_id)
#vobj = voomWithDreamWeights(res$geneCounts, form, res$meta, plot=TRUE, BPPARAM=SnowParam(24, progressbar=TRUE))


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

info = res$metadata
mat = res$geneCounts

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
dev.off()
