#library(ggplot2)
library(data.table)
library(BiocParallel)
library(limma)
library(edgeR)
library(variancePartition)
library(dplyr)
library(Matrix)
library(cowplot)

pdf(file = "figures.pdf") # The height of the plot in inches

metadata = read.csv('[PATH]/metadata.csv', row.names = 1)

geneExpr = lapply( unique(metadata$celltype), function(ct){
        fn = paste0('[PATH]/knn_diff/gene_10k/gene_expression_smoothed_k30_d50_', ct, '.rds')
        mat = readRDS(fn)
        mat
    })
geneExpr = do.call(rbind, geneExpr)

id2bc = rownames(geneExpr)
geneExpr = as.matrix(geneExpr)


metadata = metadata[id2bc, ]
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
      colSums(geneExpr[idx,,drop=FALSE])
    })
  geneCounts = do.call(cbind, geneCounts)
  colnames(geneCounts) = levels(metadata$aggr_unit)
  # simplify metadata 
  md2 = lapply( levels(metadata$aggr_unit), function(id){
      unique(metadata[metadata$aggr_unit == id,c("subject_id", "Biosample", "clinical", 
                                                 "is_covid", "pid", "tissue", "celltype")])
  })
  md2 = do.call(rbind, md2)
  list( geneCounts  = geneCounts,
        metadata    = md2) 
}

res = aggregate_counts(geneExpr, metadata)
dge = DGEList( res$geneCounts, samples=res$metadata )
dge = calcNormFactors(dge)

form = ~ (1|tissue) + (1|celltype) + (1|subject_id)
vobj = voomWithDreamWeights(dge, form, dge$samples, plot=TRUE, BPPARAM=SnowParam(24, progressbar=TRUE))
saveRDS(vobj, 'voom_obj.rds')

```
########################## ## variancePartition ########################## 
## variancePartition
form = ~ (1|tissue) + (1|celltype) + (1|is_covid) + (1|subject_id) + (1|pid)
vp = fitExtractVarPartModel(vobj, form, dge$samples, BPPARAM=SnowParam(24, progressbar=TRUE))
plotVarPart( sortCols(vp) )

## within cell type
#{r vp.aggregate.within}
vpList = lapply( levels(dge$samples$celltype), function(CT){
  idx = which(dge$samples$celltype == CT)
  form = ~ (1|is_covid) + (1|tissue) + (1|subject_id) 
  
  fitExtractVarPartModel(vobj[,idx], form, dge$samples[idx,],BPPARAM=SnowParam(6, progressbar=TRUE))
})
names(vpList) = levels(dge$samples$celltype)


#{r plot.vp.CT.summary.agg, fig.width=8, fig.height=20}
figList = lapply( names(vpList), function(ID){
  plotVarPart( vpList[[ID]] ) + ggtitle(ID)
  })
pdf('var_par_celltype.pdf')
plot_grid(plotlist=figList[1:6], ncol=3)
plot_grid(plotlist=figList[7:12], ncol=3)
plot_grid(plotlist=figList[13:15], ncol=3)

########################## ## variancePartition ########################## 

## dream
fitList = lapply( levels(dge$samples$celltype), function(CT){
  idx = which(dge$samples$celltype == CT)
  form = ~ 0 + is_covid:tissue + (1|tissue) + (1|subject_id) 
  
  L1  = getContrast(vobj[,idx], form, dge$samples[idx,], 
    c("is_covid1:tissuePFC",'is_covid0:tissuePFC'))
  L2  = getContrast(vobj[,idx], form, dge$samples[idx,], 
    c("is_covid1:tissuemedulla",'is_covid0:tissuemedulla'))
  L3  = getContrast(vobj[,idx], form, dge$samples[idx,], 
    c("is_covid1:tissuechoroidplexus",'is_covid0:tissuechoroidplexus'))
  L = data.frame(PFC = L1, medulla = L2, choroidplexus = L3)
  fit = dream(vobj[,idx], form, dge$samples[idx,], L = L, BPPARAM=SnowParam(24, progressbar=TRUE))
  fit
})
names(fitList) = levels(dge$samples$celltype)

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
write.csv(resTab, 'dream_degs.csv', row.names = F)


countUP  = resTab[,data.frame(nDE = sum((p.adjust(P.Value, 'fdr')  < 0.05) & (logFC > 0))), by=c('Tissue', 'CellType')]
countDN  = resTab[,data.frame(nDE = sum((p.adjust(P.Value, 'fdr')  < 0.05) & (logFC < 0))), by=c('Tissue', 'CellType')]
countDE = data.frame(CellType = countUP$CellType, Tissue = countUP$Tissue, UP = countUP$nDE, DOWN=countDN$nDE)

## plot DEs ##
ymax = max(countUP$nDE)*1.05
ggplot(countUP, aes(CellType, nDE, fill=Tissue)) + geom_bar(stat="identity", position = "dodge") + theme_bw() + theme(aspect.ratio=2, plot.title = element_text(hjust = 0.5)) + ggtitle("Number of UP DE genes") + coord_flip() + scale_y_continuous(expand=c(0,0), limits=c(0, ymax))

ymax = max(countDN$nDE)*1.05
ggplot(countDN, aes(CellType, nDE, fill=Tissue)) + geom_bar(stat="identity", position = "dodge") + theme_bw() + theme(aspect.ratio=2, plot.title = element_text(hjust = 0.5)) + ggtitle("Number of UP DE genes") + coord_flip() + scale_y_continuous(expand=c(0,0), limits=c(0, ymax))

## plot DEs ##
dev.off()
