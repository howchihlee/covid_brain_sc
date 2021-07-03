library(BiocParallel)
library(tradeSeq)
library(mgcv)
library(clusterExperiment)
library(slingshot)

exprs = read.csv('./data/mic_gene_expression_raw.csv', row.names = 1)
counts = as.matrix(exprs)
meta = read.csv('./data/mic_diff_meta.csv')
rd = read.csv('./mic_diff_pca.csv')

ind_gene = colSums(counts) > (0.5 * nrow(counts) )
counts = counts[,ind_gene]

sds <- slingshot(rd[,1:2], meta$leiden)
pseudotime = slingPseudotime(sds)[,1]

#aicK = evaluateK(t(counts), sds = sds, BPPARAM = BiocParallel::bpparam())



sce <- fitGAM(counts = t(counts), pseudotime = pseudotime, 
              cellWeights = rep(1, nrow(counts)), verbose = TRUE,
              parallel=TRUE, BPPARAM = BiocParallel::bpparam())

assoRes <- associationTest(sce)
assoRes$fdr = p.adjust(assoRes$pvalue, 'fdr')
head(assoRes)

saveRDS(assoRes, 'mic_sling_output/pseudotime_association.rds')

sigRes = assoRes[assoRes$fdr < 0.05,]
degs = rownames(sigRes[order(sigRes$pvalue),])
nPointsClus = 20
clusPat = clusterExpressionPatterns(sce, nPoints = nPointsClus, 
                                    genes = degs, clusterFunction = "hierarchical01", 
                                    ncores = 20,
                                    )

write.csv(clusPat$yhatScaled,'mic_sling_output/yhat_scaled.csv')

dn_curves = (which(colSums(diff(t(clusPat$yhatScaled), ) <= 0, ) == (nPointsClus - 1) ))
up_curves = (which(colSums(diff(t(clusPat$yhatScaled), ) >= 0, ) == (nPointsClus - 1) ))

not_monotone = setdiff(1:length(degs), dn_curves)
not_monotone = setdiff(not_monotone, up_curves)

not_monotone_genes = rownames(clusPat$yhatScaled)[not_monotone]
clusPat_nm = clusterExpressionPatterns(sce, nPoints = nPointsClus, 
                                    genes = not_monotone_genes, 
                                    clusterFunction = "hierarchical01", 
                                    ncores = 20,
                                    )


clusterLabels_nm <- primaryCluster(clusPat_nm$rsec)
clusterLabels_nm = clusterLabels_nm

out = data.frame(cluster = clusterLabels_nm - 1 )
rownames(out) = rownames(clusPat_nm$yhatScaled)
write.csv(out, './mic_sling_output/clusterLabel_nm.csv')
write.csv(clusPat_nm$yhatScaled, './mic_sling_output/yhat_scaled_nm.csv')

#merged_label = seq(length(degs))
merged_label = rep(-1, length(degs))
merged_label[not_monotone] = clusterLabels_nm - 1
merged_label[dn_curves] = max(clusterLabels_nm)
merged_label[up_curves] = max(clusterLabels_nm) + 1

out = data.frame(cluster = merged_label )
rownames(out) = rownames(clusPat$yhatScaled)
write.csv(out, './mic_sling_output/clusterLabel.csv')
write.csv(clusPat$yhatScaled, './mic_sling_output/yhat_scaled.csv')





pdf('fig/decreasing_genes_smooth_traj.pdf')

genes = c('MT.CO1', 'MT.ND4', 'CSMD1')
for (g in genes){
#print(plotSmoothers(sce, t(counts), gene = g))
 print(plotSmoothers(sce, t(counts), gene = g))
}
dev.off()
