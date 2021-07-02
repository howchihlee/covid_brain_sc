library(Seurat)
library(ggplot2)
library(cowplot)
library(data.table)
library(DoubletFinder)
library(pheatmap)

rm(list=ls())
gc()

scale_colour_discrete=function(...) {scale_colour_brewer(..., palette="Set2")}
scale_fill_discrete=function(...) {scale_fill_brewer(...,palette="Set2")}
color=c(RColorBrewer::brewer.pal(8,"Set2"),RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(6,"Set3"),RColorBrewer::brewer.pal(10,"RdYlGn"),RColorBrewer::brewer.pal(10,"BrBG"),RColorBrewer::brewer.pal(10,"PiYG"),RColorBrewer::brewer.pal(10,"PuOr"),RColorBrewer::brewer.pal(10,"RdBu"),RColorBrewer::brewer.pal(10,"PRGn"))

write_csv = function(vec, fn){write.table(vec, file = fn, append = FALSE, quote = FALSE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")}

input='./'
setwd(input)
exprs=read.csv(paste0(input,'data/gene_expression.csv'), row.names=1, header=T , check.names = F)
exprs=t(exprs)
meta=read.csv(paste0(input,'data/metadata.csv'), row.names=1, header=T, check.names = F)
meta=meta[colnames(exprs),]
subject=read.csv(paste0(input,'data/subject_info.csv'), header=T, check.names = F)
group=as.character(subject$)

# remove doublet
# here doubled are similated within each subject, consider whether need to be changed based on the sequencing library
output = '[PATH]/Seurat_analysis/'
output1=paste0(output,'1.doublets/')
dir.create(output1,showWarning=F)

per = 0.01  # percentage of doublets which need to be removed, may need to be changed
usample=unique(as.character(meta$subject_id))

for (i in 1:length(usample)) {

  toutput=paste0(output1,usample[i],'/')
    dir.create(toutput,showWarning=F)
    setwd(toutput)
    sdata=exprs[,which(as.character(meta$subject_id) %in% usample[i])]

  brain = CreateSeuratObject(sdata)
  brain = NormalizeData(brain, normalization.method = "LogNormalize", scale.factor = 10000)
  brain = ScaleData(brain)
  brain = FindVariableFeatures(brain, selection.method = "vst", nfeatures = 2000)
  brain = RunPCA(brain, features = VariableFeatures(object = brain), npcs = 50)
  brain = RunUMAP(brain, dims = 1:30)
  sweep.res.list = paramSweep_v3(brain, PCs=1:30, sct = FALSE)
  sweep.stats = summarizeSweep(sweep.res.list, GT = FALSE)
  pdf(paste0("1.find.pK.pdf"))
  bcmvn = find.pK(sweep.stats)
  dev.off()
  pk = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  nExp_poi = round(per*ncol(brain))

  brain = doubletFinder_v3(brain, PCs=1:30, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

  brain@meta.data[,paste0("DF_hi.lo_",per)] = brain@meta.data[,paste0('DF.classifications_0.25_',pk,'_',nExp_poi)]
  brain@meta.data[,paste0("pANN_",per)] = brain@meta.data[,paste0('pANN_0.25_',pk,'_',nExp_poi)]

  pdf(paste0("0.doublet.",per,".pdf"), width = 5, height = 5)
  print(DimPlot(object = brain, reduction = "umap", group.by=paste0("DF_hi.lo_",per), order=c("Doublet","Singlet"), cols=c("black","red")))
  print(FeaturePlot(object = brain, reduction = "umap", features = c("nCount_RNA"), cols = c("grey", "blue")))
  print(FeaturePlot(object = brain, reduction = "umap", features = c("nFeature_RNA"), cols = c("grey", "blue")))
  dev.off()
  singletname=rownames(brain@meta.data[which(brain@meta.data[,paste0("DF_hi.lo_",per)] %in% 'Singlet'),])
  cbrain.data=sdata[,singletname]
  #smeta=FetchData(object = brain, vars = c("nCount_RNA", "nFeature_RNA", "percent.mito", "pANN_0.075", "DF_hi.lo_0.075"))
  smeta=FetchData(object = brain, vars = c("nCount_RNA", "nFeature_RNA", paste0("pANN", per), paste0("DF_hi.lo_",per)))
    
  if (i==1) {
        abrain.data=cbrain.data
    } else {
        abrain.data=cbind(abrain.data,cbrain.data)
    }

    if (i==1) {
        abrain.meta=smeta
    } else {
        abrain.meta=rbind(abrain.meta,smeta)
    }
}

setwd(output1)
saveRDS(abrain.data, file=paste0('0.exprs.matrix.rm.mit.doublets',per,'.all.rds'))
write_csv(abrain.meta, 'doublets.csv')


# build the Seurat object
soutput=paste0(output,'2.seurat.analysis/')
dir.create(soutput,showWarning=F)
setwd(soutput)

usample=c("NPBB113","NPBB052","NPBB098","NPBB105","NPBB135","NPBB138","NPBB139","NPBB137","NPBB131")
group1=c('Ctrl','Ctrl','Ctrl','Ctrl','COVID','COVID','COVID','COVID','COVID_HIV')
group2=c('Ctrl','Ctrl','Ctrl','Ctrl','COVID','COVID','COVID','COVID','COVID')
names(group1)=usample
names(group2)=usample

abrain.meta.data=meta[colnames(abrain.data),]
abrain.meta.data$group1=group1[as.character(abrain.meta.data$subject_id)]
abrain.meta.data$group2=group2[as.character(abrain.meta.data$subject_id)]

abrain=CreateSeuratObject(counts = abrain.data, meta.data=abrain.meta.data)
abrain[["percent.mito"]] = PercentageFeatureSet(object = abrain, pattern = "^MT-")
abrain=NormalizeData(object = abrain, normalization.method = "LogNormalize", scale.factor = 10000)
abrain=FindVariableFeatures(object = abrain, selection.method = "vst", nfeatures = 2000)

abrain=ScaleData(object = abrain, vars.to.regress = c("percent.mito"))
abrain=RunPCA(object = abrain, features = VariableFeatures(object = abrain), npcs = 100, verbose = FALSE)
abrain=ProjectDim(object = abrain, verbose = FALSE)


npc=20
abrain=FindNeighbors(object = abrain, dims = 1:npc, k.param = 100) 
abrain=FindClusters(object = abrain, resolution = 0.2, graph.name = 'RNA_snn')
abrain=RunUMAP(object = abrain, min.dist = 1.5, graph = 'RNA_snn') 


pdf(paste0("2.seurat.analysis.pdf"))
DimPlot(object = abrain, reduction='pca', cols=color, group.by = 'group2')
ElbowPlot(object = abrain, ndims = 100)
DimPlot(abrain, cols=color)
DimPlot(abrain, cols=color, label=T)
DimPlot(abrain, cols=color, group.by = 'subject_id')
DimPlot(abrain, cols=color, group.by = 'group1')
DimPlot(abrain, cols=color, group.by = 'group2')
DimPlot(abrain, cols=color, group.by = 'tissue')
FeaturePlot(object = abrain, reduction='umap', features = c('nCount_RNA'))
FeaturePlot(object = abrain, reduction='umap', features = c('nFeature_RNA'))
dev.off()


meta.data=abrain@meta.data
saveRDS(abrain, file = "seurat.object.rds")
