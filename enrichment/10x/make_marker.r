library(Seurat)

write_csv = function(vec, fn){write.table(vec, file = fn, append = FALSE, quote = FALSE, sep = ",",
                eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                col.names = TRUE, qmethod = c("escape", "double"),
                fileEncoding = "")}

data = Read10X(data.dir = "b_cells/hg19")
x <- CreateSeuratObject(counts = data,  project = 'b_cells')

to_merge_name = c('b_cells', 'cd4_t_helper', 'cd14_monocytes', 'cytotoxic_t', 'cd56_nk', 'cd34')
to_merge_obj = c()

for (ct in to_merge_name[2:6]) {
    fn = paste0(ct, '/hg19')
    data = Read10X(data.dir = fn)
    obj <- CreateSeuratObject(counts = data, project = ct)
    to_merge_obj = cbind(to_merge_obj, obj)
    }

pbmc.big <- merge(x, y = to_merge_obj, add.cell.ids = to_merge_name, project = "PBMCbig")
pbmc.big = NormalizeData(pbmc.big, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc.big = ScaleData(pbmc.big)

markers = FindAllMarkers(pbmc.big)
write_csv(markers, 'seurat_markers.csv')