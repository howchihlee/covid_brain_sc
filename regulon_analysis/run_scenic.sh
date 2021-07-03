conda activate scenic_protocol
arboreto_with_multiprocessing.py input/gene_expression.loom resources/hs_hgnc_curated_tfs.txt \
    --method grnboost2 \
    --output output/adj.tsv \
    --num_workers 30 \
    --seed 777
        
pyscenic ctx output/adj.tsv \
    ./cisTargetdb/hg38/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather ./cisTargetdb/hg38/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather \
    --annotations_fname ./cisTargetdb/hg38/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname input/gene_expression.loom \
    --output output/reg.csv \
    --mask_dropouts \
    --num_workers 30
        
pyscenic aucell input/gene_expression.loom \
    output/reg.csv \
    --output output/pyscenic_aucell_output.loom \
    --num_workers 30
    