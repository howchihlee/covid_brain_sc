base_url=https://cf.10xgenomics.com/samples/cell-exp/1.1.0/
suffix=_filtered_gene_bc_matrices.tar.gz

for ct in b_cells cd4_t_helper cd14_monocytes cytotoxic_t cd56_nk cd34;
do
    url=$base_url$ct/$ct$suffix
    fn=$ct$suffix
    wget -nc $url
    tar -xvf $fn
    mv filtered_matrices_mex $ct
done