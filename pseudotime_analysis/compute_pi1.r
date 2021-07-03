library('limma')
df = read.csv('./pval_cluster.csv')
for (i in unique(df$cluster)){
    
    print(paste(i, 1-propTrueNull(df$pval[df$cluster == i]) ))}