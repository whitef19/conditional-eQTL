library(ggplot2)
library(corrplot)

arg = commandArgs(trailingOnly = T)

expression_file = arg[1] # RNA-seq expression matrix (normalized)
phenotype_file = arg[2] # traits and covariates associated to RNA-seq dataset
moms_genetic_pc_file = arg[3] # Genetic PCs 
offspring_genetic_pc_file = arg[4] # Genetic PCs 

outdir = arg[5] # output directory

nb_included_epc = as.numeric(arg[6]) # Number of expression PCs to include
nb_included_mom_gpc = as.numeric(arg[7]) # Number of genetic PCs to include
nb_included_offspring_gpc = as.numeric(arg[8]) # Number of genetic PCs to include

### create new directory for output files
dir.create(outdir)

### open expression matrix
df = read.csv(expression_file, sep="\t", row.names=4, check.names=F)[,-c(1,2,3)]

### run PCA with base R function "pcrcomp"
pca = prcomp(t(df))

### look variance explained by 20 first PCs
n_pc = 20
var = summary(pca)$importance[2,1:n_pc]
Var = data.frame("x"=1:length(var), "var"=as.vector(var)*100)
pdf(paste0(outdir,"/scree_plot.pdf"), width=3, height=3)
ggplot(Var, aes(x=x, y=var)) + geom_bar(stat="identity", fill="orangered2") + ggtitle("Scree plot") + xlab("PCs") + ylab("Variance (%)")
dev.off()

### extract and write 20 first PCs
PCs = data.frame(pca$x[, 1:n_pc])
PCs = cbind(ID=rownames(PCs),PCs)
write.table(PCs, paste0(outdir, "/expression_PCs.tsv"), sep="\t", row.names=F, quote=F)

### plot PCs
pdf(paste0(outdir,"/PC_plots.pdf"), width=3, height=3)
ggplot(PCs, aes(x=PC1, y=PC2)) + geom_point(color="mediumvioletred", alpha=0.5) + ggtitle("PCA") + xlab(paste0("PC1 (", Var$var[1] ," %)")) + ylab(paste0("PC2 (", Var$var[2] ," %)")) + theme(legend.position="none")
ggplot(PCs, aes(x=PC1, y=PC3)) + geom_point(color="mediumvioletred", alpha=0.5) + ggtitle("PCA") + xlab(paste0("PC1 (", Var$var[1] ," %)")) + ylab(paste0("PC3 (", Var$var[3] ," %)")) + theme(legend.position="none")
ggplot(PCs, aes(x=PC2, y=PC3)) + geom_point(color="mediumvioletred", alpha=0.5) + ggtitle("PCA") + xlab(paste0("PC2 (", Var$var[2] ," %)")) + ylab(paste0("PC3 (", Var$var[3] ," %)")) + theme(legend.position="none")
dev.off()

### look for correlation with known variables
cov_df = read.csv(phenotype_file, sep="\t")
merged_df = merge(PCs, cov_df, by="ID")
M = cor(na.omit(merged_df[,-1]))
pdf(paste0(outdir, "/corrplot.pdf"))
corrplot(M[(n_pc+1):ncol(M), 1:n_pc], tl.col="black", col=colorRampPalette(c("orangered","white","darkmagenta"))(200))
dev.off()


### re-plot PCs with coloration according to known variables (example SCT: proportion of syncytiotrophoblast)
pdf(paste0(outdir, "/colored_PC_plots.pdf"), width=3, height=3)
ggplot(merged_df, aes(x=PC1, y=PC2, color=SCT)) + geom_point(alpha=0.5) + xlab(paste0("PC1 (", Var$var[1] ," %)")) + ylab(paste0("PC2 (", Var$var[2] ," %)")) + theme(legend.position="bottom")+scale_colour_gradient(low="gold1",high="darkmagenta")
dev.off()

### we could plot the distribution of the PCs' loadings to get further insight on the dataset
#loadings = data.frame(pca$rotation)[,1:n_pc]
#loadings = cbind(ID=rownames(loadings),loadings)
#pdf(paste0(outdir, "/loadings.pdf"), width=3, height=3)
#ggplot(loadings, aes(x=PC1)) + geom_histogram(fill="darkmagenta")+ylab("Features")+xlab("Loadings PC1")
#ggplot(loadings, aes(x=PC2)) + geom_histogram(fill="darkmagenta")+ylab("Features")+xlab("Loadings PC2")
#dev.off()


#### add genetic PCs to expression PCs and other variables of interest (sex)
mom_genetic_pc = read.csv(moms_genetic_pc_file, sep="\t")
mom_genetic_pc = mom_genetic_pc[, c("ID", paste0("PC", 1:nb_included_mom_gpc))]
colnames(mom_genetic_pc) = c("ID", paste0("mom_PC",1:nb_included_mom_gpc))

off_genetic_pc = read.csv(offspring_genetic_pc_file, sep="\t")
off_genetic_pc = off_genetic_pc[, c("ID", paste0("PC", 1:nb_included_offspring_gpc))]
colnames(off_genetic_pc) = c("ID", paste0("off_PC",1:nb_included_offspring_gpc))

merged_df = merged_df[, c("ID","sex", paste0("PC",1:nb_included_epc))]
colnames(merged_df) = c("ID","sex", paste0("ePC",1:nb_included_epc))

df = merge(merged_df, mom_genetic_pc, by="ID", all.x=T)
df = merge(df, off_genetic_pc, by="ID", all.x=T)

write.table(df, paste0(outdir, "/all_covariates.tsv"), sep="\t", row.names=F, quote=F)

