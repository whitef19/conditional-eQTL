
# Test eQTLs with 3 models for genes located <= 500 kb from lead SNPs 
# - mothers genotypes
# - offspring genotypes
# - mothers genotype conditional for offspring genotypes
# Models are ajustmented for 7 expression PCs, 5 genetic PCs and sex

arg = commandArgs(trailingOnly = T)
expression_file = arg[1] # RNA-seq expression matrix (normalized)
covariates_file = arg[2] # covariates associated to RNA-seq dataset and genetic PCs
mom_genotype_file = arg[3] # mothers genotype doses
offspring_genotype_file = arg[4] # offsprings genotype doses

pairs_list_file = arg[5] # list of pairs genes-SNP to test (genes located inside 500kb from one of the lead SNPs)


expression_df = read.csv(expression_file,sep="\t",row.names=4, check.names=F)[,-c(1,2,3)]
mom_dose_df = data.frame(t(read.csv(mom_genotype_file, sep="\t", check.names=F, row.names=1)))
mom_dose_df = cbind(ID=rownames(mom_dose_df), mom_dose_df)

offspring_dose_df = data.frame(t(read.csv(offspring_genotype_file, sep="\t", check.names=F, row.names=1)))
offspring_dose_df = cbind(ID=rownames(offspring_dose_df), offspring_dose_df)

covariates_df = read.csv(covariates_file, sep="\t", check.names=F)
pairs_df = read.csv(pairs_list_file, sep="\t", header=F)
colnames(pairs_df) = c("geneID", "rsID")

### filter expression df
expression_df = data.frame(t(expression_df[unique(pairs_df$geneID), ]))
expression_df = cbind(ID=rownames(expression_df), expression_df)

df = merge(covariates_df, expression_df, by="ID")

mom_df = merge(df, mom_dose_df, by="ID")
rownames(mom_df) = mom_df$ID

offspring_df = merge(df, offspring_dose_df, by="ID")
rownames(offspring_df) = offspring_df$ID

colnames(mom_dose_df) = c("ID", paste0("mom_", colnames(mom_dose_df)[-1])) 
colnames(offspring_dose_df) = c("ID", paste0("offspring_", colnames(offspring_dose_df)[-1])) 
full_df = merge(df, mom_dose_df, by="ID")
full_df = merge(full_df, offspring_dose_df, by="ID")

report_df <- as.data.frame(matrix(data=NA, nrow=nrow(pairs_df), ncol=10))
colnames(report_df) <- c("rsID","gene","mother_beta", "mother_p", "offspring_beta","offspring_p","cond_mother_b","cond_mother_p","cond_offspring_b", "cond_offspring_p")

for ( i in  1:nrow(pairs_df)){

	geneID = as.character(pairs_df[i, "geneID"])
	rsID = as.character(pairs_df[i, "rsID"])

	if (rsID %in% colnames(mom_df)){
		fit_mom = summary(lm(mom_df[, geneID] ~ mom_df[, rsID] + sex + mom_PC1 + mom_PC2 + mom_PC3 + mom_PC4 + mom_PC5 + ePC1 + ePC2 + ePC3 + ePC4 + ePC5 + ePC6 + ePC7, data = mom_df))
		model_mom_b = signif(fit_mom$coefficients[2,1], 2)
    	model_mom_p = format(signif(fit_mom$coefficients[2,4], 2), scientific = T)
	} else {
		model_mom_b = model_mom_p = "NA"
	}
	
	if (rsID %in% colnames(offspring_df)) {
		fit_offspring = summary(lm(offspring_df[, geneID] ~ offspring_df[, rsID] + sex + off_PC1 + off_PC2 + off_PC3 + off_PC4 + off_PC5 + ePC1 + ePC2 + ePC3 + ePC4 + ePC5 + ePC6 + ePC7, data = offspring_df))
		model_off_b = signif(fit_offspring$coefficients[2,1], 2)
    	model_off_p = format(signif(fit_offspring$coefficients[2,4], 2), scientific = T)
	} else {
		model_off_b = model_off_p = "NA"
	}

	if (model_mom_b != "NA" &  model_off_b != "NA") {
		fit_conditional = summary(lm(full_df[, geneID] ~ full_df[, paste0("mom_", rsID)] + full_df[, paste0("offspring_", rsID)] + sex + mom_PC1 + mom_PC2 + mom_PC3 + mom_PC4 + mom_PC5 + ePC1 + ePC2 + ePC3 + ePC4 + ePC5 + ePC6 + ePC7, data = full_df))
		cond_mom_b = signif(fit_conditional$coefficients[2,1], 2)
	    cond_mom_p = format(signif(fit_conditional$coefficients[2,4], 2), scientific = T)
		
		cond_off_b = signif(fit_conditional$coefficients[3,1], 2)
	    cond_off_p = format(signif(fit_conditional$coefficients[3,4], 2), scientific = T)
	} else {
		cond_mom_b = cond_mom_p = cond_off_b = cond_off_p = "NA"
	}

	report_df[i, ] = c(rsID, geneID, model_mom_b, model_mom_p, model_off_b, model_off_p, cond_mom_b, cond_mom_p, cond_off_b, cond_off_p)
}
write.table(report_df, file="associations.tsv", sep="\t", row.names=F, quote=F)

