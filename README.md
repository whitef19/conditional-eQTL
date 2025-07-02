
# eQTL analysis using placental RNA-seq data from Gen3G cohort

Steps:
1) Prepare RNA-seq data
2) Compute expression PCs
3) Extract test list 
4) Test associations

### 1) Prepare RNA-seq data

Input files:
*placenta_rnaseq.v1.1.clean.maternal_side.gene_reads.gct* and *placenta_rnaseq.v1.1.clean.maternal_side.gene_tpm.gct*

Normalize RNA-seq data following GTEx V8 steps https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl
```bash
# env
source /project/def-jacquesp/whitef/venv/mirqtl/bin/activate
module restore qtl
unset R_LIBS

# variables


rnaseq_samples=${RNASEQ_DIR}/scripts_and_refs/maternal_side_samples.txt
mom_genotype_samples=${MOM_GENOTYPE_DIR}/scripts_and_refs/sample_european_unrelated_list.txt
collapsed_annotation_gtf=${ANNOTATION_PATH}/gencode.v30.annotation.collapsed.gtf
prefix=data/placenta_rnaseq.normalized
included_samples=data/sample_list_maternal_side_eur.txt

# samples included
cut -f1 ${rnaseq_samples} | cat - ${mom_genotype_samples} | sort | uniq -d > ${included_samples}

# variables
python3 eqtl_prepare_expression.py ${RNASEQ_DIR}/placenta_rnaseq.v1.1.clean.maternal_side.gene_tpm.gct ${RNASEQ_DIR}/placenta_rnaseq.v1.1.clean.maternal_side.gene_reads.gct ${collapsed_annotation_gtf} ${rnaseq_samples} ${prefix} \
	--tpm_threshold 0.1 \
	--count_threshold 6 \
	--sample_frac_threshold 0.2 \
	--normalization_method tmm \
	--sample_ids ${included_samples}

```

### 2) Compute expression PCs and add genetic PCs into one single file

```bash
Rscript run_pca.R data/placenta_rnaseq.normalized.expression.bed.gz data/known_variables.tsv ${MOM_GENOTYPE_DIR}/mothers.genetic_pcs.tsv ${OFFSPRING_GENOTYPE_DIR}/offsprings.genetic_pcs.tsv data/ 7 5 5 
```

### 3) Extract test list 

List of lead SNPs extracted from GenDip2 GDrive 

```bash
### extract lead SNP list
awk '{print $1"\t"FILENAME}' data/*leadSNP_v2.txt | sed 's/_/\t/' | cut -f1,2 | sed 's$data/$$g' | grep -v "Marker" | tee data/all_leadSNP_v2.tsv | cut -f1 | sort | uniq > data/all_leadSNP_v2.txt

ml plink/1.9b_6.21-x86_64

### extract mother's genotypes for lead SNPs
plink --bfile ${MOM_GENOTYPE_DIR}/plink/mothers_rsid_turnkey --extract data/all_leadSNP_v2.txt --recode A-transpose --out data/moms_genotypes
cut -f2,7- data/moms_genotypes.traw | sed -E 's/_[0-9]*//g' > data/moms_genotypes.dose.tsv

### extract offspring's genotypes for lead SNPs
plink --bfile ${OFFSPRING_GENOTYPE_DIR}/plink/offsprings_rsid_allvar --extract data/all_leadSNP_v2.txt --recode A-transpose --out data/offspring_genotypes
cut -f2,7- data/offspring_genotypes.traw | sed -E 's/_[0-9]*//g' > data/offspring_genotypes.dose.tsv

### extract SNP positions and create +/- 500 kb windows
awk '{OFS="\t"; print $1,$4,$2}' data/offspring_genotypes.traw > tmp
awk '{OFS="\t"; print $1,$4,$2}' data/moms_genotypes.traw >> tmp
sort tmp | uniq | sort -Vk1,2 | grep -v "CHR" | awk '{OFS="\t"; print "chr"$1, $2-500000, $2+500000,$3,$4}' > data/all_leadSNP_v2.flanked500kb.bed

### extract genes inside 500 kb window
bedtools intersect -wo -a data/all_leadSNP_v2.flanked500kb.bed -b data/placenta_rnaseq.normalized.expression.bed.gz | awk '{OFS="\t"; print $8,$4}' | sort | uniq > data/gene_SNP_pairs_list.tsv
```


### 4) Test associations

```bash

Rscript run_conditional_associations.R data/placenta_rnaseq.normalized.expression.bed.gz data/all_covariates.tsv data/moms_genotypes.dose.tsv data/offspring_genotypes.dose.tsv data/gene_SNP_pairs_list.tsv

```
