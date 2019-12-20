#!/usr/bin/env bash

# The master script to call everything else.
#
# Requirements: a machine with Ubuntu and conda.
#
# Federico Marotta (federico.marotta@edu.unito.it)
# Dec 2019

# Use stricter bash options
set -eo pipefail
IFS=$'\n\t'

# Declare parameters and resources
tissue='Cells_EBV_transformed_lymphocytes'
tissue_name='Cells - EBV-transformed lymphocytes'
epigenome='GM12878'
chr='chr22'
tmpdir='/mnt/red/fmarotta/tmp'
sort_mem='20G'
threads=5

# Data paths
rpkm='data/gtex_expr/GTEx_Data_20160115_v7_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct.gz'
cnts='data/gtex_expr/GTEx_Data_20160115_v7_RNAseq_RNASeQCv1.1.8_gene_reads.gct.gz'
geno="data/gtex_geno/${chr}.biallelic.phased.vcf.gz"
attr='data/gtex_annot/GTEx_v7_Annotations_SampleAttributesDS.txt'
cov="data/gtex_cov/${tissue}.v7.covariates.txt"
regreg="data/regreg/meta/${epigenome}.regulatory_regions.bed"
contigs='data/gtex_ref/gencode.v19.genes.v7.patched_contigs.gtf'
pwm='data/hocomoco/HOCOMOCOv11_full_HUMAN_mono_transfac_format.txt'
bkg='data/background/hg19.background_altern_format.txt'
ref="data/ucsc/${chr}.fa.gz"
tf='data/hocomoco/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv'

# Add the binaries to PATH
export PATH=${PATH}:${PWD}/bin

# Allow us to activate conda environments from within the script
eval "$(conda shell.bash hook)"


###################
### PREPARATION ###
###################

# Master conda environment
conda create -y -n master.FM_MAT0043 -c bioconda -c conda-forge -c r bioconductor-multtest r-data.table r-remotes r-docopt
# Normalizing expression
conda create -y -n normexpr.FM_MAT0043 -c conda-forge -c r -c anaconda pip r-docopt r-data.table
conda activate normexpr.FM_MAT0043
pip install numpy scipy pandas docopt
conda deactivate
# Regulatory regions
conda create -y -n regreg.FM_MAT0043 -c bioconda -c conda-forge -c r bedtools r-data.table r-docopt
# TBA
conda create -y -n tba.FM_MAT0043 -c conda-forge -c r zlib rust r-docopt r-data.table r-iotools
# TBF + TF
conda create -y -n biomart.FM_MAT0043 -c bioconda -c conda-forge -c r bioconductor-biomart r-data.table r-docopt r-xml2
# Models
conda create -y -n models.FM_MAT0043 -c conda-forge -c r r-ranger r-pls r-glmnet r-BART r-iotools r-data.table r-docopt r-metap
# TIGAR
conda create -y -n tigar.FM_MAT0043 -c bioconda -c conda-forge -c r -c anaconda tabix r-data.table r-docopt pip
conda activate tigar.FM_MAT0043
pip install dfply sklearn statsmodels
conda deactivate

# Make the directory structure
mkdir -p results/expression/{rna_processed,adjusted}
mkdir -p results/regreg/gtex/{reference,associations}
mkdir -p results/{tba,tba_tf}/gtex/${chr}
mkdir -p results/models/{ranger,bart,ridge,pcr,grlasso}/{tba,tba_plus_tf}/gtex/${chr}

# Install vcf_rider
conda activate tba.FM_MAT0043
mkdir -p src
cd src/
git clone https://github.com/vodkatad/vcf_rider
cd vcf_rider/
RUSTFLAGS="-L/usr/lib/x86_64-linux-gnu" cargo build --release
ln -s ../src/vcf_rider/target/release/vcf_rider ../../bin/
ln -s ../src/vcf_rider/target/release/indel_stats ../../bin/
cd ../../
conda deactivate

# Install TIGAR
cd results/models/
git clone https://github.com/xmeng34/TIGAR tigar
cd tigar/
find . -name "*.sh" | xargs sed -i 's|#!/usr/bin/bash|#!/usr/bin/env bash|'
find . -name "*.sh" | xargs chmod u+x
find . -name "*.py" | xargs chmod u+x
find . -name "*.py" | xargs sed -i 's|time.clock()|time.perf_counter()|'
sed -i 's|if sum(k_fold_R2)/5 < 0.01:|if False:|' ./Model_Train_Pred/DPR_Train.py
mkdir -p tba/gtex
cd ../../../

# Install additional R packages
conda activate master.FM_MAT0043
mkdir -p lib/R
tar=`which tar`
gzip=`which gzip`
Rscript -e "install.packages('remotes', lib = 'lib/R', repos = 'https://cran.r-project.org')"
Rscript -e ".libPaths(c('lib/R', .libPaths())); Sys.setenv(TAR=\"$tar\"); .Options\$unzip <- \"$gzip\"; remotes::install_github(c('fmarotta/fplyr', 'fmarotta/cvtools', 'fmarotta/byrow'))"
conda deactivate


############################
### NORMALIZE EXPRESSION ###
############################

conda activate normexpr.FM_MAT0043

echo "Normalizing expression"
python scripts/1-normalize_expression.py \
    --rpkm=$rpkm \
    --cnts=$cnts \
    --geno=$geno \
    --attr=$attr \
    --out=results/expression/rna_processed

echo "Obtaining the residuals of expression ~ covariates"
Rscript scripts/2-adjust_expression.R \
    --qnorm_std=results/expression/rna_processed/${tissue}.normalized_expression.qnorm.std.tsv \
    --covariates=$cov \
    --adj_expr=results/expression/adjusted/${tissue}.expression_residuals.tsv


##########################
### REGULATORY REGIONS ###
##########################

conda activate regreg.FM_MAT0043

echo "Getting the genes annotation"
sed '/^##/d' $contigs | \
    grep -f data/meta/genes.txt | \
    bawk '$3 == "gene" {split($9, a, "\""); print "chr"$1, $4, $5, a[2], $6, $7}' | \
    sort -k1,1V -k2,2n > results/regreg/gtex/reference/genes_annotations.bed

echo "Finding enhancers associations"
bedtools intersect -wo \
                   -f 1 \
                   -a <(bawk '$5 == 1' $regreg) \
                   -b <(bawk '{$2 = ($2 > 1000000) ? $2 - 1000000 : 0; $3 = $3 + 1000000; print}' results/regreg/gtex/reference/genes_annotations.bed) | \
    bawk '{print $1, $9, $4}' | \
    sort > results/regreg/gtex/associations/${epigenome}.enhancers_associations.tsv

echo "Finding promoters associations"
bedtools closest -id \
                 -D a \
                 -a results/regreg/gtex/reference/genes_annotations.bed \
                 -b <(bawk '$5 == 2 || $5 == 3' $regreg) | \
    bawk '$11 != -1' | \
    bawk '$12 > -5000' | \
    bawk '$12 != 0 || ($6 == "-" && $3 - $8 < 1000000 || $6 == "+" && $9 - $2 < 1000000)' | \
    bawk '{print $1, $4, $10}' | \
    sort > results/regreg/gtex/associations/${epigenome}.promoters_associations.tsv

echo "Finding default associations"
bawk '($5 == 4) {n = split($4, a, "_"); print $1, a[2], $4}' $regreg | \
    grep -f data/meta/genes.txt | \
    sort > results/regreg/gtex/associations/${epigenome}.default_associations.tsv

echo "Getting the list of associated regions"
cat results/regreg/gtex/associations/${epigenome}.{promoters,enhancers,default}_associations.tsv | \
    cut -f 3 | sort -u > results/regreg/gtex/associations/${epigenome}.associated_regions.txt


###########
### TBA ###
###########

conda activate tba.FM_MAT0043

echo "Running indel_stats"
indel_stats \
    <(vcftools --gzvcf $geno \
               --keep <(cut -f 1,7 $attr | grep "$tissue_name" | cut -d - -f 1,2) \
               --recode \
               --stdout) \
    <(grep -f results/regreg/gtex/associations/${epigenome}.associated_regions.txt $regreg) \
    > results/tba/gtex/${chr}/${epigenome}.${tissue}.indel_stats.tsv

echo "Computing the TBA"
vcf_rider \
    -b <(grep -f <(bawk '$6 == "false" {print $1}' results/tba/gtex/${chr}/${epigenome}.${tissue}.indel_stats.tsv) $regreg) \
    -p <(transfac2pcm < $pwm) \
    -f $bkg \
    -v <(vcftools --gzvcf $geno \
                  --keep <(cut -f 1,7 $attr | grep "$tissue_name" | cut -d - -f 1,2) \
                  --recode \
                  --stdout) \
    -r <(zcat $ref | sed 's/chr//1' | bawk '{print toupper($0)}') \
    -a results/tba/gtex/${chr}/${epigenome}.${tissue}.vcf_rider_snps.tsv \
    > results/tba/gtex/${chr}/${epigenome}.${tissue}.vcf_rider_tba.tsv

echo "Joining the genes and the reg. reg."
cat results/regreg/gtex/associations/${epigenome}.{promoters,enhancers,default}_associations.tsv | \
    cut -f 2-3 | \
    sort -k2,2 | \
    join -t$'\t' \
         -1 2 \
         -2 1 \
         -o 1.1,1.2,2.2,2.3,2.4,2.5,2.6,2.7 \
         - <(sort -t $'\t' -S $sort_mem -T $tmpdir --parallel=2 -k1,1 results/tba/gtex/${chr}/${epigenome}.${tissue}.vcf_rider_tba.tsv) | \
    sort -t $'\t' -S $sort_mem -T $tmpdir --parallel=2 -k1,2 > results/tba/gtex/${chr}/${epigenome}.${tissue}.sorted_tba.tsv

echo "Aggregating the TBA"
Rscript scripts/3-aggregate_tba.R \
    --snps=results/tba/gtex/${chr}/${epigenome}.${tissue}.vcf_rider_snps.tsv \
    --sorted_tba=results/tba/gtex/${chr}/${epigenome}.${tissue}.sorted_tba.tsv \
    --tba=results/tba/gtex/${chr}/${epigenome}.${tissue}.tba.tsv \
    --fraction=0.1 \
    --threads=$threads

conda activate biomart.FM_MAT0043

echo "Convoluting the TBA"
Rscript scripts/4-convolute_tba.R \
    --tba=results/tba/gtex/${chr}/${epigenome}.${tissue}.tba.tsv \
    --tf=$tf \
    --expr=results/expression/adjusted/${tissue}.expression_residuals.tsv \
    --tba_tf=results/tba_tf/gtex/${chr}/${epigenome}.${tissue}.tba_plus_tf.tsv \
    --input_function='plus' \
    --threads=$threads


##############
### MODELS ###
##############

conda activate models.FM_MAT0043

echo "Running models (tba)"
for m in "RIDGE" "BART" "PCR" "RANGER"; do
    echo "Now running $m"
    Rscript scripts/5-bygene_perf.R \
        --tba=results/tba/gtex/${chr}/${epigenome}.${tissue}.tba.tsv \
        --expr=results/expression/adjusted/${tissue}.expression_residuals.tsv \
        --genes=data/meta/genes.txt \
        --samples=data/meta/samples.txt \
        --model=${m}_model \
        --autoparams=${m}_autoparams \
        --threads=$threads \
        --out_prefix=results/models/${m,,}/tba/gtex/${chr}/${epigenome}.${tissue}
done

echo "Running models (tba plus tf)"
for m in "RIDGE" "BART"; do
    echo "Now running $m"
    Rscript scripts/5-bygene_perf.R \
        --tba=results/tba_tf/gtex/${chr}/${epigenome}.${tissue}.tba_plus_tf.tsv \
        --expr=results/expression/adjusted/${tissue}.expression_residuals.tsv \
        --genes=data/meta/genes.txt \
        --samples=data/meta/samples.txt \
        --model=${m}_model \
        --autoparams=${m}_autoparams \
        --reshape=reshape_sum \
        --threads=$threads \
        --out_prefix=results/models/${m,,}/tba_plus_tf/gtex/${chr}/${epigenome}.${tissue}
done

#############
### TIGAR ###
#############

conda activate tigar.FM_MAT0043

echo "Obtaining the expression in TIGAR format"
Rscript scripts/6-tigar_expr.R \
    --annot=$contigs \
    --expr=results/expression/adjusted/${tissue}.expression_residuals.tsv \
    --tigexpr=results/expression/adjusted/${tissue}.tigar_expression.tsv

cd results/models/tigar
echo "Running TIGAR"
./TIGAR_Model_Train.sh \
    --model DPR \
    --Gene_Exp <(bawk 'NR == 1; NR > 1 {print $0 | "grep -f ../../../data/meta/genes.txt"}' ../../expression/adjusted/${tissue}.tigar_expression.tsv) \
    --train_sample ../../../data/meta/samples.txt \
    --chr 22 \
    --genofile_dir ../../../$geno \
    --genofile_type vcf \
    --Format GT \
    --thread $threads \
    --out ${PWD}/tba/gtex

