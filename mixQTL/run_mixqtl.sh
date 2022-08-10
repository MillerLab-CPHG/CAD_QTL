#!/bin/bash

########## All the scripts contained in this file were used either unchanged or very modestly adapted from the following mixQTL github: https://github.com/hakyimlab/mixqtl
########## Data formats and examples for each script can be found here: https://github.com/hakyimlab/mixqtl/wiki/Example-and-tutorial-of-command-line-tool

###pathnames for scripts used below
parse_script="/path/to/libsize/file/scripts/mixQTL/parse_vcf.py"
gencode_script="/path/to/libsize/file/scripts/mixQTL/parse_gencode.py"
exp_parse_script="/path/to/libsize/file/scripts/mixQTL/parse_expression.py"

###pathnames for directories and files called by aforementioned scripts 
gencode_dir="/path/to/reference/file/gencode"
v37="gencode.v37.annotation.gtf.gz"

mixqtldir="/path/to/file/mixQTL"
vcfdir="/path/to/vcf/file/chr15_hg38_noindels.vcf.gz"
input_dir="/path/to/input/files/mixqtl"

############################# running prep scripts from mixQTL github

module load gcc/7.1.0 python/3.6.8

##### Parse VCF file into format mixQTL likes, finished 4/13/21
        #sbatch --mem=12g -p standard -t 4:00:00 --account=cphg-millerlab --wrap="python $parse_script -vcf ${vcfdir} -output_prefix $hapdir/converted"
        #python $parse_script -vcf ${vcfdir} -output_prefix $hapdir/converted

##### Parse gencode annotation file, finished 4/13/21
#       python $gencode_script -gencode ${gencode_dir}/${v37} -output ${input_dir}/v37

##### Parse raw counts expression data, didn't run this step for our counts bc it was already in the right format. 
#       python $exp_parse_script -expression $raw_cts -output ${hapdir}/exp_hg38_mixQTL.txt.gz   

var_annot="_variant_annotation.txt.gz"
covs="/path/to/covar/file/QTL_covariate_file_March2021.txt"
lib_size="/path/to/libsize/file/RNAseq_library_size.txt" #You need to make this file yourself using log files from RNA QC process
exp_annot="/path/to/expression/TPM/v37_gencode_mixqtl_input_filtered.txt"
exp_reads="/path/to/expression/rawreads/raw_readcounts_filtered_merged.txt.gz"

########## Run mixQTL

module load gcc/7.1.0 openmpi/3.1.4 R/4.0.0

a="groupname"
p="queue"
t="4:00:00"
m="50g"

cd /path/to/output/files/mixQTL

for i in "22"; do

hap1="converted_chr${i}_hap1.txt.gz"
hap2="converted_chr${i}_hap2.txt.gz"

        sbatch --account=$a -p $p -t $t --mem=$m --wrap="Rscript $mixqtldir/run_mixqtl.R -library_size $input_dir/$lib_size -variant_annotation $input_dir/converted_chr${i}$var_annot -window 500000 -haplotype_1 $input_dir/$hap1 -haplotype_2 $input_dir/$hap2 -covariates $covs -expression_annotation $exp_annot -expression_total_count $input_dir/../$exp_reads -expression_count_1 $input_dir/haps_expression_one.tsv.gz -expression_count_2 $input_dir/haps_expression_two.tsv.gz -output mixQTL_res_chr${i}_date.tsv.gz -trc_cutoff 10"

done
