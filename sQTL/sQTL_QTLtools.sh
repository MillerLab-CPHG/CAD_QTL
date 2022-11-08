#!/bin/bash

#### UVA-CPHG Miller Lab: The purpose of this script is to submit QTLtools jobs by chromosome to a slurm controller. 
#### The two subsets of scripts run permutation- and nominal-pass, respectively. Inputs are genotype (vcf), LeafCutter output (bed), and covariates.
#### BCFtools to make MAF-filtered VCFs can be found here: https://samtools.github.io/bcftools/bcftools.html
#### Leafcutter can be found here: https://davidaknowles.github.io/leafcutter/
#### QTLtools can be found here: https://qtltools.github.io/qtltools/

run_dir="/path/to/file/LeafCutter"

vcf_dir="/path/to/file/vcf/hg38"
bed="/path/to/file/bedfile/Tissue_type.leafcutter.modified.bed.gz"
covs="/path/to/file/covariates/Tissue_type.combined_covariates.txt"
nominal_covs="/path/to/file/nominal_analyses/covariates.txt"

#Load qtltools module
module load gcc/7.1.0 qtltools/1.3.1

for j in $(seq 1 22); do
       sbatch -A cphg-millerlab-vip -p standard -n 2 -t 1-12:00:00 --mem=25g \
               --wrap="QTLtools cis --vcf $vcf_dir/chr${j}_vcf_file_name_hg38.vcf.gz \
                       --bed $bed \
                       --permute 100000 --window 250000 --std-err \
                       --cov $run_dir/$covs --include-covariates $run_dir/$nominal_covs \
                       --out $run_dir/Tissue_type_chr${j}_sqtls_perm100k.txt.gz \
                       --region chr${j}"
done

for j in $(seq 1 22); do
        sbatch -A cphg-millerlab-vip -p standard -t 10:00:00 --mem=25g \
                --wrap="QTLtools cis --vcf $vcf_dir/chr${j}_vcf_file_name_hg38.vcf.gz --bed $bed \
                --nominal 1 --window 250000 --std-err \
                --cov $run_dir/$covs --include-covariates $run_dir/$nominal_covs\
                --out $run_dir/Tissue_type_chr${j}_sqtls_nominal.txt.gz \
                --region chr${j}"
done
