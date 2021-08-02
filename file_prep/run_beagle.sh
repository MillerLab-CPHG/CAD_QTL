#!bin/bash

########################################
## The purpose of this script is to perform phasing and imputation on VCF files 
## obtained from gencove's beta version low-pass whole genome sequencing data.
## 
## NOTE: running this on ORIGINAL B37 files before liftover. See Picard liftover script for more  
## 
## More info on Beagle can be found here: https://bioinformaticshome.com/tools/imputation/descriptions/BEAGLE.html 
## Beagle citation here:  https://doi.org/10.1016/j.ajhg.2009.01.005    PMID: 19200528
########################################

program="/path/to/program/pkg/beagle.27Apr20.b81.jar"

refdir="/path/to/refgenome"
wd="/path/to/your/vcfs"

vcf="vcf_filename.vcf.gz" # If not merged, can call separate job for each chromosome since you only run one at a time anyway

a="account-name"
t="2-12:00:00"
p="standard"

##################################

module load java/1.8.0

cd ${wd}

for c in {1..22}; do

sbatch -t ${t} -p $p --account=$a --wrap="java -Xmx20g -jar ${program} gt=${wd}/${vcf} ref=${refdir}/chr${c}.1kg.phase3.v5a.b37.bref3 out=fileinfo_phase_imp_b37_chr${c} chrom=${c} impute=true gp=true seed=149827"

done
