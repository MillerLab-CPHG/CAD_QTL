########### The purpose of this script is to run Summary-data-based 'Mendelian Randomization' (SMR) on eQTL data 
### "to test for pleiotropic association between the expression level of a gene and a complex trait of interest using 
### summary-level data from GWAS and expression quantitative trait loci (eQTL) studies". 
### Everything is formatted for slurm controller job submission.
### More info on this command-line program can be found here: https://yanglab.westlake.edu.cn/software/smr/
### Original citation (Zhihong Zhu, et al Nat Gen 2016) here: https://www.nature.com/articles/ng.3538/ 

###############################

SMR="/path/to/SMRprogram/smr_Linux"

############# STEP 1 ################## Make vcf files for LD calculation using VCFtools (https://vcftools.github.io/man_latest.html) 
### and Tabix (http://www.htslib.org/doc/tabix.html)

       module load vcftools tabix

       vcf_in="/path/to/your/referenceVCF/vcf_filename.vcf.gz" #Total reference population you're using--all 1000G for example
       keepfile="/path/to/your/referenceVCF/Ancestry_population_IDs.lst" #List of IDs for individuals to keep
       keepsnpsfile="/path/to/your/SMRfolder/mixqtl_plink_snplist.txt" #List of variants present in your results file
       output="/path/to/your/SMRfolder/hg38_ancestry_filename"
       
       a="account=slurmgroup"
       p="p slurmqueue"
       m="mem=25g"

       sbatch --$a --$m -$p -t 10:00:00 --wrap="vcftools --gzvcf $vcf_in --out h38_EURnoFIN_paintor --recode --keep $keepfile --snps $keepsnpsfile"
       bgzip $output.recode.vcf
       tabix $output.recode.vcf.gz 

       sbatch --$a --$m -$p -t 2:00:00 --wrap="vcftools --gzvcf vcf.gz --out h38_EURnoFIN_paintor --freq --snps snps_for_SMR.txt"
       bgzip $output.frq
 
############# STEP 2 ################## Make bed files for LD calculation using Plink v1.9: https://www.cog-genomics.org/plink/
       module load plink
       a="account=slurmgroup"
       p="p slurmqueue"
       m="mem=80g"
       t="t 60:00"
       vcf="/path/to/vcf_file/hg38_ancestry_purpose.vcf.gz"
       sbatch -o 1000G_EUR_bed.blog -n 4 --$a --$m -$t -$p --wrap="plink --make-bed --vcf $vcf --out hg38_ancestry_plink"


############# STEP 3 ################# mixQTL rseults files formatted like FastQTL cis nominal pass results:
### This space-delimited file is formatted to match qtltools "nominal cis pass" output. NO HEADER.
### The columns are described here: https://qtltools.github.io/qtltools/ 

       qtl="/path/to/SMRfolder/mixQTL_qtltools_SMRinput.txt.gz"
       mybesd="/path/to/SMRfolder/mixQTL_SMR"
       a="account=slurmgroup"
       p="p slurmqueue"
       m="mem=30g"
       t="t 10:00"

       sbatch --$a --$m -$p -$t --wrap="$SMR --eqtl-summary $qtl --qtltools-nominal-format --make-besd --out $mybesd"

############# STEP 3A ################# Once besd files made, if using QTLtools format need to combine with 
### allele/freq info generated in second job from step 1. This was done separately in R.
### NOTE: SNP ORDER MUST NOT BE CHANGED


############# STEP 4 ################## Job submission to slurm controller
### GWAS summary stats need to be formatted like COJO file, space-delimited
### SNP A1 A2 freq b se p N 
### rs1001 A G 0.8493 0.0024 0.0055 0.6653 129850

        a="account=cphg-millerlab-vip"
        p="p standard"
        m="mem=80g"
        t="t 60:00"

        mydata="hg38_ancestry_plink"
        qtl="/path/to/SMRfolder/mixQTL_SMR"
        GWAS1="/path/to/SMRfolder/GWAS1_sumstats.txt"
        out_GWAS1="/path/to/SMRfolder/SMR_eQTL_GWAS1"

        sbatch --$a --$m -$p -$t --wrap="$SMR --bfile $mydata --gwas-summary $GWAS1 --beqtl-summary $qtl --out $out_GWAS1 --thread-num 10 --maf 0.01" # --extract-snp ${gene}_snps.lst 

