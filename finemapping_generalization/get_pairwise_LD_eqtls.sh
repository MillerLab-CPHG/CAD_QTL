#!bin/bash

##### The purpose of this script is to generate pairwise LD (in this case in 1KG Europeans only) to better evaluate generalization between lead eQTLs for pairs of datasets. 
##### BCFtools can be found here: https://samtools.github.io/bcftools/bcftools.html
##### Plink can be found here: https://www.cog-genomics.org/plink
##### The 1000 Genomes data portal can be found here: https://www.internationalgenome.org/data-portal/sample


##### Step 1 make VCF of EUR (without FIN) from 1000G reference file including all variants 
keepsnpsfile="snps_for_generalization_ld.txt" #This file is a single column with all SNPs to be included
vcf_in="/path/to/1000G_Phase3/vcfs/1000G_hg38_phased.vcf.gz"
idlist="/path/to/1000G_Phase3/EUR_noFIN_IDs.lst"
vcf_out="/path/to/my/data/1000G3_EURnoFIN_for_generalization.vcf.gz"

module load bcftools tabix

a="account=account_name"
p="p queue_name"
m="mem=25g"
t="t 12:00:00"

##### Job submission to the slurm controller
sbatch --$a --$m -$t -$p --wrap="bcftools view -Oz -i 'ID=@${keepsnpsfile}' -o $vcf_out -S $idlist $vcf_in"


##### Step 2: run Plink to generate pairwise correlation (r2) for each pair (depends on above job being completed)

outdir="/path/to/my/data/LD_calcs"

tissues=("starnetaor" "starnetmam" "gtexaor" "gtexcor" "gtextib")

module load plink

for i in ${tissues[@]}; do
        while read -r gene var1 var2; do

        plink --ld $var1 $var2 --vcf $vcf_out --out $outdir/${var1}_${var2}
        grep "R-sq" $outdir/${var1}_${var2}.log > rsq.tmp
        sed "s/R-sq/${gene} $var1 $var2 $i  R-sq/" rsq.tmp >> $outdir/LD_all.txt #This is tab-delimited (generated using ctrl-v-tab on a Mac)
        rm $outdir/${var1}_${var2}.nosex

        done < mixqtl_${i}_egene_snps.txt #This file has three columns: gene name, my lead eQTL, other study's tissue-specific lead eQTL
done
