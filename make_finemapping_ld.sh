###################### The purpose of this script is to generate (temporary bc they are massively large) LD files in Plink to prep input files for fastenloc or paintor prep. This is necessary because both tools require multiple input files which each have the exact same SNPs in the exact same order with no identifying column or row names. If only running a couple loci on each chromosome (i.e., no genes which might both have the same SNP), you could just calculate one file per chromosome but the LD matrices are VERY large so beware it will immediately take up a lot of space. 
#### Suggested implementation is to uncomment and run each step separately by chromosome and delete temporary files in between for space efficiency.
#### Input start files are mixQTL results combined with GWAS files to identify all overlapping SNPs, a reference population for LD (I used 1000G EUR without Finns), and a list of genes by chromosome.
#### Is this the best way to do this? Almost certainly not but it works.

#### vcftools github: https://vcftools.github.io/man_latest.html 
#### bcftools github: http://samtools.github.io/bcftools/bcftools.html
#### tabix website: http://www.htslib.org/doc/tabix.html
#### Fastenloc github: https://github.com/xqwen/fastenloc/ 
#### Dap-g github: https://github.com/xqwen/dap
#### Plink 2 website: https://www.cog-genomics.org/plink2

######################

a="account=usergroupname"
p="p slurmqueue"

module load vcftools bcftools tabix plink

for chr in "22"; do

        while read gene; do #this wrapper could be gene-name based or ensg based, check which is easier for you.

########### Step 1 ###########
#       zcat /scratch/userid/mixqtl_gwas_forpaintor_chr${chr}.txt.gz | grep "${gene}$" | cut -f3 -d ' ' > /scratch/ch2um/genes/${gene}_snplist.tmp
#       vcf_in="/path/to/your/reference_vcf/hg38_ancestry_population.vcf.gz"
#       keepsnpsfile="/scratch/userid/genes/${gene}_snplist.tmp"
#       vcf_out="/scratch/userid/genes/${gene}_hg38_ancestry_population"

#       m="mem=40g"
#       t="t 2:00:00"
#       sbatch --$a --$m -$t -$p --wrap="vcftools --gzvcf $vcf_in --recode --out $vcf_out --snps $keepsnpsfile --chr chr$chr"

########### Step 2 ###########
#       vcf_out="/scratch/userid/genes/${gene}_hg38_ancestry_population"
#       bed="hg38_ancestry_population"
#       mv ${vcf_out}.recode.vcf ${vcf_out}.vcf
#       bgzip ${vcf_out}.vcf
#       tabix -f ${vcf_out}.vcf.gz

#       m="mem=50g"
#       t="t 15:00"
#       sbatch -o ${gene}_bed.blog -n 4 --$a --$m -$t -$p --wrap="plink --make-bed --vcf ${vcf_out}.vcf.gz --out /scratch/ch2um/${gene}_${bed}"

########### Step 3 ###########
#       bfile="/scratch/userid/${gene}_${bed}"
#       outfile="/scratch/userid/chr${chr}plink/${gene}"

#       m="mem=100g"
#       t="t 20:00"
#       sbatch -o ${gene}_LD.blog -n 2 --$a --$m -$t -$p --wrap="plink --bfile $bfile --r2 square spaces --chr $chr --out $outfile"  

########### Step 4 ###########
#       vcf_out="/scratch/userid/genes/${gene}_hg38_ancestry_population"
#       bfile="/scratch/userid/${gene}_${bed}"
#       keepsnpsfile="/scratch/userid/genes/${gene}_snplist.tmp"
#       rm $bfile $keepsnpsfile $vcf_out

        done < /scratch/userid/chr${chr}_genenames.lst
done
