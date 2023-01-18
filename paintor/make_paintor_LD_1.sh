########### The purpose of this script series is to generate (temporary bc they are massively large) LD files in Plink to prep input files for fastenloc or paintor prep. This is necessary because both tools require multiple input files which each have the exact same SNPs in the exact same order with no identifying column or row names. If only running a couple loci on each chromosome (i.e., no genes which might both have the same SNP), you could just calculate one file per chromosome but the LD matrices are VERY large so beware it will immediately take up a lot of space. 

#### bcftools github: http://samtools.github.io/bcftools/bcftools.html
#### tabix website: http://www.htslib.org/doc/tabix.html
#### Fastenloc github: https://github.com/xqwen/fastenloc/ 
#### Dap-g github: https://github.com/xqwen/dap
#### Plink 2 website: https://www.cog-genomics.org/plink2

eur_in="/path/to/vcf/files/hg38_EURnoFIN_paintor.vcf.gz"

########################################

module load bcftools tabix 

### First, clear out old stuff from previous chromosomes/runs
        rm /scratch/dir/genes/*
        rm /scratch/dir/*_plink_*

#### Loop set up to run one chromosome at a time because the files (and number of files) can get massive
######## Step 1: make list of SNPs and reference VCF for each eGene (can use your own LD if similar to GWAS/annotations you're using)

a="account=accountname"
p="p queue"
m="mem=40g"
t="t 30:00"

for chr in "22"; do
gene_file="/scratch/dir/chr${chr}_Sep2022_eQTL_siggenes.lst"
chromosome="chr${chr}"

        while read gene; do #this wrapper could be gene-name based or ensg based

        zcat /scratch/dir/Date_qtl_gwas_forpaintor_chr${chr}.txt.gz | grep "${gene}$" | cut -f2 -d ' ' > /scratch/dir/genes/${gene}_snplist.tmp
        keepsnpsfile="/scratch/dir/genes/${gene}_snplist.tmp"
        eur_out="/scratch/dir/dir/${gene}_hg38_EURnoFIN_paintor.vcf.gz"

        sbatch --$a --$m -$t -$p --wrap="bcftools view -Oz -i 'ID=@${keepsnpsfile}' -o $eur_out $eur_in -r '${chromosome}'"

        done < $gene_file
done
