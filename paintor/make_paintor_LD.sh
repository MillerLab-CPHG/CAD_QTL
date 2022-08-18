##### 
# The purpose of this script is to generate (temporary bc they are massively large) LD files in Plink to prep input files 
# for fastenloc or paintor prep. This is necessary because both tools require multiple input files which each have 
# the exact same SNPs in the exact same order with no identifying column or row names. If only running a couple loci 
# on each chromosome (i.e., no genes which might both have the same SNP), you could just calculate one file per 
# chromosome but the LD matrices are VERY large so beware it will immediately take up a lot of space. 
# Suggested implementation is to uncomment and run each step separately by chromosome and delete temporary files in 
# between for space efficiency.
##### 

##### vcftools github: https://vcftools.github.io/man_latest.html 
##### bcftools github: http://samtools.github.io/bcftools/bcftools.html
##### tabix website: http://www.htslib.org/doc/tabix.html
##### Fastenloc github: https://github.com/xqwen/fastenloc/ 
##### Dap-g github: https://github.com/xqwen/dap
##### Plink 2 website: https://www.cog-genomics.org/plink2

bed="1000G_EURnoFIN_plink"

a="account=cphg-millerlab-vip"
p="p standard"

module load vcftools bcftools tabix plink

for chr in {20..22}; do

       while read gene; do #this wrapper could be gene-name based or ensg based

       zcat /scratch/ch2um/Aug22_mixQTL_res.txt.gz | grep "${gene}$" | cut -f3 -d ' ' > /scratch/ch2um/genes/${gene}_snplist.tmp
       zcat /scratch/ch2um/mixqtl_gwas_forpaintor_chr${chr}.txt.gz | grep "${gene}$" | cut -f3 > /scratch/ch2um/genes/${gene}_snplist.tmp
       zcat /scratch/ch2um/Apr22_mixqtl_gwas_forpaintor_chr${chr}.txt.gz | grep "${gene}$" | cut -f3 -d ' ' > /scratch/ch2um/genes/${gene}_snplist.tmp
       eur_in="/project/cphg-millerlab/chani/1000G_vcfs/hg38_EURnoFIN_paintor.vcf.gz"
       keepsnpsfile="/scratch/ch2um/genes/${gene}_snplist.tmp"
       eur_out="/scratch/ch2um/genes/${gene}_hg38_EURnoFIN_paintor"

       m="mem=40g"
       t="t 30:00"
       sbatch --$a --$m -$t -$p --wrap="vcftools --gzvcf $eur_in --recode --out $eur_out --snps $keepsnpsfile --chr chr$chr --maf 0.01 --max-maf 0.99"

##### Step 1A for LD with monomorphic SNPs
       input_dir="/scratch/ch2um/chr${chr}plink"
       tail -n+198 -q /scratch/ch2um/genes/${gene}_hg38_EURnoFIN_paintor.recode.vcf | cut -f3 > $input_dir/new_${gene}_snplist.tmp

##### Step 2
       vcf_eur_out="/scratch/ch2um/genes/${gene}_hg38_EURnoFIN_paintor"
       bed="1000G_plink"
       mv ${vcf_eur_out}.recode.vcf ${vcf_eur_out}.vcf
       bgzip ${vcf_eur_out}.vcf
       tabix -f ${vcf_eur_out}.vcf.gz

       m="mem=50g"
       t="t 15:00"
       sbatch -o ${gene}_eur_bed.blog -n 4 --$a --$m -$t -$p --wrap="plink --make-bed --vcf ${vcf_eur_out}.vcf.gz --out /scratch/ch2um/${gene}_${bed}_EUR"

##### Step 3
       keepsnpsfile="/scratch/ch2um/chr${chr}plink/new_${gene}_snplist.tmp"
       bed="1000G_plink"
       bfile="/scratch/ch2um/${gene}_${bed}"
       outfile="/scratch/ch2um/chr${chr}plink/${gene}"

       m="mem=100g"
       t="t 20:00"

       # Make LD file for EUR minus Finns
       sbatch -o ${gene}_LD.blog -n 2 --$a --$m -$t -$p --wrap="plink --bfile ${bfile}_EUR --r2 square spaces --chr $chr --extract ${keepsnpsfile} --out ${outfile}_EUR"  

##### Step 4
       rm $bfile *tmp $vcf_out
       rm $keepsnpsfile
       mv ${outfile}_EUR.ld ${outfile}.ld

       done < /path/to/folder/chr${chr}_paintor_genenames.lst
done
