########### The purpose of this script series is to generate (temporary bc they are massively large) LD files in Plink to prep input files for fastenloc or paintor prep. This is necessary because both tools require multiple input files which each have the exact same SNPs in the exact same order with no identifying column or row names. If only running a couple loci on each chromosome (i.e., no genes which might both have the same SNP), you could just calculate one file per chromosome but the LD matrices are VERY large so beware it will immediately take up a lot of space. 

#### tabix website: http://www.htslib.org/doc/tabix.html
#### Paintor website: https://github.com/gkichaev/PAINTOR_V3.0
#### Plink 2 website: https://www.cog-genomics.org/plink2

module load tabix plink

bed="1000G_plink"

a="account=cphg-millerlab-vip"
p="p standard"
m="mem=50g"
t="t 15:00"

for chr in "22"; do
input_dir="/scratch/dir/chr${chr}plink"

        while read gene; do

######## Step 2A: remove Euro monomorphic SNPs because it breaks Paintor (files cannot have diff number of rows or infinite/missing values)
        vcf_eur_out="/scratch/ch2um/genes/${gene}_hg38_EURnoFIN_paintor.vcf.gz"
        tabix -f ${vcf_eur_out}
        zcat $vcf_eur_out | tail -n+199 -q | cut -f3 > $input_dir/new_${gene}_snplist.tmp

######## Step 2B: make gene-specific bed files

        sbatch -o ${gene}_eur_bed.blog -n 4 --$a --$m -$t -$p --wrap="plink --make-bed --vcf ${vcf_eur_out} --out /scratch/dir/${gene}_${bed}_EUR"

        done < /scratch/dir/chr${chr}_Date_eQTL_siggenes.lst
done

