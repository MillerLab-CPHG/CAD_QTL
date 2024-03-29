########### The purpose of this script series is to generate (temporary bc they are massively large) LD files in Plink to prep input files for fastenloc or paintor prep. This is necessary because both tools require multiple input files which each have the exact same SNPs in the exact same order with no identifying column or row names. If only running a couple loci on each chromosome (i.e., no genes which might both have the same SNP), you could just calculate one file per chromosome but the LD matrices are VERY large so beware it will immediately take up a lot of space. 

#### Plink 2 website: https://www.cog-genomics.org/plink2

bed="1000G_plink"

########################################

module load plink

a="account=accountname"
p="p queue"
m="mem=100g"
t="t 10:00"

for chr in "22"; do

        while read gene; do

### Step 3: make gene-specific LD files, in this case using non-Finnish Europeans from 1000G phase 3
        keepsnpsfile="/scratch/dir/chr${chr}plink/new_${gene}_snplist.tmp"
        bfile="/scratch/dir/${gene}_${bed}"
        outfile="/scratch/dir/chr${chr}plink/${gene}"

        sbatch -o ${gene}_LD.blog -n 2 --$a --$m -$t -$p --wrap="plink --bfile ${bfile}_${pop} --r2 square spaces --chr $chr --extract ${keepsnpsfile} --out ${outfile}"    

        done < /scratch/dir/chr${chr}_Date_eQTL_siggenes.lst
done

