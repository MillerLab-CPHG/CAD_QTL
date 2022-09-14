########### The purpose of this script series is to generate (temporary bc they are massively large) LD files in Plink to prep input files for fastenloc or paintor prep. This is necessary because both tools require multiple input files which each have the exact same SNPs in the exact same order with no identifying column or row names. If only running a couple loci on each chromosome (i.e., no genes which might both have the same SNP), you could just calculate one file per chromosome but the LD matrices are VERY large so beware it will immediately take up a lot of space. 
#### Suggested implementation is to uncomment and run each step separately by chromosome and delete temporary files in between for space efficiency.

#### Plink 2 website: https://www.cog-genomics.org/plink2

bed="1000G_plink"

########################################

module load plink

a="account=cphg-millerlab-vip"
p="p standard"
m="mem=100g"
t="t 20:00"

for chr in "22"; do

        while read gene; do

### Step 3: make gene-specific LD files, in this case using non-Finnish Europeans from 1000G phase 3
        keepsnpsfile="/scratch/ch2um/chr${chr}plink/new_${gene}_snplist.tmp"
        bfile="/scratch/ch2um/${gene}_${bed}"
        outfile="/scratch/ch2um/chr${chr}plink/${gene}"

        sbatch -o ${gene}_LD.blog -n 2 --$a --$m -$t -$p --wrap="plink --bfile ${bfile}_${pop} --r2 square spaces --chr $chr --extract ${keepsnpsfile} --out ${outfile}"    

        done < /scratch/ch2um/chr${chr}_Sep2022_eQTL_siggenes.lst
done

