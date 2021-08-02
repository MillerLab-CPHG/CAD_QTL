###### The purpose of this script is to run nominal and permutation eQTL analyses genome-wide using 
### QTLtools, a more efficient implementation of the MatrixQTL framework [CITE: 
### "QTLtools: a complete tool set for molecular QTL discovery and analysis. Nat Commun 8, 15452 (2017)."]
### Download/read about QTLtools here: https://qtltools.github.io/qtltools/
### Non-slurm-based system may require different parameters

#!/bin/bash


############Step 3 run nominal and then permutation passes###################

wd="/path/to/working/directory"
vcfdir="/path/to/your/vcfs"

genotypes="input_filename.vcf.gz"
phenotypes="Prefix_for_expression_bedfile"
prefix="Prefix_for_output_files"
cov_file="Name_of_covariatefile_Date.txt"

## Note that chromosomes with larger number of genes take longer to run, hence two different groups for permutation-based
chrs=("1" "2" "3" "4" "5" "6" "7" "8" "9" "17" "19")
#chrs=("10" "11" "12" "13" "14" "15" "16" "18" "20" "21" "22")

module load qtltools

cd ${wd}

#### Nominal pass
for chr in ${chrs[@]}; do

a="cphg-millerlab-vip"
p="standard"
m="8g"
t="12:00:00"

sbatch --account=$a --mem=8g -p $p -t $t --wrap="QTLtools cis --vcf ${vcfdir}/${genotypes} --bed ${phenotypes}.bed.gz --cov ${wd}/covar/${cov_file} --region chr${chr} --nominal 1 --out ../QTL_res/QTLtools/QTLtools_nominal_chr${chr}.txt --seed 29048347"

done


#### Permutation pass

a="cphg-millerlab-vip"
p="standard"
m="8g"
#t="12:00:00"
t="1-12:00:00"

for chr in ${chrs[@]}; do

sbatch --account=$a --mem=$m -p $p -t $t --wrap="QTLtools cis --vcf ${vcfdir}/${genotypes} --bed ${phenotypes}.bed.gz --cov ${wd}/covar/${cov_file} --permute 10000 --region chr${chr} --out ../QTL_res/QTLtools_UVAcoronary_July2021_perm_chr${chr}.txt --seed 29048347"

done

###### Example of running region-specific nominal pass for genes of interest
#sbatch --account=cphg-millerlab --mem=8g -p standard -t 30:00 --wrap="QTLtools cis --vcf ${vcfdir}/${genotypes} --bed ${phenotypes}.bed.gz --cov ${wd}/peer/${cov_file} --nominal 1 --out outpath/QTLtools_regioname_output.txt --region chrnum:startpos-endpos --seed 29048347"

