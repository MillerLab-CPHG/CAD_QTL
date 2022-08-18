########### The purpose of this script is to run PAINTOR v3.0, which allows for both multi-ethnic fine-mapping using multiple reference LDs, and can utilize annotations for tissue-specific expression, etc.
##### Paintor github can be found here: https://github.com/gkichaev/PAINTOR_V3.0
##### Intro to file formats for input and requirements to run paintor available in this wiki: 
        ##### https://github.com/gkichaev/PAINTOR_V3.0/wiki/2.-Input-Files-and-Formats
##### If you use Paintor, cite this paper: https://academic.oup.com/bioinformatics/article/33/2/248/2525720

######### IMPORTANT NOTE: Note: All file formats are assumed to be single space delimited. 
######### If your file is tab-delimited you can use the following command to modify: sed -i 's/\t/ /g' <filename>

########### Running paintor ###########

for chr in {1..22}; do

module load gcc/7.1.0 openmpi/3.1.4 intel/18.0 intelmpi/18.0 mvapich2/2.3.1 scalapack/2.0.2 fftw/3.3.6 R/4.0.3 gcc/9.2.0 paintor

input_dir="/scratch/ch2um/chr${chr}plink"
output_dir="/project/cphg-millerlab/chani/QTL_res/paintor"

while read gene; do

mv $input_dir/${gene}_EUR.ld $input_dir/${gene}.ld
## Make gene-specific results files
                zcat /scratch/ch2um/Apr22_mixqtl_gwas_forpaintor_chr${chr}.txt.gz | grep "${gene}$" | expand -t 1 > $input_dir/${gene}.tmp
                zcat /scratch/ch2um/Apr22_mixqtl_gwas_forpaintor_chr${chr}.txt.gz | head -1 | expand -t 1 > $input_dir/${gene}
                        while read snp; do
                        grep "$snp\ " $input_dir/${gene}.tmp >> $input_dir/${gene}
                        done < $input_dir/new_${gene}_snplist.tmp

## Make gene-specific annotation files  
                grep "${gene}\ " ENCODE_4annot_chr${chr}.txt | cut -d ' ' -f3-8 > $input_dir/${gene}.annot.tmp
                head -1 ENCODE_4annot_chr${chr}.txt | cut -d ' ' -f5-8 > $input_dir/${gene}.annotations
                        while read snp; do
                        grep "$snp\ " $input_dir/${gene}.annot.tmp | cut -d ' ' -f3-6 >> $input_dir/${gene}.annotations
                        done < $input_dir/new_${gene}_snplist.tmp
                echo "$gene" > $input_dir/${gene}_name.tmp

        list_files="/scratch/ch2um/chr${chr}plink/${gene}_name.tmp"

        a="cphg-millerlab-vip"
        m="80g"
        t="75:00"
        p="standard"

## Example submission: PAINTOR -input.files [input filename] -in [input directory] -out [output directory] -Zhead [Zscore header(s)] -LDname [LD suffix(es)] -annotations [annotation1,annotation2...] <other options>     
## For CAD GWAS
        sbatch -o paintor_logs/${gene}_H3K27AC.log --mem=$m --account=$a -p $p -t $t --wrap="PAINTOR -input $list_files \
        -in $input_dir -out $output_dir -LDname ld,ld,ld -Zhead UVAQTL_Z,GCST011365_Z,GCST005194_Z \
        -annotations H3K27AC -Gname ${gene}_H3K27AC_enrichment -Lname ${gene}_H3K27AC_LogBayes -RESname CAD_H3K27AC_results \
        -set_seed 8675309"

        sbatch -o paintor_logs/${gene}_2annot.log --mem=$m --account=$a -p $p -t $t --wrap="PAINTOR -input $list_files \
        -in $input_dir -out $output_dir -LDname ld,ld,ld -Zhead UVAQTL_Z,GCST011365_Z,GCST005194_Z \
        -annotations H3K4ME3,H3K27AC -Gname ${gene}_H3K_enrichment -Lname ${gene}_H3K_LogBayes -RESname CAD_H3K_results \
        -set_seed 8675309"

        sbatch -o paintor_logs/${gene}_ABC.log --mem=$m --account=$a -p $p -t $t --wrap="PAINTOR -input $list_files \
        -in $input_dir -out $output_dir -LDname ld,ld,ld -Zhead UVAQTL_Z,GCST011365_Z,GCST005194_Z \
        -annotations ABC_score -Gname ${gene}_ABC_enrichment -Lname ${gene}_ABC_LogBayes -RESname CAD_ABC_results\
        -set_seed 8675309"

## For BP GWAS
        sbatch -o paintor_logs/${gene}_H3K27AC.log --mem=$m --account=$a -p $p -t $t --wrap="PAINTOR -input $list_files \
        -in $input_dir -out $output_dir -LDname ld,ld,ld,ld -Zhead UVAQTL_Z,PP_Z,SBP_Z,DBP_Z \
        -annotations H3K27AC -Gname ${gene}_H3K27AC_enrichment -Lname ${gene}_H3K27AC_LogBayes -RESname BP_H3K27AC_results \
        -set_seed 8675309"

        sbatch -o paintor_logs/${gene}_2annot.log --mem=$m --account=$a -p $p -t $t --wrap="PAINTOR -input $list_files \
        -in $input_dir -out $output_dir -LDname ld,ld,ld,ld -Zhead UVAQTL_Z,PP_Z,SBP_Z,DBP_Z \
        -annotations H3K4ME3,H3K27AC -Gname ${gene}_H3K_enrichment -Lname ${gene}_H3K_LogBayes -RESname BP_H3K_results \
        -set_seed 8675309"

        sbatch -o paintor_logs/${gene}_ABC.log --mem=$m --account=$a -p $p -t $t --wrap="PAINTOR -input $list_files \
        -in $input_dir -out $output_dir -LDname ld,ld,ld,ld -Zhead UVAQTL_Z,PP_Z,SBP_Z,DBP_Z \
        -annotations ABC_score -Gname ${gene}_ABC_enrichment -Lname ${gene}_ABC_LogBayes -RESname BP_ABC_results\
        -set_seed 8675309"

## For coronary traits GWAS
        sbatch -o paintor_logs/${gene}_H3K27AC.log --mem=$m --account=$a -p $p -t $t --wrap="PAINTOR -input $list_files \
        -in $input_dir -out $output_dir -LDname ld,ld,ld,ld,ld -Zhead UVAQTL_Z,CAC_1KG_Z,Topmed_CAC_Z,IMT_Z,Plaque_Z \
        -annotations H3K27AC -Gname ${gene}_H3K27AC_enrichment -Lname ${gene}_H3K27AC_LogBayes -RESname cor_traits_H3K27AC_results \
        -set_seed 8675309"

        sbatch -o paintor_logs/${gene}_2annot.log --mem=$m --account=$a -p $p -t $t --wrap="PAINTOR -input $list_files \
        -in $input_dir -out $output_dir -LDname ld,ld,ld,ld,ld -Zhead UVAQTL_Z,CAC_1KG_Z,Topmed_CAC_Z,IMT_Z,Plaque_Z \
        -annotations H3K4ME3,H3K27AC -Gname ${gene}_H3K_enrichment -Lname ${gene}_H3K_LogBayes -RESname cor_traits_H3K_results \
        -set_seed 8675309"

        sbatch -o paintor_logs/${gene}_ABC.log --mem=$m --account=$a -p $p -t $t --wrap="PAINTOR -input $list_files \
        -in $input_dir -out $output_dir -LDname ld,ld,ld,ld,ld -Zhead UVAQTL_Z,CAC_1KG_Z,Topmed_CAC_Z,IMT_Z,Plaque_Z \
        -annotations ABC_score -Gname ${gene}_ABC_enrichment -Lname ${gene}_ABC_LogBayes -RESname cor_traits_ABC_results\
        -set_seed 8675309"

## For cholesterol GWAS
        sbatch -o paintor_logs/${gene}_H3K27AC.log --mem=$m --account=$a -p $p -t $t --wrap="PAINTOR -input $list_files \
        -in $input_dir -out $output_dir -LDname ld,ld,ld,ld,ld -Zhead UVAQTL_Z,HDL_Z,LDL_Z,TC_Z,logTG_Z \
        -annotations H3K27AC -Gname ${gene}_H3K27AC_enrichment -Lname ${gene}_H3K27AC_LogBayes -RESname cholesterol_H3K27AC_results \
        -set_seed 8675309"

        sbatch -o paintor_logs/${gene}_2annot.log --mem=$m --account=$a -p $p -t $t --wrap="PAINTOR -input $list_files \
        -in $input_dir -out $output_dir -LDname ld,ld,ld,ld,ld -Zhead UVAQTL_Z,HDL_Z,LDL_Z,TC_Z,logTG_Z \
        -annotations H3K4ME3,H3K27AC -Gname ${gene}_H3K_enrichment -Lname ${gene}_H3K_LogBayes -RESname cholesterol_H3K_results \
        -set_seed 8675309"

        sbatch -o paintor_logs/${gene}_ABC.log --mem=$m --account=$a -p $p -t $t --wrap="PAINTOR -input $list_files \
        -in $input_dir -out $output_dir -LDname ld,ld,ld,ld,ld -Zhead UVAQTL_Z,HDL_Z,LDL_Z,TC_Z,logTG_Z \
        -annotations ABC_score -Gname ${gene}_ABC_enrichment -Lname ${gene}_ABC_LogBayes -RESname cholesterol_ABC_results\
        -set_seed 8675309"

# grep dumped paintor_logs/${gene}_H3K27AC.log >> /scratch/ch2um/Segfault_genes.lst

        done < /scratch/ch2um/paintor_genes_chr${chr}_June22.lst
done


