## Originally based off of ______.
## This file gives an overview of the process for running RFMix and preparing data for ancestry specific PCA.
## Some customization will be needed. This pipeline is for a Unix system, with the bash shell .

#!bin/bash

# Start by preparing the vcf file, which must have the admixed individuals first, 
# then the first reference, then the second reference, and so on.

###### Make a classes file that has nadmixed 0's, nref1 1's, nref2 2's, etc.  
# The numbers here are numbers of haplotypes (twice the number of individuals).
# The resulting file is space-delimited on one line.
# Before running these commands you will need to manually set the numbers of haplotypes (nadmix, nref1, etc).
nadmix="344"
nref1="314"
nref2="694"
nref3="414"
nref4="412"
for i in `seq $nadmix`; do echo "0 " >> classes; done
for i in `seq $nref1`; do echo "1 " >> classes; done
for i in `seq $nref2`; do echo "2 " >> classes; done
for i in `seq $nref3`; do echo "3 " >> classes; done
let "nref4minus = $nref4 -1"
for i in `seq $nref4minus`; do echo "4 " >> classes; done
echo "4" >> classes

### We will apply the remaining rfmix commands separately to each chromosome, then combine before running PCA.

#####  Phase the vcf file using Beagle v4.0 (this can be updated for use with Beagle v4.1 or later). 
for chr in `seq 22`; do
  java -jar beagle.r1399.jar gt=chr${chr}_unphased.vcf.gz out=chr${chr}_phased gprobs=false
done


##### Convert phased vcf to alleles file (note maintenance of allele coding)
for chr in `seq 22`; do 
  zcat chr${chr}.vcf.gz | grep -v '#' | cut -f10- | tr -d '\t|' > chr${chr}.alleles
done


##### Make the snploc file.
# For this you need a map file, which you can obtain from HapMap (be sure to obtain the correct genome build).
# The base2genetic.jar program can be obtained from the Beagle utilities page:
# https://faculty.washington.edu/browning/beagle_utilities/utilities.html
module load beagle

for chr in `seq 22`; do 
  zcat chr${chr}.vcf.gz | grep -v '^#' | cut -f2 | java -jar base2genetic.jar 1 mapfile > chr${chr}.snploc
done


##### Run RFMix (the forward-backward option generates posterior probabilities) 
dir=`pwd`
module load gcc/9.2.0 openmpi python/3.7.7
cd /home/ch2um/rfmix/rfmix-master
for chr in `seq 22`; do 
  sbatch --mem=5g -p standard --account=cphg-millerlab --wrap="/home/ch2um/rfmix/rfmix-master/rfmix PopPhased ${dir}/chr${chr}.alleles ${dir}/myfile.classes ${dir}/chr${chr}.snploc --forward-backward -o ${dir}/chr${chr} --disable-parallel -n 5"
done
cd $dir


#########################
######################### 
# NOTE: Everything below this section is part of the ____ local ancestry pipeline. RFmix results from previous step were used
# to generate regional-ancestry haploptype files for each study participant that were then used as site-specific covariates
# in local-ancestry-adjusted analyses run using R script LAmatrix.R
######################### 
#########################



##### Convert RFMix output using rephasing.
# This will replace the first $nadmix columns of the alleles file with the rephased file.
for chr in `seq 22`; do 
  cut -c `echo "${nadmix}+1" | bc`- chr${chr}.alleles > chr${chr}.temp1
  paste -d' ' chr${chr}.allelesRephased0.txt chr${chr}.temp1 | tr -d ' ' > chr${chr}.rephased
  rm chr${chr}.temp1
done


##### Apply masking, and put back the header and snp ids by using the vcf file.
for chr in `seq 22`; do 
  for ancestry in 1 2 3 4 5 6; do
    python viterbi2maskedbgl.py chr${chr}.vcf.gz chr${chr}.rephased chr${chr}.0.Viterbi.txt $ancestry | gzip > chr${chr}_anc${ancestry}.bgl.gz
  done
done


###### At this point it is a good idea to thin down the markers for the PCA (to reduce LD) - for example to 10-50 thousand markers genomewide. The program filterlines.jar can be found on the Beagle utilities website reference above
#> Here the file chr_${chr}.keepsnps has the ids of the snps to be retained. You will need to create this file.
for chr in `seq 22`; do 
  zcat chr${chr}_anc${ancestry}.bgl.gz | java -jar filterlines.jar '#' 2 chr_${chr}.keepsnps | gzip > chr${chr}_thinned_anc${ancestry}.bgl.gz
done

# Now combine the chromosomes into a single file.
for ancestry in 1 2 3 4 5 6; do
  zcat chr1_thinned_anc${ancestry}.bgl.gz | head -1 > allchr_anc${ancestry}.bgl
  for chr in `seq 22`; do
    zcat chr${chr}_thinned_anc${ancestry}.bgl.gz |  tail -n +2 >> allchr_anc${ancestry}.bgl
  done
done


##### Calculate ancestry proportions. Note that reference individuals will appear to have 100% ancestry even if not from the correct reference (since they are not masked).
# Determine which individuals have at least 50% of desired ancestry (other proportions may be substituted for 0.5).
minprop=0.5
for ancestry in 1 2 3 4 5 6; do
 cat allchr_thinned_anc${ancestry}.bgl | python ancestry_prop_minindiv.py NA > allchr_thinned_anc${ancestry}.aprop
 cat allchr_thinned_anc${ancestry}.aprop | java -jar filterlines.jar  '#' 2 ${minprop} 0.9999 | cut -d' ' -f1 > allchr_thinned_anc${ancestry}.apropgt${minprop}
done


##### Figure out which individuals are which reference using the classes file and vcf file.
##### The transpose.jar file can be obtained from the Beagle utilities website
#for ancestry in 1 2 3 4 5 6; do
#(cat allchr_thinned_anc${ancestry}.bgl | head -50 | grep "^I" | cut -d' ' -f3-; cat classes) | java -jar transpose.jar | java -jar filterlines.jar '#' 2 $ancestry $ancestry | cut -d' ' -f1 > ref$ancestry
cat ref1 ref2 ref3 > allref 


##### Filter out individuals with high levels of masking.
# Other unwanted individuals can also be filtered out at this point, including reference from other groups if desired
# Relatives should be filtered out prior to analysis: it is assumed that the file "relatives" has one id per line for the related individuals that need to be removed
# filtercolumns.jar and filterlines.jar can be obtained from the Beagle utilities website
#for ancestry in 1 2 3 4 5 6; do
# cat allchr_thinned_anc${ancestry}.apropgt${minprop} ref$ancestry | java -jar filterlines.jar '#' -1 relatives > ids_anc${ancestry}_norel
# (echo "I"; echo "id"; cat ids_anc${ancestry}) > temp${ancestry}
# cat allchr_thinned_anc${ancestry}.apropgt${minprop} allref | java -jar filterlines.jar '#' -1 relatives > ids_anc${ancestry}_allref_norel 
# (echo "I"; echo "id"; cat ids_anc${ancestry}_allref) > temp${ancestry}_allref
# cat allchr_thinned_anc${ancestry}.bgl | java -jar filtercolumns.jar '#' 1 temp${ancestry} > allchr_thinned_anc${ancestry}_apropgt${minprop}_singleref.bgl
# cat allchr_thinned_anc${ancestry}.bgl | java -jar filtercolumns.jar '#' 1 temp${ancestry}_allref > allchr_thinned_anc${ancestry}_apropgt${minprop}_allref.bgl
done
rm temp*


##### Run the PCA and produce plots.
##### The R statistical program can be obtained from https://www.r-project.org/
#R CMD BATCH mds_rcode_pipeline.R


