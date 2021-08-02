##### The purpose of this script is to call the accompanying R script to interpolate continental ancestry estimates
### from RFmix to locations in each chromosome corresponding to genetic variants, so that position-specific 
### "local ancestry" can be utilized as a covariate in eQTL analyses. Non-slurm-based systems may require different parameters.

#!bin/bash

a="account-name"
p="standard"
t="2-"
m="20g"

module load gcc/7.1.0 openmpi/3.1.4 R/4.0.3

### Submitted each chromosome separately to avoid getting put in HPC jail because these take a while to run.

for chrom in "1"; do

        while read indiv; do

        sbatch --account=$a -p $p -t $t --mem=$m --wrap="Rscript interpolate_ancestry.R $chrom $indiv"

        done < /path/to/list/of/sampleIDs/IDs_for_QTL.lst

done
