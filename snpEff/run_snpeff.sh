#!bin/bash

### The purpose of this file is to annotate lead QTLs according to their function for respective genes
### The verbose comment is necessary because regular results file doesn't include variant names to match to genes
### The output file you care about will be your slurm log file (see 'annotate_snpEff.R')
### To cut out the chaff from the log file you can grep "rs[1-9]" into a new file

snpeff="/path/to/pkg/snpEff/snpEff.jar"

in_file="/path/to/vcfs/hg38/eQTL_leadsnp_genotypes.vcf" #This has to be unzipped
db="hg38kg"

############
a="account=cphg-millerlab-vip"
p="p standard"
t="t 3:00:00"
m="mem=50g"

module load java/1.11.0

sbatch --$a --$m -$p -$t --wrap="java -jar $snpeff ann -v $db $in_file"
