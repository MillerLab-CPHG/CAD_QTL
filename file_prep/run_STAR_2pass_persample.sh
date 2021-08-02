ref="/path/to/reference/annotation/file/gencode"
fastq_dir="/path/to/your/QCd/expressiondata/Fastq-Trimmed-Files"
vcf_dir="/path/to/your/vcfs"
out_dir="/path/you/want/output/files/stored/in/STAR_hg19"

module load star/2.7.2b

for i in {001..500}; do
        sample_name="SampleID${i}" #Change this to your sample ID naming scheme

    sbatch -A cphg-millerlab \
           -p parallel \
           -t 24:00:00 \
           --mem=80g \
           -N 2 \
           -n 20 \
           --wrap="STAR --genomeDir $ref/genome_index_hg19 \
                        --readFilesIn ${fastq_dir}/${sample_name}_1_val_1.fq.gz ${fastq_dir}/${sample_name}_2_val_2.fq.gz \
                        --readFilesCommand zcat \
                        --twopassMode Basic \
                        --sjdbGTFfile $ref/gencode.v37lift37.annotation.gtf \
                        --outFileNamePrefix ${out_dir}/${sample_name}_ \
                        --quantMode TranscriptomeSAM GeneCounts \
                        --outSAMtype None \
                        --varVCFfile ${vcf_dir}/${sample_name}_hg19_coronary.vcf \
                        --outSAMattributes NH HI AS NM MD vA vG \
                        --outFilterType BySJout \
                        --outFilterMultimapNmax 20 \
                        --outFilterMismatchNmax 999 \
                        --outFilterMismatchNoverReadLmax 0.04 \
                        --alignIntronMin 20 \
                        --alignIntronMax 1000000 \
                        --alignSJoverhangMin 8 \
                        --alignSJDBoverhangMin 1 \
                        --sjdbScore 1 \
                        --limitBAMsortRAM 50000000000"
done
