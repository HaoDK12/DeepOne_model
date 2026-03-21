#!/bin/bash
sample=$1
cellline=$2
ref="refs/ref_on_144.fa"

cat << EOF
#!/bin/bash
#SBATCH --account DREAM_CRISPR
#SBATCH -c 10
#SBATCH --mem 30g
#SBATCH --time 24:00:00

cd \$PBS_O_WORKDIR
## clean reads
fastp -w 10 -i rawData/${cellline}/${cellline}_on_${sample}_1.fq.gz -I rawData/${cellline}/${cellline}_on_${sample}_2.fq.gz -o cleanData/${cellline}_on_${sample}_clean_1.fq.gz -O cleanData/${cellline}_on_${sample}_clean_2.fq.gz -j cleanData/${cellline}_on_${sample}.json -h cleanData/${cellline}_on_${sample}.html
## read fastqc
fastqc -t 10 -o fastQC/ cleanData/${cellline}_on_${sample}_clean_1.fq.gz cleanData/${cellline}_on_${sample}_clean_2.fq.gz

## overlapMerge
flash -t 10 -m 10 -z -o ${cellline}_on_${sample} -d Flash/${cellline}_on_${sample} cleanData/${cellline}_on_${sample}_clean_1.fq.gz cleanData/${cellline}_on_${sample}_clean_2.fq.gz

## Align to reference
bwa mem -t 10 ${ref} Flash/${cellline}_on_${sample}/${cellline}_on_${sample}.extendedFrags.fastq.gz | samtools sort -@ 10 -o alignData/${cellline}_on_${sample}.sort.bam
##index bam
samtools index alignData/${cellline}_on_${sample}.sort.bam

## extract library
python3 scripts/get_seq.py -i alignData/${cellline}_on_${sample}.sort.bam -r ${ref} -o raw_reads/${cellline}_on_${sample} -s raw_reads/${cellline}_on_${sample}.statistic.txt -m raw_reads/${cellline}_on_${sample}.metrics.txt --leftS 0 --leftE 103 -t 37
EOF


