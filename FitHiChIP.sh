#!/bin/bash
#SBATCH --job-name=KO_r1 --output=KO_r1.out --error=KO_r1.err --time=10:00:00 --qos=normal --nodes=1 --mem-per-cpu=32G

#Written by Sadhana Gaddam

ml python
ml py-scipy/1.1.0_py27
ml biology
ml samtools/1.8
ml bedtools
ml htslib/1.8
ml hic-pro
ml bowtie2
ml py-macs2
ml R/3.6.1

PATH=$PATH:/home/groups/oro/software/FitHiChIP
sample=AHDC1_KO_rep1_10000

p="/oak/stanford/groups/oro/anncoll/HiChIP_NovaSeq/AHDC1_KO_HiChIP/hic_results/matrix/rep1/raw/10000/"
q="/oak/stanford/groups/oro/anncoll/HiChIP_NovaSeq/AHDC1_KO_HiChIP/hic_results/loops/AHDC1_KO_r1"
r="/oak/stanford/groups/oro/anncoll/HiChIP_NovaSeq/AHDC1_KO_HiChIP/hic_results/loops/"

awk '($2-$1)<201' ${p}/rep1_10000.matrix > ${q}/${sample}.max2mb.matrix
python /home/groups/oro/software/scripts/HiChIP/run_Mat2Bed.py -c all -k ${p}/rep1_10000_abs.bed   -m ${q}/${sample}.max2mb.matrix

ints=${q}/${sample}.max2mb.matrix.bed

cat $ints | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' > ${q}/temp_contact.bed
cat ${q}/temp_contact.bed | awk '(int(($5+$6)/2)-int(($3+$2)/2))>20000 && (int(($5+$6)/2)-int(($3+$2)/2))<2000000' > ${q}/${sample}.FitHiC_contact_counts.bed
gzip ${q}/${sample}.FitHiC_contact_counts.bed


cat ${r}/AHDC1_WT_r1/out_macs2_peaks.narrowPeak  ${r}/AHDC1_WT_r2/out_macs2_peaks.narrowPeak ${r}/AHDC1_WT_r3/out_macs2_peaks.narrowPeak > ${r}/${sample}_peaks.narrowPeak
sortBed -i ${r}/${sample}_peaks.narrowPeak > ${r}/${sample}_sort_peaks.narrowPeak
mergeBed -i ${r}/${sample}_sort_peaks.narrowPeak > ${r}/${sample}_merge.narrowPeak

bash FitHiChIP_HiCPro.sh -C ${q}/configfile_BiasCorrection_CoverageBias
